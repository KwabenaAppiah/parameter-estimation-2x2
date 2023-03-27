import datetime
import math
import os
import time
import matplotlib.pyplot as plt
import numpy as np


class LineGraph:
    def __init__(self, ev_type, pp_type):
        # FOR INDV AND COMP
        self._ev_type = ev_type
        self._pp_type = pp_type

        # FOR COMP ONLY
        self._true_vals = []
        self._bad_matrices_exist = False
        self._is_there_a_plottable_mtrx = False
        self.set_subplots_comp()

    # INDV AND COMP GETTERS
    def get_ev_type(self):
        return self._ev_type

    def get_pp_type(self):
        return self._pp_type

    def get_date_str(self):
        raw_date = datetime.datetime.now()
        t = time.localtime()
        # current_time = time.strftime("%H.%M", t)
        # date_str = str(raw_date.year) + "." + str(time.strftime('%m')) + "." + str(time.strftime('%d')) + current_time
        date_str = str(raw_date.year) + "."+ str(time.strftime("%m"))+ "." + str(time.strftime("%d"))
        return date_str

    # def is_there_a_plottable_mtrx(self):
    #     return self._is_there_a_plottable_mtrx

    # COMP GETTERS ONLY
    def get_subplots_comp(self):
        return self._fig, self._ax

    def get_true_vals_comp(self):
        return self._true_vals

    def has_bad_matrices_comp(self):
        return self._bad_matrices_exist


    # COMP SETTERS
    def set_has_bad_matrices_comp(self, bool):
        self._bad_matrices_exist = bool


    # def set_subplots_comp(self, name):
    def set_subplots_comp(self):
        true_vals, ev_type, pp_type = self.get_true_vals_comp(), self.get_ev_type(), self.get_pp_type()
        self._fig, self._ax = plt.subplots()
        title = ev_type.upper() + " | " + pp_type.upper() + " - Matrices w/ an Avg. Err. > 1e-3"
        self._ax.set_yscale('log')
        self._ax.set_title(label = title, pad = 25)
        self._ax.set_xlabel("Iterations")
        self._ax.set_ylabel("Avg. Rel. Errors")
        self._fig.set_size_inches(16, 8)

    # def set_is_there_a_plottable_mtrx(self, bool):
    #     self._is_there_a_plottable_mtrx = bool


    def init(self, guesses, true_vals, cycle_num, set_boundry_val):
        a11s, a12s, a21s, a22s = guesses
        a11, a12, a21, a22 = true_vals
        total_err_steps = []
        a21_rel_err = []
        a22_rel_err = []
        #is_there_a_plottable_mtrx = False

        if a21 != 0 and a22 != 0:
            for i in range(len(a21s)):
                a21_rel_err.append(abs(a21s[i] - a21) / abs(a21))
                a22_rel_err.append(abs(a22s[i] - a22) / abs(a22))
                total_err_steps.append(abs(a21_rel_err[i] + a22_rel_err[i]) / 2)

            # Select color based on final val (i.e. like Tr / Det plane)
            a21_nth_rel_err = abs(a21s[-1] - a21) / abs(a21)
            a22_nth_rel_err = abs(a22s[-1] - a22) / abs(a22)
            a21_a22_nth_avg_rel_err = abs(a21_nth_rel_err + a22_nth_rel_err) / 2


            # if self.is_there_a_plottable_mtrx() == True:  # all matricies
            #     print("YES!")
            #     self.write_to_file("mtrx_" + str(cycle_num) + "|" + str(a21_a22_nth_avg_rel_err) + "|" + str(true_vals), "all")
            #     return True


            if a21_a22_nth_avg_rel_err >= set_boundry_val: # bad matricies
                self.display_indv(true_vals, cycle_num, total_err_steps)
                self.write_to_file("mtrx_" + str(cycle_num) + "|" + str(a21_a22_nth_avg_rel_err) + "|" + str(true_vals), "bad")
                self.write_to_file("mtrx_" + str(cycle_num) + "|" + str(a21_a22_nth_avg_rel_err) + "|" + str(true_vals), "all")
                self.plot_graph_comp(cycle_num, total_err_steps)
                self.set_has_bad_matrices_comp(True)
                # self.set_is_there_a_plottable_mtrx(True)
                return True

            elif a21_a22_nth_avg_rel_err < set_boundry_val: # good matricies
                self.write_to_file("mtrx_" + str(cycle_num) + "|" + str(a21_a22_nth_avg_rel_err) + "|" + str(true_vals), "good")
                self.write_to_file("mtrx_" + str(cycle_num) + "|" + str(a21_a22_nth_avg_rel_err) + "|" + str(true_vals), "all")
                # self.set_is_there_a_plottable_mtrx(True)
                return True

            else:
                # print("a21:", abs(a21), "and a22:", abs(a22) )
                # print("Mtrx.", true_vals, "is not plottable.")
                return False


    def plot_graph_comp(self, cycle_num, total_err_steps):
        fig, ax = self.get_subplots_comp()
        ax.plot(list(range(len(total_err_steps))), total_err_steps, label = "MX" + " " + str(cycle_num))


    def write_to_file(self, line_str, mtrx_file_type):
        output = line_str + "\n" + "\n"
        ev_pp_type = self.get_ev_type() + "_" + self.get_pp_type()
        filename = ev_pp_type + "_" + mtrx_file_type + "_matrices.txt"
        subdir = "../output/" + ev_pp_type + "_" + self.get_date_str() + "/txt_files/"
        os.makedirs(subdir, exist_ok = True)
        f = open(subdir + filename, "a")
        f.write(output)
        f.close()

    def format_fl_vals(self, val):
        return "{:.3f}".format(val)


    def display_indv(self, true_vals, cycle_num, total_err_steps):
        ev_type = self.get_ev_type()
        pp_type = self.get_pp_type()
        fig, ax = plt.subplots()
        fig.set_size_inches(7, 5)
        ax.plot(list(range(len(total_err_steps))), total_err_steps, label = "MX" + " " + str(cycle_num))
        a11, a12, a21, a22 = true_vals
        a11_str, a12_str = self.format_fl_vals(a11), self.format_fl_vals(a12)
        a21_str, a22_str = self.format_fl_vals(a21), self.format_fl_vals(a22)
        true_vals_str = "[" + a11_str +", " + a12_str + ", "+ a21_str +", " + a22_str + "]"
        # title = ev_type.upper()+ " | " + pp_type.upper() + " - MATRIX " + str(cycle_num) + " : " + str(true_vals)
        title = ev_type.upper()+ " | " + pp_type.upper() + " - MTRX " + str(cycle_num) + " : " + true_vals_str
        ax.set_title(label = title, pad = 20)
        ax.set_xlabel("Iterations")
        ax.set_ylabel("Avg. Rel. Errors")
        ax.set_yscale("log")

        # Output IMG files
        output_type = "indv"
        ev_pp_type = self.get_ev_type() + "_" + self.get_pp_type()
        subdir = "../output/" + ev_pp_type + "_" + self.get_date_str() + "/line_graph" + "_" + output_type+ "/"
        filename = "line_graph" + "_" + ev_pp_type
        os.makedirs(subdir, exist_ok = True)
        fig.savefig(subdir + filename + "_mtrx_" + str(cycle_num) + ".png", dpi = 300)
        plt.close(fig)


    def display_comp(self):
        if self.has_bad_matrices_comp() != False:
            true_vals, ev_type, pp_type = self.get_true_vals_comp(), self.get_ev_type(), self.get_pp_type()
            fig, ax = self.get_subplots_comp()
            ax.legend(loc = 'center left', bbox_to_anchor = (1, 0.5))

            #Output files
            output_type = "comp"
            ev_pp_type = self.get_ev_type() + '_' + self.get_pp_type()
            # sub_dir =   "../output/line_graph" + "_" + output_type + "_" + ev_type + "_" + pp_type
            subdir =   "../output/" + ev_pp_type + "_" + self.get_date_str() + "/"
            filename =  "line_graph" + "_" + ev_pp_type + "_" + output_type
            os.makedirs(subdir, exist_ok = True)
            fig.savefig(subdir + filename + ".png", dpi = 300)
            plt.close(fig)

        else:
            fig, ax = self.get_subplots_comp()
            plt.close(fig)
            print("RESULT: No bad matrices were found.")
