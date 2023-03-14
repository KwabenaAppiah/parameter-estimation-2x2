import matplotlib.pyplot as plt
import numpy as np
import datetime
import time
import math
import os


class LineGraphComp:
    def __init__(self, ev_type, pp_type):
        self._ev_type = ev_type
        self._pp_type = pp_type
        self._true_vals = []
        self._bad_matrices_exist = False
        self.set_subplots("Composite of Bad Matrices")


    #GETTERS
    # def has_bad_matrices(self):
    #     return self._mtrx_dict

    def has_bad_matrices(self):
        return self._bad_matrices_exist


    def get_ev_type(self):
        return self._ev_type

    def get_pp_type(self):
        return self._pp_type

    def get_subplots(self):
        return self._fig, self._ax

    def get_true_vals(self):
        return self._true_vals


#SETTERS
    def set_has_bad_matrices(self, bool):
        self._bad_matrices_exist = bool

    def set_subplots(self, name):
        true_vals, ev_type, pp_type = self.get_true_vals(), self.get_ev_type(), self.get_pp_type()
        self._fig, self._ax = plt.subplots()
        title = ev_type.upper() + " | " + pp_type.upper() + " - Matrices w/ an Avg. Err. > 1e-3"
        self._ax.set_yscale('log')
        self._ax.set_title(label = title, pad = 25)
        self._ax.set_xlabel("Iterations")
        self._ax.set_ylabel("Avg. Rel. Errors")
        self._fig.set_size_inches(16, 8)


    def organize_data(self, guesses, true_vals, cycle_num):
        a11s, a12s, a21s, a22s = guesses
        a11, a12, a21, a22  = true_vals
        total_err_steps = []
        a21_rel_err = []
        a22_rel_err = []

        if a21 != 0 and a22 != 0:
            j = 0
            for i in range(len(a21s)):
                a21_rel_err.append(abs(a21s[i] - a21) / abs(a21))
                a22_rel_err.append(abs(a22s[i] - a22) / abs(a22))
                total_err_steps.append( abs(a21_rel_err[i] + a22_rel_err[i]) / 2)

            #Select based on final val (i.e. like Tr / Det plane)
            a21_nth_rel_err = abs(a21s[-1] - a21) / abs(a21)
            a22_nth_rel_err = abs(a22s[-1] - a22) / abs(a22)
            a21_a22_nth_avg_rel_err = abs(a21_nth_rel_err +  a22_nth_rel_err) / 2

            if a21_a22_nth_avg_rel_err  > 1e-3:
                self.plot_graph(cycle_num, total_err_steps)
                self.set_has_bad_matrices(True)
                return True

        else:
            # print("a21:", abs(a21), "and a22:", abs(a22) )
            # print("Mtrx.", true_vals, "is not plottable.")
            return False


    def plot_graph(self, cycle_num, total_err_steps):
        fig, ax = self.get_subplots()
        ax.plot(list(range(len(total_err_steps))), total_err_steps, label = "MX" + " " + str(cycle_num))

    def get_date_str(self):
        raw_date = datetime.datetime.now()
        t = time.localtime()
        # current_time = time.strftime("%H.%M", t)
        # date_str = str(raw_date.year) + "." + str(time.strftime('%m')) + "." + str(time.strftime('%d')) + current_time
        date_str = str(raw_date.year) + "." + str(time.strftime('%m')) + "." + str(time.strftime('%d'))
        return date_str


    def display(self):
        if self.has_bad_matrices() != False:
            true_vals, ev_type, pp_type = self.get_true_vals(), self.get_ev_type(), self.get_pp_type()

            fig, ax = self.get_subplots()
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
            fig, ax = self.get_subplots()
            plt.close(fig)
            print("RESULT: No bad matrices were found.")
