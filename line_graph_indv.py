import datetime
import math
import os
import time
import matplotlib.pyplot as plt
import numpy as np


class LineGraphIndv:
    def __init__(self, ev_type, pp_type):
        self._ev_type = ev_type
        self._pp_type = pp_type

    def get_ev_type(self):
        return self._ev_type

    def get_pp_type(self):
        return self._pp_type

    def plot_graph(self, guesses, true_vals, cycle_num):
        a11s, a12s, a21s, a22s = guesses
        a11, a12, a21, a22 = true_vals
        total_err_steps = []
        a21_rel_err = []
        a22_rel_err = []

        if a21 != 0 and a22 != 0:
            for i in range(len(a21s)):
                a21_rel_err.append(abs(a21s[i] - a21) / abs(a21))
                a22_rel_err.append(abs(a22s[i] - a22) / abs(a22))
                total_err_steps.append(abs(a21_rel_err[i] + a22_rel_err[i]) / 2)

            # Select color based on final val (i.e. like Tr / Det plane)
            a21_nth_rel_err = abs(a21s[-1] - a21) / abs(a21)
            a22_nth_rel_err = abs(a22s[-1] - a22) / abs(a22)
            a21_a22_nth_avg_rel_err = abs(a21_nth_rel_err + a22_nth_rel_err) / 2

            if a21_a22_nth_avg_rel_err > 1e-3:
                self.display(true_vals, cycle_num, total_err_steps)
                self.write_to_file("Mtrx. " + str(cycle_num) + ": " + str(true_vals), "Avg. Err.: " + str(a21_a22_nth_avg_rel_err))
                # self.set_bad_matrices("MTRX: " + str(cycle_num) + str(true_vals) +"\n" + "AVG. ERR: " + str(a21_a22_nth_avg_rel_err))
                # self.plot_fig_comp(cycle_num, total_err_steps)
                return True

            else:
                # print("a21:", abs(a21), "and a22:", abs(a22) )
                # print("Mtrx.", true_vals, "is not plottable.")
                return False

    def write_to_file(self, str_1, str_2):
        output = str_1 + "\n" + str_2 + "\n" + "\n"
        ev_pp_type = self.get_ev_type() + "_" + self.get_pp_type()
        filename = "bad_matrices" + "_" + ev_pp_type + ".txt"
        subdir = "../output/" + ev_pp_type + "_" + self.get_date_str() + "/"
        os.makedirs(subdir, exist_ok=True)
        f = open(subdir + filename, "a")
        f.write(output)
        f.close()

    def get_date_str(self):
        raw_date = datetime.datetime.now()
        t = time.localtime()
        # current_time = time.strftime("%H.%M", t)
        # date_str = str(raw_date.year) + "." + str(time.strftime('%m')) + "." + str(time.strftime('%d')) + current_time
        date_str = str(raw_date.year) + "."+ str(time.strftime("%m"))+ "." + str(time.strftime("%d"))
        return date_str

    def display(self, true_vals, cycle_num, total_err_steps):
        ev_type = self.get_ev_type()
        pp_type = self.get_pp_type()
        fig, ax = plt.subplots()
        fig.set_size_inches(7, 5)
        ax.plot(list(range(len(total_err_steps))), total_err_steps, label="MX" + " " + str(cycle_num))
        title = ev_type.upper()+ " | " + pp_type.upper() + " - MATRIX " + str(cycle_num) + " : " + str(true_vals)

        ax.set_title(label = title, pad = 20)
        ax.set_xlabel("Iterations")
        ax.set_ylabel("Avg. Rel. Errors")
        ax.set_yscale("log")

        # Output IMG files
        output_type = "indv"
        ev_pp_type = self.get_ev_type() + "_" + self.get_pp_type()
        subdir = "../output/" + ev_pp_type + "_" + self.get_date_str() + "/line_graph" + "_" + output_type+ "/"
        filename = "line_graph" + "_" + ev_pp_type
        os.makedirs(subdir, exist_ok=True)
        fig.savefig(subdir + filename + "_mtrx_" + str(cycle_num) + ".png", dpi=300)
        plt.close(fig)
