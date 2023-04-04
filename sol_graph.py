import datetime
import math
import os
import time
import matplotlib.pyplot as plt
import numpy as np


class SolGraph:
    def __init__(self, ev_type, pp_type):
        # FOR INDV AND COMP
        self._ev_type = ev_type
        self._pp_type = pp_type

    # INDV AND COMP GETTERS
    def get_ev_type(self):
        return self._ev_type

    def get_pp_type(self):
        return self._pp_type

    def get_date_str(self):
        raw_date = datetime.datetime.now()
        t = time.localtime()
        date_str = str(raw_date.year) + "."+ str(time.strftime("%m"))+ "." + str(time.strftime("%d"))
        return date_str

    # def get_solutions

    def init(self, sol, guesses, true_vals, cycle_num, set_boundry_val):
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

            if a21_a22_nth_avg_rel_err >= set_boundry_val: # bad matricies
                self.display_sol(sol, true_vals, cycle_num, "bad")
                return True

            elif a21_a22_nth_avg_rel_err < set_boundry_val: # good matricies
                self.display_sol(sol, true_vals, cycle_num, "good")
                return True

            else:
                return False


    def format_fl_vals(self, val):
        return "{:.3f}".format(val)

    def display_sol(self, sol, true_vals, cycle_num, sol_type):
        t = np.array(sol.t)
        y0, y1, y2, y3 = sol.y[0], sol.y[1], sol.y[2], sol.y[3]

        ev_type = self.get_ev_type()
        pp_type = self.get_pp_type()
        fig, ax = plt.subplots()
        fig.set_size_inches(16, 8)

        ax.plot(t, y0, label = "Sol 1" )
        ax.plot(t, y1, label = "Sol 2" )
        ax.plot(t, y2, label = "Sol 3" )
        ax.plot(t, y3, label = "Sol 4" )

        # For file title
        a11, a12, a21, a22 = true_vals
        a11_str, a12_str = self.format_fl_vals(a11), self.format_fl_vals(a12)
        a21_str, a22_str = self.format_fl_vals(a21), self.format_fl_vals(a22)
        true_vals_str = "[" + a11_str +", " + a12_str + ", "+ a21_str +", " + a22_str + "]"
        title = ev_type.upper()+ " | " + pp_type.upper() + " - MTRX " + str(cycle_num) + " : " + true_vals_str
        ax.set_title(label = title, pad = 20)

        ax.set_xlabel("T", loc = "right")
        ax.set_ylabel("Y", loc = "top")
        ax.set_xscale("log")
        ax.set_yscale("log")

        # Output IMG files
        ev_pp_type = self.get_ev_type() + "_" + self.get_pp_type()
        subdir = "../output/" + ev_pp_type + "_" + self.get_date_str() + "/sol_graph" + "_" + sol_type + "/"
        filename = "sol_graph" + "_" + ev_pp_type
        ax.legend(loc="best", bbox_to_anchor=(1, 0.5))
        os.makedirs(subdir, exist_ok = True)
        fig.savefig(subdir + filename + "_mtrx_" + str(cycle_num) + ".png", dpi = 300)
        plt.close(fig)
