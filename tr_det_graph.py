import datetime
import math
import os
import time
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np


class TrDetGraph:
    def __init__(self, ev_type, pp_type, loop_limit, mu_val):
        self._ev_type = ev_type
        self._pp_type = pp_type
        self._max_trace = 0
        self.set_subplots(loop_limit, mu_val)

    def get_ev_type(self):
        return self._ev_type

    def get_pp_type(self):
        return self._pp_type

    def get_subplots(self):
        return self._fig, self._ax

    def get_max_trace(self):
        return self._max_trace

    def set_max_trace(self, new_max_trace):
        new_max_trace = abs(new_max_trace)

        if new_max_trace > self._max_trace:
            self._max_trace = new_max_trace

        # if new_max_trace > self._max_trace:
        #  print("Old Max Trace:", self._max_trace)
        #  self._max_trace = new_max_trace
        #  print("New Max Trace:", self._max_trace)
        # else:
        #  print("Max trace is unchanged:", self._max_trace)

    def set_subplots(self, loop_limit, mu_val):
        self._fig, self._ax = plt.subplots()
        # cycles_abrev = "cc"
        # mu_formatted = str("{:.2e}".format(mu_val))
        # title = self.get_ev_type().upper() + " | " + self.get_pp_type().upper() + " - " + str(loop_limit) + " cc | µ: " + mu_formatted
        title = self.get_ev_type().upper() + " | " + self.get_pp_type().upper() + " - " + str(loop_limit) + "cc"
        self._ax.set_title(label = title, pad = 20)
        self._ax.set_xlabel("Tr", loc = "right")
        self._ax.set_ylabel("Det", loc = "top")
        self._fig.set_size_inches(16, 8)
        self._ax = plt.gca()

        # Hide two spines
        self._ax.spines["right"].set_color("none")
        self._ax.spines["top"].set_color("none")

        # Move bottom and left spine to 0, 0
        self._ax.spines["bottom"].set_position(("data", 0))
        self._ax.spines["left"].set_position(("data", 0))

        # Move ticks positions
        self._ax.xaxis.set_ticks_position("bottom")
        self._ax.yaxis.set_ticks_position("left")

        self._ax.plot(1, 0, ">k", transform = self._ax.get_yaxis_transform(), clip_on = False)
        self._ax.plot(0, 1, "^k", transform = self._ax.get_xaxis_transform(), clip_on = False)

    def organize_data(self, guesses, true_vals):
        a11, a12, a21, a22 = true_vals
        a11s, a12s, a21s, a22s = guesses
        a21_rel_err = abs(a21s[-1] - a21) / abs(a21)
        a22_rel_err = abs(a22s[-1] - a22) / abs(a22)
        a21_a22_avg_rel_err = abs(a21_rel_err + a22_rel_err) / 2
        data_was_plotted = False
        # a21_str = 'a\N{SUBSCRIPT TWO}\N{SUBSCRIPT ONE}'
        # a22_str = 'a\N{SUBSCRIPT TWO}\N{SUBSCRIPT TWO}'
        a21_str = "a21"
        a22_str = "a22"

        if a21 != 0 and a22 != 0 and math.isinf(a21_a22_avg_rel_err) != True:
            trace = a11 + a22
            det = a11 * a22 - a12 * a21
            print("")
            print(np.matrix([[a11, a12], [a21, a22]]), "\n")
            print("Tr:", trace, "", "Det:", det, "\n")
            print(a21_str, "and", a22_str, "Avg. rel. err:", a21_a22_avg_rel_err, "\n")
            print(a21_str, "Avg. rel. err:", a21_rel_err)
            print(a22_str, "Avg. rel. err:", a22_rel_err, "\n")
            self.plot_points(trace, det, a21_a22_avg_rel_err)
            self.set_max_trace(trace)
            data_was_plotted = True

        else:
            print(a21_str + ":", a21, a22_str + ":", a22)
            print("This matrix is not plottable.", "\n")
            data_was_plotted = False

        return data_was_plotted

    def plot_points(self, x, y, avg_err):
        fig, ax = self.get_subplots()

        if avg_err > 1e-1:
            ax.plot(x, y, ".", color="red", alpha=1, zorder=10, markersize=10)
        elif avg_err > 1e-3:
            ax.plot(x, y, ".", color="red", alpha=0.3, zorder=10, markersize=10)
        elif avg_err > 1e-5:
            ax.plot(x, y, ".", color="#800080", alpha=0.3, zorder=10, markersize=10)
        elif avg_err > 1e-8:
            ax.plot(x, y, ".", color="blue", alpha=0.3, zorder=10, markersize=10)
        else:
            ax.plot(x, y, ".", color="blue", alpha=1, zorder=10, markersize=10)

    def graph_parabola(self):
        fig, ax = self.get_subplots()
        max_trace = self.get_max_trace()

        if abs(max_trace) == 0:
            max_trace = 10

        x = np.linspace(-(max_trace), max_trace, 100)
        y = x**2 / 4
        ax.plot(x, y, linewidth=1.5, c="orange")

    def get_date_str(self):
        raw_date = datetime.datetime.now()
        t = time.localtime()
        # current_time = time.strftime("%H.%M", t)
        # date_str = str(raw_date.year) + "." + str(time.strftime('%m')) + "." + str(time.strftime('%d')) + "_" + str(current_time)
        date_str = str(raw_date.year) + "." + str(time.strftime("%m")) + "." + str(time.strftime("%d"))
        return date_str

    def display(self, ev_type, pp_type, loop_limit):
        fig, ax = self.get_subplots()
        self.graph_parabola()

        custom_handles = [
            Line2D([0],[0], marker="o", markerfacecolor="r", color="w", alpha=1, markersize=7, label="$x̄_{Err}$ > 1e-1"),
            Line2D([0],[0], marker="o", markerfacecolor="r", color="w", alpha=0.3, markersize=7, label="$x̄_{Err}$ > 1e-3"),
            Line2D([0], [0], marker="o", markerfacecolor="#800080", color="w", alpha=0.3, markersize=7, label="$x̄_{Err}$ > 1e-5"),
            Line2D([0],[0], marker="o", markerfacecolor="b", color="w", alpha=0.3, markersize=7, label="$x̄_{Err}$ > 1e-8"),
            Line2D([0],[0], marker="o", markerfacecolor="b", color="w", markersize=7, label="Below 1e-8"),
            Line2D([0], [0], color="orange", alpha=1, lw=3, label="T\N{SUPERSCRIPT TWO} - 4D = 0")]

        ax.legend(handles=custom_handles, loc="best", bbox_to_anchor=(1, 0.5))
        ev_pp_type = self.get_ev_type() + "_" + self.get_pp_type()
        filename = "tr_det_graph" + "_" + ev_pp_type + "_" + str(loop_limit) + "_" + "cc"
        subdir = "../output/" + ev_pp_type + "_" + self.get_date_str() + "/"
        os.makedirs(subdir, exist_ok=True)
        fig.savefig(subdir + filename + ".png", dpi=300)
        plt.close(fig)
