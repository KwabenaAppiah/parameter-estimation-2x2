from matplotlib.lines import Line2D
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import numpy as np
import datetime
import math
import time
import os


class EigenvaluesGraph:
    def __init__(self, ev_type, pp_type, case_type, bounds, loop_limit ):
        self._static_vars_dict = {
            'ev_type': ev_type, 'pp_type': pp_type, 'bounds': str(bounds), 'loop_limit': loop_limit, "case_type": case_type, "param_label_1": "", "param_label_2": ""
        }
        self.set_param_labels()
        self._err_count_dict = {'under_1e-8': 0, 'over_1e-8': 0, 'over_1e-5': 0, 'over_1e-3': 0,  'over_1e-1': 0 }
        self.set_subplots(ev_type, pp_type, case_type, str(bounds), str(loop_limit))



    def get_err_count_dict_elt(self, key):
        return self._err_count_dict[key]

    def set_err_count_dict_elt(self, key, value):
        self._err_count_dict[key] = value

    def get_subplots(self):
        return self._fig, self._ax

    def set_param_labels(self):
        if self.get_static_vars_dict_elt("case_type") == "main_diagonal":
            self.set_static_vars_dict_elt("param_label_1", "$a_{11}$")
            self.set_static_vars_dict_elt("param_label_2", "$a_{22}$")

        elif self.get_static_vars_dict_elt("case_type") == "anti-diagonal":
            self.set_static_vars_dict_elt("param_label_1", "$a_{12}$")
            self.set_static_vars_dict_elt("param_label_2", "$a_{21}$")

        elif self.get_static_vars_dict_elt("case_type") == "left_column":
            self.set_static_vars_dict_elt("param_label_1", "$a_{11}$")
            self.set_static_vars_dict_elt("param_label_2", "$a_{21}$")

        elif self.get_static_vars_dict_elt("case_type") == "right_column":
            self.set_static_vars_dict_elt("param_label_1", "$a_{12}$")
            self.set_static_vars_dict_elt("param_label_2", "$a_{22}$")


    def set_subplots(self, ev_type, pp_type, case_type, bounds, loop_limit):
        self._fig, self._ax = plt.subplots()
        param_label_1, param_label_2 = self.get_static_vars_dict_elt("param_label_1"), self.get_static_vars_dict_elt("param_label_2")
        graph_description = "Avg. Relative Error of " + param_label_1 + " and " + param_label_2
        title = ev_type.upper() + " | " + pp_type.upper() + " | " + case_type + " | BNDS " + bounds + " | " + loop_limit + " cc | " + graph_description
        self._ax.set_title(label = title, pad = 30, fontsize = 15)
        self._ax.set_xlabel("$\u03BB_{1}$", loc = "right", fontsize = 14)
        self._ax.set_ylabel("$\u03BB_{2}$", loc = "top", fontsize = 14)
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

    def get_eignevalues(self, A):
        eigenvalues, eigenvectors = np.linalg.eig(A)
        ev_1, ev_2 = eigenvalues
        ev_type = self.get_static_vars_dict_elt("ev_type")

        # if(ev_type == "ce"):
        # if isinstance(ev_1, complex) and isinstance(ev_2, complex):
        #     ev_1 = math.sqrt((ev_1.real) ** 2 + (ev_1.imag) ** 2 )
        #     ev_2 = math.sqrt((ev_2.real) ** 2 + (ev_2.imag) ** 2 )
        if isinstance(ev_1, complex):
            ev_1 = math.sqrt((ev_1.real) ** 2 + (ev_1.imag) ** 2 )

        if isinstance(ev_2, complex):
            ev_2 = math.sqrt((ev_2.real) ** 2 + (ev_2.imag) ** 2 )


        return ev_1, ev_2

    # def get_avg_list(self, list_1, list_2):
    #     avg_list = []
    #     for i in range(len(list_1)):
    #         avg_list.append(list_1[i] + list_2[i] / 2)
    #     return avg_list


    def organize_data(self, A, param_estimates, avg_param_errors, case_type):

        param_labels = (self.get_static_vars_dict_elt("param_label_1"), self.get_static_vars_dict_elt("param_label_2"))

        # Parameter stuff
        a11, a12, a21, a22     = A[0,0], A[0,1], A[1,0], A[1,1]
        param_1_estimates, param_2_estimates = [], []
        param_1_estimates, param_2_estimates = param_estimates
        param_label_1, param_label_2 = (self.get_static_vars_dict_elt("param_label_1"), self.get_static_vars_dict_elt("param_label_2"))
        avg_abs_param_err, avg_rel_param_err =  avg_param_errors
        nth_avg_rel_param_err = avg_rel_param_err[-1]
        data_was_plotted = False
        param_1 = param_2 = 0

        if case_type == "main_diagonal":
            param_1, param_2 = a11, a22

        elif(case_type == "anti-diagonal"):
            param_1, param_2 = a12, a21

        elif(case_type == "left_column"):
            param_1, param_2 = a11, a21

        elif(case_type == "right_column"):
            param_1, param_2 = a12, a22



        if param_1 != 0 and param_2 != 0 and math.isinf(nth_avg_rel_param_err) != True:
            ev_1, ev_2 = self.get_eignevalues(A)
            self.plot_points(ev_1, ev_2, nth_avg_rel_param_err)
            data_was_plotted = True

        else:
            print(param_label_1 + ":", param_1, param_label_2 + ":", param_2)
            data_was_plotted = False

        return data_was_plotted


    def plot_points(self, x, y, nth_avg_rel_param_err):
        fig, ax = self.get_subplots()

        if 1e-1 < nth_avg_rel_param_err and nth_avg_rel_param_err < math.inf:
            ax.plot(x, y, ".", color = "red", alpha = 1, zorder = 10, markersize = 10)
            self.set_err_count_dict_elt("over_1e-1", self.get_err_count_dict_elt("over_1e-1") + 1)

        elif 1e-3 < nth_avg_rel_param_err and nth_avg_rel_param_err <= 1e-1:
            ax.plot(x, y, ".", color = "red", alpha = 0.3, zorder = 10, markersize = 10)
            self.set_err_count_dict_elt("over_1e-3", self.get_err_count_dict_elt("over_1e-3") + 1)

        elif 1e-5 < nth_avg_rel_param_err and nth_avg_rel_param_err <= 1e-3:
            ax.plot(x, y, ".", color = "#800080", alpha = 0.3, zorder = 10, markersize = 10)
            self.set_err_count_dict_elt("over_1e-5", self.get_err_count_dict_elt("over_1e-5") + 1)

        elif 1e-8 < nth_avg_rel_param_err and nth_avg_rel_param_err <= 1e-5:
            ax.plot(x, y, ".", color = "blue", alpha = 0.3, zorder = 10, markersize = 10)
            self.set_err_count_dict_elt("over_1e-8", self.get_err_count_dict_elt("over_1e-8") + 1)

        elif 0 <= nth_avg_rel_param_err and nth_avg_rel_param_err <= 1e-8:
            ax.plot(x, y, ".", color = "blue", alpha = 1, zorder = 10, markersize = 10)
            self.set_err_count_dict_elt("under_1e-8", self.get_err_count_dict_elt("under_1e-8") + 1)


    def get_date_str(self):
        raw_date = datetime.datetime.now()
        t = time.localtime()
        date_str = str(raw_date.year) + "." + str(time.strftime("%m")) + "." + str(time.strftime("%d"))
        return date_str

    def display(self, ev_type, pp_type, case_type, loop_limit ):
        loop_limit = str(loop_limit)
        self.display_ev_graph(ev_type, pp_type, case_type, loop_limit )



    def display_ev_graph(self, ev_type, pp_type, case_type, loop_limit):
        fig, ax = self.get_subplots()
        bbox_x = bbox_y = 0
        legend_loc = "center left"

        custom_handles = [
            Line2D([0], [0], marker = "o", markerfacecolor = "r", color = "w", alpha = 1, markersize = 7, label = "1e-1 < $x̄_{Err}$ < ∞"),
            Line2D([0], [0], marker = "o", markerfacecolor = "r", color = "w", alpha = 0.3, markersize = 7, label = "1e-3 < $x̄_{Err}$ <= 1e-1"),
            Line2D([0], [0], marker = "o", markerfacecolor = "#800080", color = "w", alpha = 0.3, markersize = 7, label = "1e-5 < $x̄_{Err}$ <= 1e-3"),
            Line2D([0], [0], marker = "o", markerfacecolor = "b", color = "w", alpha = 0.3, markersize = 7, label = "1e-8 < $x̄_{Err}$ <= 1e-5"),
            Line2D([0], [0], marker = "o", markerfacecolor = "b", color = "w", markersize = 7, label = "0 <= $x̄_{Err}$ <= 1e-8")]

        if(ev_type == "rde" and pp_type == "saddle"):
            legend_loc, bbox_x, bbox_y = "upper left", -0.15, 1

        elif(ev_type == "rde" and pp_type == "sink"):
            legend_loc, bbox_x, bbox_y = "upper left", -0.16, .9

        elif(ev_type == "rde" and pp_type == "source"):
            legend_loc, bbox_x, bbox_y = "upper left", -0.16, .9

        elif(ev_type == "re" and pp_type == "sink" ):
            legend_loc, bbox_x, bbox_y = "lower right", .91, -.06

        elif(ev_type == "re" and pp_type == "source" ):
            legend_loc, bbox_x, bbox_y = "lower right", 1, .1

        elif(ev_type == "ce"):
            legend_loc, bbox_x, bbox_y = "lower right", 1, .1


        ax.legend(handles = custom_handles, loc = legend_loc, bbox_to_anchor = (bbox_x, bbox_y), borderpad = 1.1)
        ev_and_pp_type = ev_type + "_" + pp_type
        filename = "ev_graph" + "_" + ev_and_pp_type + "_" + loop_limit + "_" + "cc"
        subdir = "../output/" + ev_and_pp_type + "_" + case_type + "_" + self.get_date_str() + "/"
        os.makedirs(subdir, exist_ok = True)
        fig.savefig(subdir + filename + ".png", dpi = 300)
        plt.close(fig)



    def get_static_vars_dict(self):
        return self._static_vars_dict

    def get_static_vars_dict_elt(self, key):
        return self._static_vars_dict[key]

    def set_static_vars_dict_elt(self, key, value):
        self._static_vars_dict[key] = value
