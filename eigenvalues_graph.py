from matplotlib.lines import Line2D
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import numpy as np
import datetime
import math
import time
import os


class EigenvaluesGraph:
    def __init__(self, *args):
        ev_type, pp_type, integration_method, bounds, loop_limit, case_type = args
        self._sim_vars_dict = {
            'ev_type': ev_type, 'pp_type': pp_type, 'integration_method': integration_method, 'bounds': bounds, 'loop_limit': loop_limit, 'case_type': case_type
        }
        self._err_count_dict = {'under_1e-8': 0, 'over_1e-8': 0, 'over_1e-5': 0, 'over_1e-3': 0,  'over_1e-1': 0 }
        self.set_subplots(ev_type, pp_type, integration_method, bounds, loop_limit, case_type)


    def get_err_count_dict_elt(self, key):
        return self._err_count_dict[key]

    def set_err_count_dict_elt(self, key, value):
        self._err_count_dict[key] = value

    def get_sim_vars_dict(self):
        return self._sim_vars_dict

    def get_sim_vars_dict_elt(self, key):
        return self._sim_vars_dict[key]

    def get_subplots(self):
        return self._fig, self._ax


    def set_subplots(self, *args):
        self._fig, self._ax = plt.subplots()
        ev_type, pp_type, integration_method, bounds, loop_limit, case_type = args
        loop_limit = str(loop_limit)
        bounds = str(bounds)

        graph_description = "Avg. Relative Error of $a_{11}$ and $a_{22}$"
        title = ev_type.upper() + " | " + pp_type.upper() + " | " + integration_method + " | BNDS " + bounds + " | " + loop_limit + " cc | " + graph_description
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

    def get_eignevalues(self, mtrx):
        eigenvalues, eigenvectors = np.linalg.eig(mtrx)
        ev_1, ev_2 = eigenvalues
        ev_type = self.get_sim_vars_dict_elt("ev_type")

        # if(ev_type == "ce"):
        # if isinstance(ev_1, complex) and isinstance(ev_2, complex):
        #     ev_1 = math.sqrt((ev_1.real) ** 2 + (ev_1.imag) ** 2 )
        #     ev_2 = math.sqrt((ev_2.real) ** 2 + (ev_2.imag) ** 2 )
        if isinstance(ev_1, complex):
            ev_1 = math.sqrt((ev_1.real) ** 2 + (ev_1.imag) ** 2 )

        if isinstance(ev_2, complex):
            ev_2 = math.sqrt((ev_2.real) ** 2 + (ev_2.imag) ** 2 )


        return ev_1, ev_2



    def organize_data(self, mtrx, guesses, true_vals, labels, case_type):
        a11s, a12s, a21s, a22s = guesses
        a11, a12, a21, a22     = true_vals
        param_1_rel_err, param_2_rel_err = [], []
        data_was_plotted = False
        param_1 = param_2 = 0
        param_l_list, param_2_list = [], []


        if case_type == "MAIN_DIAG_2x2":
            param_1, param_2 = a11, a22
            label_1, label_2 = labels
            param_1_list, param_2_list = a11s, a22s
            param_1_rel_err = abs(param_1_list[-1] - param_1) / abs(param_1)
            param_2_rel_err = abs(param_2_list[-1] - param_2) / abs(param_2)
            avg_rel_err = abs(param_1_rel_err + param_2_rel_err) / 2

        elif(case_type == "OFF_DIAG_2x2"):
            print("Off Diagonal 2x2 - This option is not available yet.")
            exit()

        elif(case_type == "BTM_ROW_2x2"):
            print("Bottom Row 2x2 - This option is not available yet.")
            exit()


        if param_1 != 0 and param_2 != 0 and math.isinf(avg_rel_err) != True:
            ev_1, ev_2 = self.get_eignevalues(mtrx)
            self.plot_points(ev_1, ev_2, avg_rel_err)
            data_was_plotted = True

        else:
            print(label_2 + ":", param_1, label_2 + ":", param_2)
            data_was_plotted = False

        return data_was_plotted


    def plot_points(self, x, y, avg_rel_err):
        fig, ax = self.get_subplots()

        if 1e-1 < avg_rel_err and avg_rel_err < math.inf:
            ax.plot(x, y, ".", color = "red", alpha = 1, zorder = 10, markersize = 10)
            self.set_err_count_dict_elt("over_1e-1", self.get_err_count_dict_elt("over_1e-1") + 1)

        elif 1e-3 < avg_rel_err and avg_rel_err <= 1e-1:
            ax.plot(x, y, ".", color = "red", alpha = 0.3, zorder = 10, markersize = 10)
            self.set_err_count_dict_elt("over_1e-3", self.get_err_count_dict_elt("over_1e-3") + 1)

        elif 1e-5 < avg_rel_err and avg_rel_err <= 1e-3:
            ax.plot(x, y, ".", color = "#800080", alpha = 0.3, zorder = 10, markersize = 10)
            self.set_err_count_dict_elt("over_1e-5", self.get_err_count_dict_elt("over_1e-5") + 1)

        elif 1e-8 < avg_rel_err and avg_rel_err <= 1e-5:
            ax.plot(x, y, ".", color = "blue", alpha = 0.3, zorder = 10, markersize = 10)
            self.set_err_count_dict_elt("over_1e-8", self.get_err_count_dict_elt("over_1e-8") + 1)

        elif 0 <= avg_rel_err and avg_rel_err <= 1e-8:
            ax.plot(x, y, ".", color = "blue", alpha = 1, zorder = 10, markersize = 10)
            self.set_err_count_dict_elt("under_1e-8", self.get_err_count_dict_elt("under_1e-8") + 1)


    def get_date_str(self):
        raw_date = datetime.datetime.now()
        t = time.localtime()
        # current_time = time.strftime("%H.%M", t)
        # date_str = str(raw_date.year) + "." + str(time.strftime('%m')) + "." + str(time.strftime('%d')) + "_" + str(current_time)
        date_str = str(raw_date.year) + "." + str(time.strftime("%m")) + "." + str(time.strftime("%d"))
        return date_str

    def display(self, ev_type, pp_type, integration_method, loop_limit ):
        integration_method = str(integration_method)
        loop_limit = str(loop_limit)
        self.display_ev_graph(ev_type, pp_type, integration_method, loop_limit )



    def display_ev_graph(self, ev_type, pp_type, integration_method, loop_limit):
        fig, ax = self.get_subplots()
        # legend_loc, bbox_x, bbox_y = self.get_legend_loc()
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
        subdir = "../output/" + ev_and_pp_type + "_" + self.get_date_str() + "/"
        os.makedirs(subdir, exist_ok = True)
        fig.savefig(subdir + filename + ".png", dpi = 300)
        plt.close(fig)
