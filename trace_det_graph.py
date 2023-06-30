from matplotlib.lines import Line2D
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import numpy as np
import datetime
import time
import math
import os


class TraceDetGraph:
    def __init__(self, *args):
        ev_type, pp_type, integration_method, bounds, loop_limit, case_type = args
        self._static_vars_dict = {
            'ev_type': ev_type, 'pp_type': pp_type, 'integration_method': integration_method, 'bounds': bounds, 'loop_limit': loop_limit, "case_type":case_type
        }
        self._max_trace = 0
        self._err_count_dict = {'under_1e-8': 0, 'over_1e-8': 0, 'over_1e-5': 0, 'over_1e-3': 0,  'over_1e-1': 0 }
        self.set_subplots(str(loop_limit), integration_method, str(bounds))

    def get_static_vars_dict(self):
        return self._static_vars_dict

    def get_static_vars_dict_elt(self, key):
        return self._static_vars_dict[key]

    def get_subplots(self):
        return self._fig, self._ax
    #
    def get_max_trace(self):
        return self._max_trace

    def get_err_count_dict_elt(self, key):
        return self._err_count_dict[key]

    def set_err_count_dict_elt(self, key, value):
        self._err_count_dict[key] = value

    def set_max_trace(self, new_max_trace):
        new_max_trace = abs(new_max_trace)

        if new_max_trace > self._max_trace:
            self._max_trace = new_max_trace


    def set_subplots(self, loop_limit, integration_method, bounds):
        self._fig, self._ax = plt.subplots()
        ev_type               = self.get_static_vars_dict_elt("ev_type")
        pp_type               = self.get_static_vars_dict_elt("pp_type")
        bounds                = str(self.get_static_vars_dict_elt("bounds"))


        graph_description = "Avg. Relative Error of $a_{11}$ and $a_{22}$"
        title = ev_type.upper() + " | " + pp_type.upper() + " | " + integration_method + " | BNDS " + bounds + " | " + loop_limit + " cc | " + graph_description
        self._ax.set_title(label = title, pad = 30, fontsize = 15)
        self._ax.set_xlabel("Tr", loc = "right", fontsize = 14)
        self._ax.set_ylabel("Det", loc = "top", fontsize = 14)
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



    def organize_data(self, guesses, true_vals, param_labels, true_trace, true_det, case_type):

        a11s, a12s, a21s, a22s = guesses
        a11, a12, a21, a22     = true_vals
        param_1_rel_err, param_2_rel_err = [], []
        data_was_plotted = False
        param_1 = param_2 = 0
        param_l_list, param_2_list = [], []

        if case_type == "MAIN_DIAG_2x2":
            param_1, param_2 = a11, a22
            label_1, label_2 = param_labels
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
            print("Tr:", true_trace)
            print("Det:", true_det, "\n")
            print(label_1, "Avg. rel. err:", param_1_rel_err)
            print(label_2, "Avg. rel. err:", param_2_rel_err)
            print("Avg. Rel. Err:", avg_rel_err, "\n")
            self.plot_points(true_trace, true_det, avg_rel_err)
            self.set_max_trace(true_trace)
            data_was_plotted = True

        else:
            print(lable_1 + ":", param_1, label_2 + ":", param_2)
            print("This matrix is not plottable.", "\n")
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


    def graph_parabola(self):
        fig, ax = self.get_subplots()
        max_trace = self.get_max_trace()

        if abs(max_trace) == 0:
            max_trace = 10
        x = np.linspace(-(max_trace), max_trace, 100)
        y = x**2 / 4
        ax.plot(x, y, linewidth = 1.5, c = "orange")

    def get_date_str(self):
        raw_date = datetime.datetime.now()
        t = time.localtime()
        date_str = str(raw_date.year) + "." + str(time.strftime("%m")) + "." + str(time.strftime("%d"))
        return date_str

    def display(self, ev_type, pp_type, integration_method, loop_limit ):

        bounds                = str(self.get_static_vars_dict_elt("bounds"))
        loop_limit            = str(self.get_static_vars_dict_elt("loop_limit"))
        self.display_tr_det_graph(ev_type, pp_type, integration_method, loop_limit )
        self.display_bar_graph(ev_type, pp_type, integration_method, loop_limit, bounds)
        self.display_pie_graph(ev_type, pp_type, integration_method, loop_limit, bounds)

    def display_tr_det_graph(self, ev_type, pp_type, integration_method, loop_limit):
        fig, ax = self.get_subplots()
        self.graph_parabola()

        custom_handles = [
            Line2D([0], [0], marker = "o", markerfacecolor = "r", color = "w", alpha = 1, markersize = 7, label = "1e-1 < $x̄_{Err}$ < ∞"),
            Line2D([0], [0], marker = "o", markerfacecolor = "r", color = "w", alpha = 0.3, markersize = 7, label = "1e-3 < $x̄_{Err}$ <= 1e-1"),
            Line2D([0], [0], marker = "o", markerfacecolor = "#800080", color = "w", alpha = 0.3, markersize = 7, label = "1e-5 < $x̄_{Err}$ <= 1e-3"),
            Line2D([0], [0], marker = "o", markerfacecolor = "b", color = "w", alpha = 0.3, markersize = 7, label = "1e-8 < $x̄_{Err}$ <= 1e-5"),
            Line2D([0], [0], marker = "o", markerfacecolor = "b", color = "w", markersize = 7, label = "0 <= $x̄_{Err}$ <= 1e-8"),
            Line2D([0], [0], color = "orange", alpha = 1, lw = 3, label = "T\N{SUPERSCRIPT TWO} - 4D = 0")]


        if(pp_type == "saddle"):
            ax.legend(handles = custom_handles, loc = "upper left", bbox_to_anchor = (-.15, 1.15), borderpad = 1)


        # #mid-left
        elif(pp_type == "source" or pp_type == "sp_source" ):
            ax.legend(handles = custom_handles, loc = "center left", bbox_to_anchor = (-.13, .4), borderpad = 1)

        #mid-right
        elif(pp_type == "sink" or pp_type == "sp_sink" or pp_type == "center"):
            ax.legend(handles = custom_handles, loc = "center right", bbox_to_anchor = (1.1, .4), borderpad = 1)

        ev_and_pp_type = ev_type + "_" + pp_type
        filename = "tr_det_graph" + "_" + ev_and_pp_type + "_" + loop_limit + "_" + "cc"
        subdir = "../output/" + ev_and_pp_type + "_" + self.get_date_str() + "/"
        os.makedirs(subdir, exist_ok = True)
        fig.savefig(subdir + filename + ".png", dpi = 300)
        plt.close(fig)

    def display_bar_graph(self, ev_type, pp_type, integration_method, loop_limit, bounds):
        fig, ax = plt.subplots(figsize = (16, 8))
        ax.set_axisbelow(True)
        ax.grid(color = '#cccccc', linestyle = 'dashed')
        ax.set_ylabel("Average Relative Error")
        ax.set_xlabel("Frequency")
        ax.xaxis.set_major_locator(ticker.MaxNLocator(integer = True))

        param_labels = np.array(["0 <= $x̄_{Err}$ <= 1e-8", "1e-8 < $x̄_{Err}$ <= 1e-5",  "1e-5 < $x̄_{Err}$ <= 1e-3", "1e-3 < $x̄_{Err}$ <= 1e-1", "1e-1 < $x̄_{Err}$ < ∞"])
        v5 = self.get_err_count_dict_elt("over_1e-1")
        v4 = self.get_err_count_dict_elt("over_1e-3")
        v3 = self.get_err_count_dict_elt("over_1e-5")
        v2 = self.get_err_count_dict_elt("over_1e-8")
        v1 = self.get_err_count_dict_elt("under_1e-8")
        values = np.array([v1, v2, v3, v4, v5])

        total_bars = plt.barh(param_labels, values)
        total_bars[0].set_color('blue')
        total_bars[1].set_color('#b2b2ff') #The HEX color value ....
        total_bars[2].set_color('#d8b2d8') #The HEX color value ....
        total_bars[3].set_color('#ffb2b2') #The HEX color value w/o any alpha adjustment
        total_bars[4].set_color('red')

        # For Output

        graph_description = "Avg. Relative Error of $a_{11}$ and $a_{22}$"
        title = ev_type.upper() + " | " + pp_type.upper() + " | " + integration_method + " | BNDS " + bounds + " | " + loop_limit + " cc | " + graph_description
        fig.suptitle(title, fontsize = 15)
        ev_and_pp_type = ev_type + "_" + pp_type
        subdir = "../output/" + ev_and_pp_type + "_" + self.get_date_str() + "/"
        filename = 'bar_graph'+ "_" + ev_and_pp_type + "_" + loop_limit + "_" + "cc"
        plt.savefig(subdir + filename + '.png', dpi = 300)
        plt.close(fig)


    def display_pie_graph(self, ev_type, pp_type, integration_method, loop_limit, bounds):
        v1 = self.get_err_count_dict_elt("over_1e-1")
        v2 = self.get_err_count_dict_elt("over_1e-3")
        v3 = self.get_err_count_dict_elt("over_1e-5")
        v4 = self.get_err_count_dict_elt("over_1e-8")
        v5 = self.get_err_count_dict_elt("under_1e-8")
        values_ph = [v1, v2, v3, v4, v5]
        slice_colors_ph = ['red', '#ffb2b2', '#d8b2d8', '#b2b2ff', 'blue' ]
        values, slice_colors, names = [], [], []
        i = 0
        while i < len(values_ph):
            if values_ph[i] != 0:
                values.append(values_ph[i])
                slice_colors.append(slice_colors_ph[i])

            i = i + 1
        values       = np.array(values)
        slice_colors = np.array(slice_colors)


        custom_handles = [
            Line2D([0], [0], marker = "o", markerfacecolor = "r", color = "w", alpha = 1, markersize = 7, label = "1e-1 < $x̄_{Err}$ < ∞"),
            Line2D([0], [0], marker = "o", markerfacecolor = "r", color = "w", alpha = 0.3, markersize = 7, label = "1e-3 < $x̄_{Err}$ <= 1e-1"),
            Line2D([0], [0], marker = "o", markerfacecolor = "#800080", color = "w", alpha = 0.3, markersize = 7, label = "1e-5 < $x̄_{Err}$ <= 1e-3"),
            Line2D([0], [0], marker = "o", markerfacecolor = "b", color = "w", alpha = 0.3, markersize = 7, label = "1e-8 < $x̄_{Err}$ <= 1e-5"),
            Line2D([0], [0], marker = "o", markerfacecolor = "b", color = "w", markersize = 7, label = "0 <= $x̄_{Err}$ <= 1e-8")]

        fig, ax = plt.subplots(figsize = (13, 10))
        plt.pie(values, autopct = '%1.1f%%', colors = slice_colors, radius = 1.5, textprops = {'color': 'white', 'weight': 'bold', 'fontsize': 15})
        fig.subplots_adjust(right = 0.7, top = 0.85)
        plt.legend(handles = custom_handles, loc = "lower right", bbox_to_anchor = (1.51, 0), borderpad = 1, fontsize = "15")

        # For file output
        graph_description = "Avg. Relative Error of $a_{11}$ and $a_{22}$"
        title = ev_type.upper() + " | " + pp_type.upper() + " | " + integration_method + " | BNDS " + bounds + " | " + loop_limit + " cc | " + graph_description
        fig.suptitle(title, fontsize = 15)
        ev_and_pp_type = ev_type + "_" + pp_type
        subdir = "../output/" + ev_and_pp_type + "_" + self.get_date_str() + "/"
        filename = 'pie_graph' + "_" + ev_and_pp_type + "_" + loop_limit + "_" + "cc"
        fig.savefig(subdir + filename + ".png", dpi = 300)
        plt.close(fig)
