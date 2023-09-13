from matplotlib.lines import Line2D
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import numpy as np
import datetime
import time
import math
import os


class NonLineGraph:
    def __init__(self, *args):
        ev_type, pp_type, bounds, loop_limit, case_type, threshold_vals = args
        self._static_vars_dict = {
            'ev_type': ev_type, 'pp_type': pp_type, 'bounds': bounds, 'loop_limit': loop_limit, "case_type":case_type, "param_label_1": "", "param_label_2": ""
        }
        self._max_trace = 0
        self._threshold_vals = threshold_vals #I.e. [1e-12, 1e-8, 1e-4, 1e-1]
        self._threshold_vals_str = self.convert_to_strings(self._threshold_vals)
        self._err_count_dict = {'under_range_1': 0, 'over_range_1': 0, 'over_range_2': 0, 'over_range_3': 0,  'over_range_4': 0 }
        self.set_param_labels()
        self.set_tr_det_subplots(str(loop_limit), str(bounds), self.get_static_vars_dict_elt("case_type"))
        self.set_ev_subplots(ev_type, pp_type, case_type, str(bounds), str(loop_limit))

    def get_threshold_vals_elt(self, index):
        return self._threshold_vals[index]

    def get_threshold_vals_str_elt(self, index):
        return self._threshold_vals_str[index]

    def convert_to_strings(self, input_list):
        string_list = ["{:.1e}".format(item).replace("e-0", "e-").replace(".0e", "e") for item in input_list]
        return string_list

    def get_eignevalues(self, A):
        eigenvalues, eigenvectors = np.linalg.eig(A)
        ev_1, ev_2 = eigenvalues
        ev_type = self.get_static_vars_dict_elt("ev_type")

        if isinstance(ev_1, complex):
            ev_1 = math.sqrt((ev_1.real) ** 2 + (ev_1.imag) ** 2 )

        if isinstance(ev_2, complex):
            ev_2 = math.sqrt((ev_2.real) ** 2 + (ev_2.imag) ** 2 )
        return ev_1, ev_2


    def get_err_count_dict_elt(self, key):
        return self._err_count_dict[key]

    def get_ev_subplots(self):
        return self._ev_fig, self._ev_ax

    def get_max_trace(self):
        return self._max_trace

    def get_static_vars_dict(self):
        return self._static_vars_dict

    def get_static_vars_dict_elt(self, key):
        return self._static_vars_dict[key]

    def get_tr_det_subplots(self):
        return self._td_fig, self._td_ax

    ##############################################

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

    def set_static_vars_dict_elt(self, key, value):
        self._static_vars_dict[key] = value

    def set_err_count_dict_elt(self, key, value):
        self._err_count_dict[key] = value

    def set_max_trace(self, new_max_trace):
        new_max_trace = abs(new_max_trace)

        if new_max_trace > self._max_trace:
            self._max_trace = new_max_trace


    def set_ev_subplots(self, ev_type, pp_type, case_type, bounds, loop_limit):
        self._ev_fig, self._ev_ax = plt.subplots()
        param_label_1, param_label_2 = self.get_static_vars_dict_elt("param_label_1"), self.get_static_vars_dict_elt("param_label_2")
        graph_description = "Avg. Relative Error of " + param_label_1 + " and " + param_label_2
        title = ev_type.upper() + " | " + pp_type.upper() + " | " + case_type.upper() + " | BNDS " + bounds + " | " + loop_limit + " CC | " + graph_description
        self._ev_ax.set_title(label = title, pad = 30, fontsize = 15)
        self._ev_ax.set_xlabel("$\u03BB_{1}$", loc = "right", fontsize = 14)
        self._ev_ax.set_ylabel("$\u03BB_{2}$", loc = "top", fontsize = 14)
        self._ev_fig.set_size_inches(16, 8)
        self._ev_ax = plt.gca()

        # Hide two spines
        self._ev_ax.spines["right"].set_color("none")
        self._ev_ax.spines["top"].set_color("none")

        # Move bottom and left spine to 0, 0
        self._ev_ax.spines["bottom"].set_position(("data", 0))
        self._ev_ax.spines["left"].set_position(("data", 0))

        # Move ticks positions
        self._ev_ax.xaxis.set_ticks_position("bottom")
        self._ev_ax.yaxis.set_ticks_position("left")

        self._ev_ax.plot(1, 0, ">k", transform = self._ev_ax.get_yaxis_transform(), clip_on = False)
        self._ev_ax.plot(0, 1, "^k", transform = self._ev_ax.get_xaxis_transform(), clip_on = False)


    def set_tr_det_subplots(self, loop_limit, bounds, case_type):
        self._td_fig, self._td_ax = plt.subplots()
        ev_type               = self.get_static_vars_dict_elt("ev_type")
        pp_type               = self.get_static_vars_dict_elt("pp_type")
        bounds                = str(self.get_static_vars_dict_elt("bounds"))
        param_label_1 = self.get_static_vars_dict_elt("param_label_1")
        param_label_2 = self.get_static_vars_dict_elt("param_label_2")

        graph_description = "Avg. Relative Error of " + param_label_1 + " and " + param_label_2
        title = ev_type.upper() + " | " + pp_type.upper() + " | " + case_type.upper() + " | BNDS " + bounds + " | " + loop_limit + " CC | " + graph_description
        self._td_ax.set_title(label = title, pad = 30, fontsize = 15)
        self._td_ax.set_xlabel("Tr", loc = "right", fontsize = 14)
        self._td_ax.set_ylabel("Det", loc = "top", fontsize = 14)
        self._td_fig.set_size_inches(16, 8)
        self._td_ax = plt.gca()

        # Hide two spines
        self._td_ax.spines["right"].set_color("none")
        self._td_ax.spines["top"].set_color("none")

        # Move bottom and left spine to 0, 0
        self._td_ax.spines["bottom"].set_position(("data", 0))
        self._td_ax.spines["left"].set_position(("data", 0))

        # Move ticks positions
        self._td_ax.xaxis.set_ticks_position("bottom")
        self._td_ax.yaxis.set_ticks_position("left")

        self._td_ax.plot(1, 0, ">k", transform = self._td_ax.get_yaxis_transform(), clip_on = False)
        self._td_ax.plot(0, 1, "^k", transform = self._td_ax.get_xaxis_transform(), clip_on = False)

    def organize_data(self, A, param_estimates, avg_param_errors, trace_A, det_A, case_type):
        param_1_estimates, param_2_estimates = param_estimates
        # param_1_abs_err, param_2_abs_err, param_1_rel_err, param_2_rel_err = param_errors
        a11, a12, a21, a22  = A[0,0], A[0, 1], A[1,0], A[1,1]
        param_1_rel_err, param_2_rel_err = [], []
        data_was_plotted = False
        param_1 = param_2 = 0
        param_label_1 = self.get_static_vars_dict_elt("param_label_1")
        param_label_2 = self.get_static_vars_dict_elt("param_label_2")
        avg_abs_param_err, avg_rel_param_err =  avg_param_errors
        nth_avg_rel_param_err = avg_rel_param_err[-1]

        if case_type == "main_diagonal":
            param_1, param_2 = a11, a22


        elif(case_type == "anti-diagonal"):
            param_1, param_2 = a12, a21


        elif(case_type == "left_column"):
            param_1, param_2 = a11, a21

        elif(case_type == "right_column"):
            param_1, param_2 = a12, a22

        if param_1 != 0 and param_2 != 0 and math.isinf(nth_avg_rel_param_err) != True:
            td_fig, td_ax =  self.get_tr_det_subplots()
            self.plot_points(trace_A, det_A, nth_avg_rel_param_err, (td_fig, td_ax))
            self.set_max_trace(trace_A)

            ev_fig, ev_ax = self.get_ev_subplots()
            ev_1, ev_2 = self.get_eignevalues(A)
            self.update_error_count(nth_avg_rel_param_err)
            self.plot_points(ev_1, ev_2, nth_avg_rel_param_err, (ev_fig, ev_ax))
            data_was_plotted = True

        else:
            print(param_label_1 + ":", param_1, param_label_2 + ":", param_2)
            print("This matrix is not plottable.", "\n")
            data_was_plotted = False

        return data_was_plotted

    def update_error_count(self, nth_avg_rel_param_err):
        # if 1e-1 < ... < math.inf:
        if self.get_threshold_vals_elt(3) < nth_avg_rel_param_err and nth_avg_rel_param_err < math.inf:
            self.set_err_count_dict_elt("over_range_4", self.get_err_count_dict_elt("over_range_4") + 1) # over_1e-1

        # elif 1e-4 < ... <= 1e-1:
        elif self.get_threshold_vals_elt(2) < nth_avg_rel_param_err and nth_avg_rel_param_err <= self.get_threshold_vals_elt(2):
            self.set_err_count_dict_elt("over_range_3", self.get_err_count_dict_elt("over_range_3") + 1) # over_1e-4

        # elif 1e-8 < ... <= 1e-4:
        elif self.get_threshold_vals_elt(1) < nth_avg_rel_param_err and nth_avg_rel_param_err <= self.get_threshold_vals_elt(2):
            self.set_err_count_dict_elt("over_range_2", self.get_err_count_dict_elt("over_range_2") + 1) # over_1e-8

        #elif 1e-12 < ...  <= 1e-8:
        elif self.get_threshold_vals_elt(0) < nth_avg_rel_param_err and nth_avg_rel_param_err <= self.get_threshold_vals_elt(1):
            self.set_err_count_dict_elt("over_range_1", self.get_err_count_dict_elt("over_range_1") + 1) # over_1e-12

        #elif 0 <= ...  <= 1e-12:
        elif 0 <= nth_avg_rel_param_err and nth_avg_rel_param_err <= self.get_threshold_vals_elt(0):
            self.set_err_count_dict_elt("under_range_1", self.get_err_count_dict_elt("under_range_1") + 1) #under_1e-12



    def plot_points(self, x, y, nth_avg_rel_param_err, subplots):
        fig, ax = subplots

        # if 1e-1 < ... < math.inf:
        if self.get_threshold_vals_elt(3) < nth_avg_rel_param_err and nth_avg_rel_param_err < math.inf:
            ax.plot(x, y, ".", color = "red", alpha = 1, zorder = 10, markersize = 10)
            # self.set_err_count_dict_elt("over_range_4", self.get_err_count_dict_elt("over_range_4") + 1) # over_1e-1

        # elif 1e-4 < ... <= 1e-1:
        elif self.get_threshold_vals_elt(2) < nth_avg_rel_param_err and nth_avg_rel_param_err <= self.get_threshold_vals_elt(2):
            ax.plot(x, y, ".", color = "red", alpha = 0.3, zorder = 10, markersize = 10)
            # self.set_err_count_dict_elt("over_range_3", self.get_err_count_dict_elt("over_range_3") + 1) # over_1e-4

        # elif 1e-8 < ... <= 1e-4:
        elif self.get_threshold_vals_elt(1) < nth_avg_rel_param_err and nth_avg_rel_param_err <= self.get_threshold_vals_elt(2):
            ax.plot(x, y, ".", color = "#800080", alpha = 0.3, zorder = 10, markersize = 10)
            # self.set_err_count_dict_elt("over_range_2", self.get_err_count_dict_elt("over_range_2") + 1) # over_1e-8

        #elif 1e-12 < ...  <= 1e-8:
        elif self.get_threshold_vals_elt(0) < nth_avg_rel_param_err and nth_avg_rel_param_err <= self.get_threshold_vals_elt(1):
            ax.plot(x, y, ".", color = "blue", alpha = 0.3, zorder = 10, markersize = 10)
            # self.set_err_count_dict_elt("over_range_1", self.get_err_count_dict_elt("over_range_1") + 1) # over_1e-12

        #elif 0 <= ...  <= 1e-12:
        elif 0 <= nth_avg_rel_param_err and nth_avg_rel_param_err <= self.get_threshold_vals_elt(0):
            ax.plot(x, y, ".", color = "blue", alpha = 1, zorder = 10, markersize = 10)
            # self.set_err_count_dict_elt("under_range_1", self.get_err_count_dict_elt("under_range_1") + 1) #under_1e-12


    def graph_parabola(self):
        fig, ax = self.get_tr_det_subplots()
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


    def display(self, ev_type, pp_type, case_type, loop_limit ):
        bounds                = str(self.get_static_vars_dict_elt("bounds"))
        loop_limit            = str(self.get_static_vars_dict_elt("loop_limit"))

        self.display_tr_det_graph(ev_type, pp_type, case_type, loop_limit )
        self.display_ev_graph(ev_type, pp_type, case_type, loop_limit )
        self.display_bar_graph(ev_type, pp_type, case_type, loop_limit, bounds)
        self.display_pie_graph(ev_type, pp_type, case_type, loop_limit, bounds)



    def display_ev_graph(self, ev_type, pp_type, case_type, loop_limit):
        fig, ax = self.get_ev_subplots()
        bbox_x = bbox_y = 0
        legend_loc = "center left"

        label_1 = "0 <= $x̄_{Err}$ <="  + self.get_threshold_vals_str_elt(0)
        label_2 = self.get_threshold_vals_str_elt(0) + " < $x̄_{Err}$ <= " + self.get_threshold_vals_str_elt(1)
        label_3 = self.get_threshold_vals_str_elt(1) + " < $x̄_{Err}$ <= " + self.get_threshold_vals_str_elt(2)
        label_4 = self.get_threshold_vals_str_elt(2) + " < $x̄_{Err}$ <=" + self.get_threshold_vals_str_elt(3)
        label_5 = self.get_threshold_vals_str_elt(3)+ " < $x̄_{Err}$ < ∞"

        custom_handles = [
            Line2D([0], [0], marker = "o", markerfacecolor = "r", color = "w", alpha = 1, markersize = 7, label = label_5),
            Line2D([0], [0], marker = "o", markerfacecolor = "r", color = "w", alpha = 0.3, markersize = 7, label = label_4),
            Line2D([0], [0], marker = "o", markerfacecolor = "#800080", color = "w", alpha = 0.3, markersize = 7, label = label_3),
            Line2D([0], [0], marker = "o", markerfacecolor = "b", color = "w", alpha = 0.3, markersize = 7, label = label_2),
            Line2D([0], [0], marker = "o", markerfacecolor = "b", color = "w", markersize = 7, label = label_1)]


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
        fig.savefig(subdir + filename + ".jpg", dpi = 300)
        plt.close(fig)


    def display_tr_det_graph(self, ev_type, pp_type, case_type, loop_limit):
        fig, ax = self.get_tr_det_subplots()
        self.graph_parabola()

        label_1 = '0 <= $x̄_{Err}$ <='  + self.get_threshold_vals_str_elt(0)
        label_2 = self.get_threshold_vals_str_elt(0) + " < $x̄_{Err}$ <= " + self.get_threshold_vals_str_elt(1)
        label_3 =  self.get_threshold_vals_str_elt(1) + " < $x̄_{Err}$ <= " + self.get_threshold_vals_str_elt(2)
        label_4 =  self.get_threshold_vals_str_elt(2) + " < $x̄_{Err}$ <=" + self.get_threshold_vals_str_elt(3)
        label_5 =   self.get_threshold_vals_str_elt(3)+ " < $x̄_{Err}$ < ∞"

        custom_handles = [
            Line2D([0], [0], marker = "o", markerfacecolor = "r", color = "w", alpha = 1, markersize = 7, label = label_5),
            Line2D([0], [0], marker = "o", markerfacecolor = "r", color = "w", alpha = 0.3, markersize = 7, label = label_4),
            Line2D([0], [0], marker = "o", markerfacecolor = "#800080", color = "w", alpha = 0.3, markersize = 7, label = label_3),
            Line2D([0], [0], marker = "o", markerfacecolor = "b", color = "w", alpha = 0.3, markersize = 7, label = label_2),
            Line2D([0], [0], marker = "o", markerfacecolor = "b", color = "w", markersize = 7, label = label_1),
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
        subdir = "../output/" + ev_and_pp_type + "_" + case_type + "_" + self.get_date_str() + "/"
        os.makedirs(subdir, exist_ok = True)
        fig.savefig(subdir + filename + ".jpg", dpi = 300)
        plt.close(fig)


    def display_pie_graph(self, ev_type, pp_type, case_type, loop_limit, bounds):
        v1 = self.get_err_count_dict_elt("under_range_1") # under_1e-12
        v2 = self.get_err_count_dict_elt("over_range_1") # over_1e-12
        v3 = self.get_err_count_dict_elt("over_range_2") # over_1e-8
        v4 = self.get_err_count_dict_elt("over_range_3") # over_1e-4
        v5 = self.get_err_count_dict_elt("over_range_4") # over_1e-1

        values_ph = [v1, v2, v3, v4, v5]
        # slice_colors_ph = ['red', '#ffb2b2', '#d8b2d8', '#b2b2ff', 'blue' ]
        slice_colors_ph = [ 'blue',  '#b2b2ff', '#d8b2d8', '#ffb2b2', 'red']
        values, slice_colors, names = [], [], []
        param_label_1 = self.get_static_vars_dict_elt("param_label_1")
        param_label_2 = self.get_static_vars_dict_elt("param_label_2")

        i = 0
        while i < len(values_ph):
            if values_ph[i] != 0:
                values.append(values_ph[i])
                slice_colors.append(slice_colors_ph[i])
            i = i + 1
        values       = np.array(values)
        slice_colors = np.array(slice_colors)

        label_1 = '0 <= $x̄_{Err}$ <='  + self.get_threshold_vals_str_elt(0)
        label_2 = self.get_threshold_vals_str_elt(0) + " < $x̄_{Err}$ <= " + self.get_threshold_vals_str_elt(1)
        label_3 = self.get_threshold_vals_str_elt(1) + " < $x̄_{Err}$ <= " + self.get_threshold_vals_str_elt(2)
        label_4 = self.get_threshold_vals_str_elt(2) + " < $x̄_{Err}$ <=" + self.get_threshold_vals_str_elt(3)
        label_5 = self.get_threshold_vals_str_elt(3)+ " < $x̄_{Err}$ < ∞"


        custom_handles = [
            Line2D([0], [0], marker = "o", markerfacecolor = "r", color = "w", alpha = 1, markersize = 7, label = label_5),
            Line2D([0], [0], marker = "o", markerfacecolor = "r", color = "w", alpha = 0.3, markersize = 7, label = label_4),
            Line2D([0], [0], marker = "o", markerfacecolor = "#800080", color = "w", alpha = 0.3, markersize = 7, label = label_3),
            Line2D([0], [0], marker = "o", markerfacecolor = "b", color = "w", alpha = 0.3, markersize = 7, label = label_2),
            Line2D([0], [0], marker = "o", markerfacecolor = "b", color = "w", markersize = 7, label = label_1)]

        fig, ax = plt.subplots(figsize = (13, 10))
        plt.pie(values, autopct = '%1.1f%%', colors = slice_colors, radius = 1.5, textprops = {'color': 'white', 'weight': 'bold', 'fontsize': 15})
        fig.subplots_adjust(right = 0.7, top = 0.85)
        plt.legend(handles = custom_handles, loc = "lower right", bbox_to_anchor = (1.51, 0), borderpad = 1, fontsize = "15")

        # For file output
        graph_description = "Avg. Relative Error of " + param_label_1 + " and " + param_label_2
        title = ev_type.upper() + " | " + pp_type.upper() + " | " + case_type.upper() + " | BNDS " + bounds + " | " + loop_limit + " CC | " + graph_description
        fig.suptitle(title, fontsize = 15)
        ev_and_pp_type = ev_type + "_" + pp_type
        subdir = "../output/" + ev_and_pp_type + "_" + case_type + "_" +  self.get_date_str() + "/"
        filename = 'pie_graph' + "_" + ev_and_pp_type + "_" + loop_limit + "_" + "cc"
        fig.savefig(subdir + filename + ".jpg", dpi = 300)
        plt.close(fig)


    def display_bar_graph(self, ev_type, pp_type, case_type, loop_limit, bounds):
        fig, ax = plt.subplots(figsize = (16, 8))
        ax.set_axisbelow(True)
        ax.grid(color = '#cccccc', linestyle = 'dashed')
        ax.set_ylabel("Average Relative Error")
        ax.set_xlabel("Frequency")
        ax.xaxis.set_major_locator(ticker.MaxNLocator(integer = True))

        param_label_1 = self.get_static_vars_dict_elt("param_label_1")
        param_label_2 = self.get_static_vars_dict_elt("param_label_2")

        # Set-up labels
        label_1 = '0 <= $x̄_{Err}$ <='  + self.get_threshold_vals_str_elt(0)
        label_2 = self.get_threshold_vals_str_elt(0) + " < $x̄_{Err}$ <= " + self.get_threshold_vals_str_elt(1)
        label_3 = self.get_threshold_vals_str_elt(1) + " < $x̄_{Err}$ <= " + self.get_threshold_vals_str_elt(2)
        label_4 = self.get_threshold_vals_str_elt(2) + " < $x̄_{Err}$ <=" + self.get_threshold_vals_str_elt(3)
        label_5 = self.get_threshold_vals_str_elt(3)+ " < $x̄_{Err}$ < ∞"
        param_labels = np.array([label_1, label_2,  label_3, label_4, label_5])

        # Set-up values
        v1 = self.get_err_count_dict_elt("under_range_1") # under_1e-12
        v2 = self.get_err_count_dict_elt("over_range_1") # over_1e-12
        v3 = self.get_err_count_dict_elt("over_range_2") # over_1e-8
        v4 = self.get_err_count_dict_elt("over_range_3") # over_1e-4
        v5 = self.get_err_count_dict_elt("over_range_4") # over_1e-1

        values = np.array([v1, v2, v3, v4, v5])
        total_bars = plt.barh(param_labels, values)
        total_bars[0].set_color('blue')
        total_bars[1].set_color('#b2b2ff') #The HEX color value ....
        total_bars[2].set_color('#d8b2d8') #The HEX color value ....
        total_bars[3].set_color('#ffb2b2') #The HEX color value w/o any alpha adjustment
        total_bars[4].set_color('red')

        # For Output
        graph_description = "Avg. Relative Error of " + param_label_1 + " and " + param_label_2
        title = ev_type.upper() + " | " + pp_type.upper() + " | " + case_type.upper() + " | BNDS " + bounds + " | " + loop_limit + " CC | " + graph_description
        fig.suptitle(title, fontsize = 15)
        ev_and_pp_type = ev_type + "_" + pp_type
        subdir = "../output/" + ev_and_pp_type + "_" + case_type + "_" +  self.get_date_str() + "/"
        filename = 'bar_graph'+ "_" + ev_and_pp_type + "_" + loop_limit + "_" + "cc"
        plt.savefig(subdir + filename + '.jpg', dpi = 300)
        plt.close(fig)
