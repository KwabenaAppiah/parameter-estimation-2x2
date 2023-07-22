import matplotlib.pyplot as plt
import numpy as np
import datetime
import time
import math
import os


class LineGraph:
    def __init__(self, *args):
        ev_type, pp_type, bounds, loop_limit, case_type = args

        # FOR COMP ONLY
        self._true_params = []
        self._has_bad_matrices = False
        self._static_vars_dict = {
            'ev_type': ev_type, 'pp_type': pp_type, 'bounds': str(bounds), 'loop_limit': loop_limit, "case_type": case_type, "param_label_1": "", "param_label_2": ""
        }
        self.set_param_labels()
        self.set_subplots_comp()



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


    def organize_data(self, t, update_times, A, param_estimates, param_errors, avg_param_errors, curr_indx, threshold, S_lists, U_lists):

        # Init. static vars
        ev_type   = self.get_static_vars_dict_elt('ev_type')
        pp_type   = self.get_static_vars_dict_elt('pp_type')
        bounds    = self.get_static_vars_dict_elt('bounds')
        case_type = self.get_static_vars_dict_elt('case_type')
        loop_limit = self.get_static_vars_dict_elt('loop_limit')
        param_labels = (self.get_static_vars_dict_elt("param_label_1"), self.get_static_vars_dict_elt("param_label_2"))
        static_args = (ev_type, pp_type, bounds, loop_limit, param_labels, case_type)

        curr_indx = str(curr_indx)
        matrix_type = ""
        was_content_plotted = False

        # Parameter stuff
        # a11_estimates, a12_estimates, a21_estimates, a22_estimates = param_estimates
        param_1_estimates, param_2_estimates = param_estimates
        a11, a12, a21, a22  = A[0,0], A[0, 1], A[1,0], A[1,1]
        param_1_abs_err, param_2_abs_err, param_1_rel_err, param_2_rel_err = param_errors
        true_params = (a11, a12, a21, a22)
        true_params_str   = str(true_params)
        true_params_title = self.get_true_params_title(true_params)
        param_1 = param_2 = 0
        param_l_list, param_2_estimates = [], []
        # avg_abs_param_err = self.get_avg_list(param_1_abs_err, param_2_abs_err)
        # avg_rel_param_err = self.get_avg_list(param_1_rel_err, param_2_rel_err)
        avg_abs_param_err, avg_rel_param_err =  avg_param_errors
        nth_avg_rel_param_err = avg_rel_param_err[-1]
        linthresh_value = .000000000000001


        # Solution stuff
        # x_estimates, y_estimates, xt_estimates, yt_estimates = S_lists
        # x_sig_err_list, y_sig_err_list   = U_lists

        if case_type == "main_diagonal":
            param_1, param_2 = a11, a22
            # param_1_estimates, param_2_estimates = a11_estimates, a22_estimates

        elif(case_type == "anti-diagonal"):
            param_1, param_2 = a12, a21


        elif(case_type == "left_column"):
            param_1, param_2 = a11, a21

        elif(case_type == "right_column"):
            param_1, param_2 = a12, a22



        # if param_1 != 0 and param_2 != 0:
        #     if(pp_type == "sink" or pp_type == "sp_sink"):
        #         linthresh_value = .000000000000001
        #     else:
        #         linthresh_value = .1

        if nth_avg_rel_param_err >= threshold:
            matrix_type = "bad"
            # self.plot_graph_comp(total_avg_rel_param_err_list)
            self.set_has_bad_matrices_comp(True)
            data_was_plotted = True

        elif nth_avg_rel_param_err < threshold:
            matrix_type = "good"
            data_was_plotted = True

        else:
            data_was_plotted = False

        # # Output regardless
        self.write_to_file("mtrx_" + curr_indx + "|" + true_params_str, case_type, "all")
        self.write_to_file("mtrx_" + curr_indx + "|" + true_params_str, case_type, matrix_type)

        # Parameter eror-based graphs
        self.display_avg_param_err_indv(t, static_args, true_params_title, avg_rel_param_err, curr_indx, linthresh_value, matrix_type,"Relative")
        self.display_avg_param_err_indv(t, static_args, true_params_title, avg_abs_param_err, curr_indx, linthresh_value, matrix_type, "Absolute")
        self.display_param_err(t, static_args, param_1_rel_err, param_2_rel_err, true_params_title, curr_indx, linthresh_value, matrix_type, "Relative")
        self.display_param_err(t, static_args, param_1_rel_err, param_2_rel_err, true_params_title, curr_indx, linthresh_value, matrix_type, "Absolute")

        # # # Solutions-based graphs
        self.display_sol_xy_over_t(t, static_args, S_lists, true_params_title, curr_indx, linthresh_value, matrix_type)
        self.display_sol_x(t, static_args, S_lists, true_params_title, curr_indx, linthresh_value, matrix_type)
        self.display_sol_y(t, static_args, S_lists, true_params_title, curr_indx, linthresh_value,matrix_type)
        self.display_sol_xy(static_args, S_lists, true_params_title, curr_indx, linthresh_value, matrix_type)
        self.display_sol_signal_err(t, static_args, U_lists, true_params_title, curr_indx, linthresh_value, matrix_type)
        self.display_sol_signal_err_split(t, static_args, U_lists, true_params_title, curr_indx, linthresh_value, matrix_type)


        return was_content_plotted

    def write_to_file(self, line_str, case_type, matrix_type):
        output = line_str + "\n" + "\n"
        ev_and_pp_type = self.get_static_vars_dict_elt("ev_type") + "_" + self.get_static_vars_dict_elt("pp_type")
        filename = ev_and_pp_type + "_" + matrix_type + "_matrices.txt"
        subdir = "../output/" + ev_and_pp_type + "_" + case_type + "_"+ self.get_date_str() + "/text_files/"
        os.makedirs(subdir, exist_ok = True)
        f = open(subdir + filename, "a")
        f.write(output)
        f.close()

    def format_fl_vals(self, val):
        return "{:.5f}".format(val)


    def plot_graph_comp(self, total_avg_param_err_list):
        fig, ax = self.get_subplots_comp()
        ax.plot(list(range(len(total_avg_param_err_list))), total_avg_param_err_list, label = "MX" + " " + str(self.get_static_vars_dict_elt("loop_limit")))

    def display_avg_param_err_indv(self, t, static_args, true_params_title, avg_param_err, curr_indx, linthresh_value, matrix_type, error_type):
        err_label_1 = err_label_2 = ""

        if(error_type == "Absolute"):
            err_label_1, err_label_2 = "Abs.", "abs"
            line_color = "brown"

        else:
            err_label_1, err_label_2 = "Rel.", "rel"
            line_color = "green"

        ev_type, pp_type, bounds, loop_limit,  param_labels, case_type = static_args
        param_label_1, param_label_2 = param_labels

        fig, ax = plt.subplots()
        fig.set_size_inches(16, 8)
        ax.plot(t, avg_param_err, color = line_color, label = "Avg. " + err_label_1 + " Err: " + "\n"+ param_label_1 + " & " + param_label_2 )
        ax.set_xlabel("Time")
        ax.set_ylabel("Avg. " + error_type + " Parameter Error - " + param_label_1 + " and " +  param_label_2)
        ax.set_xscale("log")
        # ax.set_yscale("log")
        ax.set_yscale("symlog", linthresh = linthresh_value)

        title = ev_type.upper()+ " | " + pp_type.upper() + " | " + case_type  + " | BNDS " + bounds + " | MTRX " + curr_indx + " : " + true_params_title
        ax.set_title(label = title, pad = 20, fontsize = 15)

        # Output image files
        output_type = "indv"
        ev_and_pp_type = ev_type + "_" + pp_type
        subdir = "../output/" + ev_and_pp_type + "_" + case_type + "_"+ self.get_date_str() + "/avg_param_err_graphs" + "/avg_" + err_label_2 + "_param_err_graphs/avg_"+ err_label_2 +"_param_err_graph" + "_" + output_type + "/"
        # subdir = "../output/" + ev_and_pp_type + "_" + case_type + "_"+ self.get_date_str() + "/avg_" + err_label_2 + "_param_err_graphs/avg_"+ err_label_2 +"_param_err_graph" + "_" + output_type + "/"
        filename = "avg_" + err_label_2 + "_param_err_graph" + "_" + ev_and_pp_type
        ax.legend(loc = 'center left', bbox_to_anchor = (1, 0.5))
        os.makedirs(subdir, exist_ok = True)
        fig.savefig(subdir + filename + "_mtrx_" + curr_indx + ".png", dpi = 300)
        plt.close('all')


    def display_param_err(self, t, static_args, param_1_err, param_2_err, true_params_title, curr_indx,  linthresh_value, matrix_type,  error_type):
        ev_type, pp_type, bounds, loop_limit,  param_labels, case_type = static_args
        param_label_1, param_label_2 = param_labels
        err_label_1 = err_label_2 = ""
        if(error_type == "Absolute"):
            err_label =  "abs"
        else:
            err_label = "rel"

        fig, ax = plt.subplots()
        fig.set_size_inches(16, 8)
        ax.set_xlabel("Time")
        ax.set_ylabel(param_label_1 + " and " + param_label_2 + " - " + error_type + " Parameter Error")
        ax.set_xscale("log")
        # ax.set_yscale("log")
        ax.set_yscale("symlog", linthresh = linthresh_value)
        ax.plot(t, param_1_err, label = param_label_1 + " " + err_label + ". err." )
        ax.plot(t, param_2_err, label = param_label_2 + " " + err_label + ". err." )

        # For title
        title = ev_type.upper()+ " | " + pp_type.upper() + " | " + case_type + " | BNDS " + bounds + " | MTRX " + curr_indx + " : " + true_params_title
        ax.set_title(label = title, pad = 20, fontsize = 15)

        # Output img files
        ev_and_pp_type = ev_type + "_" + pp_type
        subdir = "../output/" + ev_and_pp_type + "_" + case_type + "_" + self.get_date_str() + "/param_err_graphs" + "/"+ err_label + "_param_err_graphs/" +  err_label + "_param_err_graph" + "_" + matrix_type + "/"
        # subdir = "../output/" + ev_and_pp_type + "_" + case_type + "_"+ self.get_date_str() + "/"+ err_label + "_param_err_graphs/" +  err_label + "_param_err_graph" + "_" + matrix_type + "/"
        filename = err_label + "_param_err_graph" + "_" + ev_and_pp_type
        ax.legend(loc = 'center left', bbox_to_anchor = (1, 0.5))
        os.makedirs(subdir, exist_ok = True)
        fig.savefig(subdir + filename + "_mtrx_" + curr_indx + ".png", dpi = 300)
        plt.close('all')



    # def display_avg_param_err_comp(self, error_type):
    #     if self.has_bad_matrices_comp() != False:
    #         err_label_1 = err_label_2 = ""
    #         if(error_type == "Absolute"):
    #             err_label = "abs"
    #         else:
    #             err_label =  "rel"
    #
    #         true_params, ev_type, pp_type = self.get_true_params_comp(), self.get_static_vars_dict_elt("ev_type"), self.get_static_vars_dict_elt("pp_type")
    #         fig, ax = self.get_subplots_comp()
    #         ax.set_xscale("log")
    #         ax.set_yscale("log")
    #         ax.legend(loc = 'center left', bbox_to_anchor = (1, 0.5))
    #
    #         # Output files
    #         output_type = "comp"
    #         ev_and_pp_type = ev_type + '_' + pp_type
    #         subdir =   "../output/" + ev_and_pp_type + "_" + case_type + "_"+ self.get_date_str() + "/"
    #         filename =  "avg_rel" + err_label +"_param_err_graph" + "_" + ev_and_pp_type + "_" + output_type
    #         os.makedirs(subdir, exist_ok = True)
    #         fig.savefig(subdir + filename + ".png", dpi = 300)
    #         plt.close(fig)
    #
    #     else:
    #         fig, ax = self.get_subplots_comp()
    #         plt.close(fig)
    #         print("RESULT: No bad matrices were found.")



    def display_sol_xy(self, static_args, S_lists, true_params_title, curr_indx, linthresh_value, matrix_type):
        ev_type, pp_type, bounds, loop_limit,  param_labels, case_type = static_args
        x_estimates, y_estimates, xt_estimates, yt_estimates = S_lists
        fig, (ax_1, ax_2) = plt.subplots(1, 2)
        fig.set_size_inches(30, 15)
        label_font_size = 18

        ax_1.set_xlabel("X", fontsize = label_font_size)
        ax_1.set_ylabel("Y", fontsize = label_font_size)
        ax_1.plot(x_estimates, y_estimates, label = "True Sol." , color = "green")
        ax_1.set_xscale("symlog", linthresh = linthresh_value)
        ax_1.set_yscale("symlog", linthresh = linthresh_value)
        ax_1.legend(loc = "upper right")


        ax_2.set_xlabel("XT", fontsize = label_font_size)
        ax_2.set_ylabel("YT", fontsize = label_font_size)
        ax_2.plot(xt_estimates, yt_estimates, label = "Est. Sol.", color = "red" )
        ax_2.set_xscale("symlog", linthresh = linthresh_value)
        ax_2.set_yscale("symlog", linthresh = linthresh_value)
        ax_2.legend(loc = "upper right")

        # self.animation_output_1x2(x_estimates, y_estimates, xt_estimates, yt_estimates, linthresh_value)


        # For file title
        title = ev_type.upper()+ " | " + pp_type.upper() + " | " + case_type + " | BNDS " + bounds + " | MTRX " + curr_indx + " : " + true_params_title
        # ax.set_title(label = title, pad = 20, fontsize = 15)
        fig.suptitle(title, fontsize = label_font_size)

        # Output IMG files
        ev_and_pp_type = ev_type + "_" + pp_type
        subdir = "../output/" + ev_and_pp_type + "_" + case_type + "_"+ self.get_date_str() + "/sol_xy_graphs/sol_xy_graph" + "_" + matrix_type + "/"
        filename = "sol_xy_graph" + "_" + ev_and_pp_type

        os.makedirs(subdir, exist_ok = True)
        fig.savefig(subdir + filename + "_mtrx_" + curr_indx + ".png", dpi = 300)
        plt.close('all')



    def display_sol_xy_over_t(self, t, static_args, S_lists, true_params_title, curr_indx, linthresh_value, matrix_type):
        ev_type, pp_type, bounds, loop_limit,  param_labels, case_type = static_args
        x_estimates, y_estimates, xt_estimates, yt_estimates = S_lists
        fig, (ax_1, ax_2) = plt.subplots(2, 1)
        fig.set_size_inches(16, 12)

        # For x and xt
        ax_1.set_ylabel("X and XT")
        ax_1.set_xlabel("Time")
        ax_1.plot(t, x_estimates, label = "X")
        ax_1.plot(t, xt_estimates, label = "XT")
        ax_1.legend(loc = "best")
        # ax_1.set_xscale("symlog", linthresh = linthresh_value)
        ax_1.set_xscale("log")
        ax_1.set_yscale("symlog", linthresh = linthresh_value)

        # For y and yt
        ax_2.set_ylabel("Y and YT")
        ax_2.set_xlabel("Time")
        ax_2.plot(t, y_estimates, label = "Y")
        ax_2.plot(t, yt_estimates, label = "YT")
        ax_2.legend(loc = "best")

        ax_2.set_xscale("log")
        ax_2.set_yscale("symlog", linthresh = linthresh_value)
        # ax_2.set_xscale("symlog", linthresh = linthresh_value)

        #For title
        title = ev_type.upper()+ " | " + pp_type.upper() + " | " + case_type + " | BNDS " + bounds + " | MTRX " + curr_indx + " : " + true_params_title
        fig.suptitle(title, fontsize = 15)

        # Output IMG files
        ev_and_pp_type = ev_type + "_" + pp_type
        subdir = "../output/" + ev_and_pp_type + "_" + case_type + "_"+ self.get_date_str() + "/sol_xy_over_t_graphs/sol_xy_over_t_graph" + "_" + matrix_type + "/"
        filename = "sol_xy_over_t_graph" + "_" + ev_and_pp_type
        os.makedirs(subdir, exist_ok = True)
        fig.savefig(subdir + filename + "_mtrx_" + curr_indx + ".png", dpi = 300)
        plt.close('all')


    def display_sol_x(self, t, static_args, S_lists, true_params_title, curr_indx, linthresh_value, matrix_type):
        ev_type, pp_type, bounds, loop_limit,  param_labels, case_type = static_args
        x_estimates, y_estimates, xt_estimates, yt_estimates = S_lists
        fig, ax = plt.subplots()
        fig.set_size_inches(16, 8)

        # For x and xt
        ax.set_ylabel("X and XT")
        ax.set_xlabel("Time")
        ax.plot(t, x_estimates, label = 'X')
        ax.plot(t, xt_estimates, label = "XT")
        ax.set_xscale('log')
        ax.set_yscale('symlog', linthresh = linthresh_value)
        ax.legend(loc = "best")

        title = ev_type.upper()+ " | " + pp_type.upper() + " | " + case_type + " | BNDS " + bounds + " | MTRX " + curr_indx + " : " + true_params_title
        fig.suptitle(title, fontsize = 15)

        # Output IMG files
        ev_and_pp_type = ev_type + "_" + pp_type
        subdir = "../output/" + ev_and_pp_type + "_" + case_type + "_"+ self.get_date_str() + "/sol_x_graphs/sol_x_graph" + "_" + matrix_type + "/"
        filename = "sol_x_graph" + "_" + ev_and_pp_type
        os.makedirs(subdir, exist_ok = True)
        fig.savefig(subdir + filename + "_mtrx_" + curr_indx + ".png", dpi = 300)
        plt.close('all')


    def display_sol_y(self, t, static_args, S_lists, true_params_title, curr_indx, linthresh_value, matrix_type):
        ev_type, pp_type, bounds, loop_limit,  param_labels, case_type = static_args
        x_estimates, y_estimates, xt_estimates, yt_estimates = S_lists
        fig, ax = plt.subplots()
        fig.set_size_inches(16, 8)
        # https://matplotlib.org/2.0.2/api/pyplot_api.html

        # For x and xt
        ax.set_ylabel("Y and YT")
        ax.set_xlabel("Time")
        ax.plot(t, y_estimates, label = 'Y')
        ax.plot(t, yt_estimates, label = "YT")
        ax.set_xscale('log')
        ax.set_yscale('symlog', linthresh = linthresh_value)
        ax.legend(loc = "best")

        title = ev_type.upper()+ " | " + pp_type.upper() + " | " + case_type + " | BNDS " + bounds + " | MTRX " + curr_indx + " : " + true_params_title
        fig.suptitle(title, fontsize = 15)

        # Output IMG files
        ev_and_pp_type = ev_type + "_" + pp_type
        subdir = "../output/" + ev_and_pp_type + "_" + case_type + "_"+ self.get_date_str() + "/sol_y_graphs/sol_y_graph" + "_" + matrix_type + "/"
        filename = "sol_y_graph" + "_" + ev_and_pp_type
        os.makedirs(subdir, exist_ok = True)
        fig.savefig(subdir + filename + "_mtrx_" + curr_indx + ".png", dpi = 300)
        plt.close('all')


    def display_sol_signal_err(self, t, static_args, U_lists, true_params_title, curr_indx, linthresh_value, matrix_type):
        ev_type, pp_type, bounds, loop_limit,  param_labels, case_type = static_args
        x_sig_err_estimates, y_sig_err_estimates = U_lists
        fig, ax = plt.subplots()
        fig.set_size_inches(16, 8)
        ax.set_ylabel("Signal Error")
        ax.set_xlabel("Time")
        ax.plot(t, x_sig_err_estimates, label = "X Signal Err.")
        ax.plot(t, y_sig_err_estimates, label = "Y Signal Err.")
        ax.legend(loc = "best")
        ax.set_xscale("log")
        ax.set_yscale("symlog", linthresh = linthresh_value) #signal err

        #For title
        title = ev_type.upper()+ " | " + pp_type.upper() + " | " + case_type + " | BNDS " + bounds + " | MTRX " + curr_indx + " : " + true_params_title
        fig.suptitle(title, fontsize = 15)

        # Output IMG files
        ev_and_pp_type = ev_type + "_" + pp_type
        subdir = "../output/" + ev_and_pp_type + "_" + case_type + "_"+ self.get_date_str() + "/sol_signal_err_graphs/sol_signal_err_graph" + "_" + matrix_type + "/"
        filename = "sol_signal_err_graph" + "_" + ev_and_pp_type
        os.makedirs(subdir, exist_ok = True)
        fig.savefig(subdir + filename + "_mtrx_" + curr_indx + ".png", dpi = 300)
        plt.close('all')



    def display_sol_signal_err_split(self, t, static_args, U_lists, true_params_title, curr_indx, linthresh_value, matrix_type):
        ev_type, pp_type, bounds, loop_limit,  param_labels, case_type = static_args
        x_sig_err_estimates, y_sig_err_estimates = U_lists
        fig, (ax_1, ax_2) = plt.subplots(2, 1)
        fig.set_size_inches(16, 12)

        ax_1.set_ylabel("Signal Error")
        ax_1.set_xlabel("Time")
        ax_1.plot(t, x_sig_err_estimates, label = "X Signal Err.")

        ax_1.legend(loc = "best")
        ax_1.set_xscale("log")
        ax_1.set_yscale("symlog", linthresh = linthresh_value) #signal err

        ax_2.set_ylabel("Signal Error")
        ax_2.set_xlabel("Time")
        ax_2.plot(t, y_sig_err_estimates, label = "Y Signal Err.", color = "orange")
        ax_2.legend(loc = "best")
        ax_2.set_xscale("log")
        ax_2.set_yscale("symlog", linthresh = linthresh_value) #signal err

        #For title
        title = ev_type.upper()+ " | " + pp_type.upper() + " | " + case_type + " | BNDS " + bounds + " | MTRX " + curr_indx + " : " + true_params_title
        fig.suptitle(title, fontsize = 15)

        # Output IMG files
        ev_and_pp_type = ev_type + "_" + pp_type
        subdir = "../output/" + ev_and_pp_type + "_" + case_type + "_"+ self.get_date_str() + "/sol_signal_err_split_graphs/sol_signal_err_split_graph" + "_" + matrix_type + "/"
        filename = "sol_signal_err_split_graph" + "_" + ev_and_pp_type
        os.makedirs(subdir, exist_ok = True)
        fig.savefig(subdir + filename + "_mtrx_" + curr_indx + ".png", dpi = 300)
        plt.close('all')


#-------------- GETTERS ---------------------------------------------


    def get_date_str(self):
        raw_date = datetime.datetime.now()
        t = time.localtime()
        date_str = str(raw_date.year) + "."+ str(time.strftime("%m"))+ "." + str(time.strftime("%d"))
        return date_str

    def get_subplots_comp(self):
        return self._fig, self._ax

    def get_true_params_comp(self):
        return self._true_params

    def has_bad_matrices_comp(self):
        return self._has_bad_matrices

    def get_static_vars_dict(self):
        return self._static_vars_dict

    def get_static_vars_dict_elt(self, key):
        return self._static_vars_dict[key]


#-------------- SETTERS -------------------------------------
    def set_has_bad_matrices_comp(self, val):
        self._has_bad_matrices = val

    def set_ev_type(self, val):
        self._ev_type = val

    def set_pp_type(self, val):
        self._pp_type = val

    def set_static_vars_dict_elt(self, key, value):
        self._static_vars_dict[key] = value

    def set_subplots_comp(self):
        true_params             = self.get_true_params_comp()
        case_type   = self.get_static_vars_dict_elt("case_type")
        ev_type               = self.get_static_vars_dict_elt("ev_type")
        pp_type               = self.get_static_vars_dict_elt("pp_type")
        bounds                = self.get_static_vars_dict_elt("bounds")
        self._fig, self._ax = plt.subplots()
        title = ev_type.upper() + " | " + pp_type.upper() + " | " + case_type + " | BNDS: " + bounds + " | Matrices w/ an Avg. Err. > 1e-5"
        self._ax.set_xscale('log') # Iterations
        self._ax.set_yscale('log') # Avg. Rel. Err
        self._ax.set_title(label = title, pad = 25, fontsize = 15)
        self._ax.set_xlabel("Iterations")
        self._ax.set_ylabel("Avg. Error")
        self._fig.set_size_inches(16, 8)


    def get_true_params_title(self, true_params):
        a11, a12, a21, a22 = true_params
        a11_str, a12_str = self.format_fl_vals(a11), self.format_fl_vals(a12)
        a21_str, a22_str = self.format_fl_vals(a21), self.format_fl_vals(a22)
        return "[" + a11_str + ", " + a12_str + ", " + a21_str + ", " + a22_str + "]"
