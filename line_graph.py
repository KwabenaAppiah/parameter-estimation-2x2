import matplotlib.pyplot as plt
import numpy as np
import datetime
import time
import math
import os


class LineGraph:
    def __init__(self, *args):
        ev_type, pp_type, integration_method, bounds, loop_limit, case_type = args

        # FOR COMP ONLY
        self._true_values = []
        self._has_bad_matrices = False
        self._static_vars_dict = {
            'ev_type': ev_type, 'pp_type': pp_type, 'integration_method': integration_method, 'bounds': bounds,  'loop_limit': loop_limit, 'case_type': case_type
        }
        self.set_subplots_comp()


    def get_date_str(self):
        raw_date = datetime.datetime.now()
        t = time.localtime()
        # current_time = time.strftime("%H.%M", t)
        # date_str = str(raw_date.year) + "." + str(time.strftime('%m')) + "." + str(time.strftime('%d')) + current_time
        date_str = str(raw_date.year) + "."+ str(time.strftime("%m"))+ "." + str(time.strftime("%d"))
        return date_str

    def get_subplots_comp(self):
        return self._fig, self._ax

    def get_true_values_comp(self):
        return self._true_values

    def has_bad_matrices_comp(self):
        return self._has_bad_matrices

    def get_static_vars_dict(self):
        return self._static_vars_dict

    def get_static_vars_dict_elt(self, key):
        return self._static_vars_dict[key]


    # COMP SETTERS
    def set_has_bad_matrices_comp(self, val):
        self._has_bad_matrices = val

    def set_ev_type(self, val):
        self._ev_type = val

    def set_pp_type(self, val):
        self._pp_type = val

    def set_static_vars_dict_elt(self, key, value):
        self._static_vars_dict[key] = value

    def set_subplots_comp(self):
        true_values             = self.get_true_values_comp()
        integration_method    = self.get_static_vars_dict_elt("integration_method")
        ev_type               = self.get_static_vars_dict_elt("ev_type")
        pp_type               = self.get_static_vars_dict_elt("pp_type")
        bounds                = self.get_static_vars_dict_elt("bounds")
        self._fig, self._ax = plt.subplots()
        title = ev_type.upper() + " | " + pp_type.upper() + " | " + integration_method + " | BNDS: " + bounds + " | Matrices w/ an Avg. Err. > 1e-5"
        self._ax.set_xscale('log') # Iterations
        self._ax.set_yscale('log') # Avg. Rel. Err
        self._ax.set_title(label = title, pad = 25, fontsize = 15)
        self._ax.set_xlabel("Iterations")
        self._ax.set_ylabel("Avg. Error")
        self._fig.set_size_inches(16, 8)


    def get_true_values_title(self, true_values):
        a11, a12, a21, a22 = true_values
        a11_str, a12_str = self.format_fl_vals(a11), self.format_fl_vals(a12)
        a21_str, a22_str = self.format_fl_vals(a21), self.format_fl_vals(a22)
        return "[" + a11_str + ", " + a12_str + ", " + a21_str +", " + a22_str + "]"


    def get_param_error(self, param_1, param_2, param_1_list, param_2_list, error_type):
        total_avg_param_err_list, param_1_err_list, param_2_err_list = [], [], []

        for i in range(len(param_1_list)):
            # Relative Error
            if(error_type == "Relative"):
                param_1_err_list.append(abs(param_1_list[i] - param_1) / abs(param_1))
                param_2_err_list.append(abs(param_2_list[i] - param_2) / abs(param_2))

            # Absolute Error
            elif(error_type == "Absolute"):
                param_1_err_list.append(abs(param_1_list[i] - param_1))
                param_2_err_list.append(abs(param_2_list[i] - param_2))

            else:
                print(type + "is not a valid error type.")
                quit()

            total_avg_param_err_list.append(abs(param_1_err_list[i] + param_2_err_list[i]) / 2)

        nth_avg_param_err = abs(param_1_err_list[-1] + param_2_err_list[-1]) / 2

        return param_1_err_list, param_2_err_list, total_avg_param_err_list, nth_avg_param_err


    def organize_data(self, static_args, sol_time, time_events, guesses, true_values, param_labels, curr_idx, threshold, solutions, signal_err):
        case_type = self.get_static_vars_dict_elt('case_type')
        pp_type   = self.get_static_vars_dict_elt('pp_type')
        curr_idx = str(curr_idx)
        a11s, a12s, a21s, a22s = guesses
        a11, a12, a21, a22     = true_values
        param_1_rel_err_list, param_2_rel_err_list, total_avg_rel_param_err_list, nth_avg_rel_param_err = [], [], [], 0
        param_1_abs_err_list, param_2_abs_err_list, total_avg_abs_param_err_list, nth_avg_abs_param_err = [], [], [], 0
        true_values_str   = str(true_values)
        true_values_title = self.get_true_values_title(true_values)
        was_content_plotted = False
        param_1 = param_2 = 0
        param_l_list, param_2_list = [], []
        matrix_type = ""
        label_1, label_2 = param_labels

        if case_type == "MAIN_DIAG_2x2":
            param_1, param_2 = a11, a22
            param_1_list, param_2_list = a11s, a22s


        elif(case_type == "OFF_DIAG_2x2"):
            print("Off Diagonal 2x2 - This option is not available yet.")
            exit()

        elif(case_type == "BTM_ROW_2x2"):
            print("Bottom Row 2x2 - This option is not available yet.")
            exit()

        if param_1 != 0 and param_2 != 0:
            if(pp_type == "sink" or pp_type == "sp_sink"):
                linthresh_value = .000000000000001
            else:
                linthresh_value = .1

            param_1_rel_err_list, param_2_rel_err_list, total_avg_rel_param_err_list, nth_avg_rel_param_err = self.get_param_error(param_1, param_2, param_1_list, param_2_list, "Relative")
            param_1_abs_err_list, param_2_abs_err_list, total_avg_abs_param_err_list, nth_avg_abs_param_err = self.get_param_error(param_1, param_2, param_1_list, param_2_list, "Absolute")

            if nth_avg_rel_param_err >= threshold:
                matrix_type = "bad"
                self.plot_graph_comp(total_avg_rel_param_err_list)
                self.set_has_bad_matrices_comp(True)
                data_was_plotted = True

            elif nth_avg_rel_param_err < threshold:
                matrix_type = "good"
                data_was_plotted = True

            else:
                data_was_plotted = False

            # # Output regardless
            self.write_to_file("mtrx_" + curr_idx + "|" + true_values_str, "all")
            self.write_to_file("mtrx_" + curr_idx + "|" + true_values_str, matrix_type)
            #
            # # Parameter Error-based graphs
            self.display_avg_param_err_indv(static_args, true_values_title, param_labels, total_avg_rel_param_err_list, curr_idx, "Relative")
            self.display_param_err(static_args, sol_time, param_1_rel_err_list, param_2_rel_err_list, true_values_title, param_labels, curr_idx, matrix_type, "Relative")

            self.display_avg_param_err_indv(static_args, true_values_title, param_labels, total_avg_rel_param_err_list, curr_idx, "Absolute")
            self.display_param_err(static_args, sol_time, param_1_rel_err_list, param_2_rel_err_list, true_values_title, param_labels, curr_idx, matrix_type, "Absolute")

            # # Solutions-based graphs
            self.display_sol_all(static_args, time_events, solutions, true_values_title, curr_idx, linthresh_value, matrix_type)
            self.display_sol_x(static_args, time_events, solutions, true_values_title, curr_idx, linthresh_value, matrix_type)
            self.display_sol_y(static_args, time_events, solutions, true_values_title, curr_idx, linthresh_value,matrix_type)
            self.display_sol_xy(static_args, solutions, true_values_title, curr_idx, linthresh_value, matrix_type)
            self.display_signal_err(static_args, time_events, signal_err, true_values_title, curr_idx, linthresh_value, matrix_type)

        return was_content_plotted



    def plot_graph_comp(self, total_avg_param_err_list):
        fig, ax = self.get_subplots_comp()
        ax.plot(list(range(len(total_avg_param_err_list))), total_avg_param_err_list, label = "MX" + " " + str(self.get_static_vars_dict_elt("loop_limit")))


    def write_to_file(self, line_str, mtrx_type):
        output = line_str + "\n" + "\n"
        ev_and_pp_type = self.get_static_vars_dict_elt("ev_type") + "_" + self.get_static_vars_dict_elt("pp_type")
        filename = ev_and_pp_type + "_" + mtrx_type + "_matrices.txt"
        subdir = "../output/" + ev_and_pp_type + "_" + self.get_date_str() + "/txt_files/"
        os.makedirs(subdir, exist_ok = True)
        f = open(subdir + filename, "a")
        f.write(output)
        f.close()

    def format_fl_vals(self, val):
        return "{:.5f}".format(val)


    def display_avg_param_err_indv(self, static_args, true_values_title, param_labels, total_avg_param_err_list, curr_idx, error_type):
        err_label_1 = err_label_2 = ""
        if(error_type == "Absolute"):
            err_label_1, err_label_2 = "Abs.", "abs"
        else:
            err_label_1, err_label_2 = "Rel.", "rel"

        param_label_1, param_label_2 = param_labels
        ev_type, pp_type, integration_method, bounds, loop_limit, case_type  = static_args
        fig, ax = plt.subplots()
        fig.set_size_inches(16, 8)
        ax.plot(list(range(len(total_avg_param_err_list))), total_avg_param_err_list, color = "green", label = "Avg. " + err_label_1 + " Err: " + param_label_1 + " & " + param_label_2 )
        ax.set_xlabel("Iterations")
        ax.set_ylabel("Avg. " + error_type + " Parameter Error - " + param_label_1 + " and " +  param_label_2)
        ax.set_yscale("log")
        ax.set_xscale("log")

        title = ev_type.upper()+ " | " + pp_type.upper() + " | " + integration_method + " | BNDS " + bounds + " | MTRX " + curr_idx + " : " + true_values_title
        ax.set_title(label = title, pad = 20, fontsize = 15)

        # Output IMG files
        output_type = "indv"
        ev_and_pp_type = ev_type + "_" + pp_type
        subdir = "../output/" + ev_and_pp_type + "_" + self.get_date_str() + "/avg_" + err_label_2 + "_param_err_graphs/avg_"+ err_label_2 +"_param_err_graph" + "_" + output_type + "/"
        filename = "avg_" + err_label_2 + "_param_err_graph" + "_" + ev_and_pp_type
        ax.legend(loc = "upper right")
        os.makedirs(subdir, exist_ok = True)
        fig.savefig(subdir + filename + "_mtrx_" + curr_idx + ".png", dpi = 300)
        plt.close('all')


    def display_avg_param_err_comp(self, error_type):
        if self.has_bad_matrices_comp() != False:
            err_label_1 = err_label_2 = ""
            if(error_type == "Absolute"):
                err_label = "abs"
            else:
                err_label =  "rel"

            true_values, ev_type, pp_type = self.get_true_values_comp(), self.get_static_vars_dict_elt("ev_type"), self.get_static_vars_dict_elt("pp_type")
            fig, ax = self.get_subplots_comp()
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.legend(loc = 'center left', bbox_to_anchor = (1, 0.5))

            # Output files
            output_type = "comp"
            ev_and_pp_type = ev_type + '_' + pp_type
            subdir =   "../output/" + ev_and_pp_type + "_" + self.get_date_str() + "/"
            filename =  "avg_rel" + err_label +"_param_err_graph" + "_" + ev_and_pp_type + "_" + output_type
            os.makedirs(subdir, exist_ok = True)
            fig.savefig(subdir + filename + ".png", dpi = 300)
            plt.close(fig)

        else:
            fig, ax = self.get_subplots_comp()
            plt.close(fig)
            print("RESULT: No bad matrices were found.")


    def display_param_err(self, static_args, sol_time, param_1_err, param_2_err, true_values_title, param_labels, curr_idx,  mtrx_type, error_type):
        ev_type, pp_type, integration_method, bounds, loop_limit, case_type  = static_args
        param_label_1, param_label_2 = param_labels
        err_label_1 = err_label_2 = ""
        if(error_type == "Absolute"):
            err_label =  "abs"
        else:
            err_label = "rel"
        # time = sol.t
        time = sol_time
        fig, ax = plt.subplots()
        fig.set_size_inches(16, 8)
        ax.set_xlabel("Time")
        ax.set_ylabel(param_label_1 + " and " + param_label_2 + " - " + error_type + " Parameter Error")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.plot(time, param_1_err, label = param_label_1 + " " + err_label + ". err" )
        ax.plot(time, param_2_err, label = param_label_2 + " " + err_label + ". err" )

        # For title
        title = ev_type.upper()+ " | " + pp_type.upper() + " | " + integration_method + " | BNDS " + bounds + " | MTRX " + curr_idx + " : " + true_values_title
        ax.set_title(label = title, pad = 20, fontsize = 15)

        # Output img files
        ev_and_pp_type = ev_type + "_" + pp_type
        subdir = "../output/" + ev_and_pp_type + "_" + self.get_date_str() + "/"+ err_label + "_param_err_graphs/" +  err_label + "_param_err_graph" + "_" + mtrx_type + "/"
        filename = err_label + "_param_err_graph" + "_" + ev_and_pp_type
        ax.legend(loc = "upper right")
        os.makedirs(subdir, exist_ok = True)
        fig.savefig(subdir + filename + "_mtrx_" + curr_idx + ".png", dpi = 300)
        plt.close('all')


    def display_sol_xy(self, static_args, solutions, true_values_title, curr_idx, linthresh_value, mtrx_type):
        ev_type, pp_type, integration_method, bounds, loop_limit, case_type  = static_args
        x, y, xt, yt = solutions
        fig, ax = plt.subplots()
        fig.set_size_inches(16, 8)
        ax.set_xlabel("X and XT")
        ax.set_ylabel("Y and YT")
        ax.plot(x, y, label = "True Sol" , color = "green")
        ax.plot(xt, yt, label = "Est. Sol", color = "red" )
        print("x:", x[-1])
        print("y:", y[-1], "\n")
        print("xt", xt[-1])
        print("yt:", yt[-1], "\n")
        ax.set_xscale("symlog", linthresh = linthresh_value)
        ax.set_yscale("symlog", linthresh = linthresh_value)

        # For file title
        title = ev_type.upper()+ " | " + pp_type.upper() + " | " + integration_method + " | BNDS " + bounds + " | MTRX " + curr_idx + " : " + true_values_title
        ax.set_title(label = title, pad = 20, fontsize = 15)

        # Output IMG files
        ev_and_pp_type = ev_type + "_" + pp_type
        subdir = "../output/" + ev_and_pp_type + "_" + self.get_date_str() + "/sol_graphs_xy/sol_graph_xy" + "_" + mtrx_type + "/"
        filename = "sol_graph_all" + "_" + ev_and_pp_type
        ax.legend(loc = "best", bbox_to_anchor = (1, 0.5))
        os.makedirs(subdir, exist_ok = True)
        fig.savefig(subdir + filename + "_mtrx_" + curr_idx + ".png", dpi = 300)
        plt.close('all')



    def display_sol_all(self, static_args, time_events, solutions, true_values_title, curr_idx, linthresh_value, mtrx_type):
        ev_type, pp_type, integration_method, bounds, loop_limit, case_type  = static_args
        x, y, xt, yt = solutions
        fig, (ax_1, ax_2) = plt.subplots(2, 1)
        fig.set_size_inches(16, 8)

        # For x and xt
        ax_1.set_ylabel("X and XT")
        ax_1.set_xlabel("Time")
        ax_1.plot(time_events[0], x, label = "X")
        ax_1.plot(time_events[0], xt, label = "XT")
        ax_1.legend(loc = "best")
        ax_1.set_yscale("symlog", linthresh = linthresh_value) # x and xt

        # For y and yt
        ax_2.set_ylabel("Y and YT")
        ax_2.set_xlabel("Time")
        ax_2.plot(time_events[0], y, label = "Y")
        ax_2.plot(time_events[0], yt, label = "YT")
        ax_2.legend(loc = "best")
        ax_2.set_yscale("symlog", linthresh = linthresh_value) # y and yt

        #For title
        title = ev_type.upper()+ " | " + pp_type.upper() + " | " + integration_method + " | BNDS " + bounds + " | MTRX " + curr_idx + " : " + true_values_title
        fig.suptitle(title, fontsize = 15)

        # Output IMG files
        ev_and_pp_type = ev_type + "_" + pp_type
        subdir = "../output/" + ev_and_pp_type + "_" + self.get_date_str() + "/sol_graph_all/sol_graph_all" + "_" + mtrx_type + "/"
        filename = "sol_graph_all" + "_" + ev_and_pp_type
        os.makedirs(subdir, exist_ok = True)
        fig.savefig(subdir + filename + "_mtrx_" + curr_idx + ".png", dpi = 300)
        plt.close('all')


    def display_sol_x(self, static_args, time_events, solutions, true_values_title, curr_idx, linthresh_value, mtrx_type):
        ev_type, pp_type, integration_method, bounds, loop_limit, case_type  = static_args
        x, y, xt, yt = solutions
        fig, (ax_1) = plt.subplots(1, 1)
        fig.set_size_inches(16, 8)

        # For x and xt
        ax_1.set_ylabel("X and XT")
        ax_1.set_xlabel("Time")
        ax_1.plot(time_events[0], x, label = 'X')
        ax_1.plot(time_events[0], xt, label = "XT")
        plt.yscale('symlog', linthresh = linthresh_value)
        ax_1.legend(loc = "best")

        title = ev_type.upper()+ " | " + pp_type.upper() + " | " + integration_method + " | BNDS " + bounds + " | MTRX " + curr_idx + " : " + true_values_title
        fig.suptitle(title, fontsize = 15)

        # Output IMG files
        ev_and_pp_type = ev_type + "_" + pp_type
        subdir = "../output/" + ev_and_pp_type + "_" + self.get_date_str() + "/sol_graph_x/sol_graph_x" + "_" + mtrx_type + "/"
        filename = "sol_graph_x" + "_" + ev_and_pp_type
        os.makedirs(subdir, exist_ok = True)
        fig.savefig(subdir + filename + "_mtrx_" + curr_idx + ".png", dpi = 300)
        plt.close('all')


    def display_sol_y(self, static_args, time_events, solutions, true_values_title, curr_idx, linthresh_value, mtrx_type):
        ev_type, pp_type, integration_method, bounds, loop_limit, case_type  = static_args
        x, y, xt, yt = solutions
        fig, (ax_1) = plt.subplots(1, 1)
        fig.set_size_inches(16, 8)
        # https://matplotlib.org/2.0.2/api/pyplot_api.html

        # For x and xt
        ax_1.set_ylabel("Y and YT")
        ax_1.set_xlabel("Time")
        ax_1.plot(time_events[0], y, label = 'Y')
        ax_1.plot(time_events[0], yt, label = "YT")
        plt.yscale('symlog', linthresh = linthresh_value)

        ax_1.legend(loc = "best")

        title = ev_type.upper()+ " | " + pp_type.upper() + " | " + integration_method + " | BNDS " + bounds + " | MTRX " + curr_idx + " : " + true_values_title
        fig.suptitle(title, fontsize = 15)

        # Output IMG files
        ev_and_pp_type = ev_type + "_" + pp_type
        subdir = "../output/" + ev_and_pp_type + "_" + self.get_date_str() + "/sol_graph_y/sol_graph_y" + "_" + mtrx_type + "/"
        filename = "sol_graph_y" + "_" + ev_and_pp_type
        os.makedirs(subdir, exist_ok = True)
        fig.savefig(subdir + filename + "_mtrx_" + curr_idx + ".png", dpi = 300)
        plt.close('all')


    def display_signal_err(self, static_args, time_events, signal_err, true_values_title, curr_idx, linthresh_value, mtrx_type):
        ev_type, pp_type, integration_method, bounds, loop_limit, case_type  = static_args
        x_signal_err_list, y_signal_err_list = signal_err
        fig, ax = plt.subplots()
        fig.set_size_inches(16, 8)
        ax.set_ylabel("Signal Error")
        ax.set_xlabel("Time")
        ax.plot(time_events[0], x_signal_err_list, label = "X Signal Err")
        ax.plot(time_events[0], y_signal_err_list, label = "Y Signal Err")
        ax.legend(loc = "best")
        ax.set_yscale("symlog", linthresh = linthresh_value) #signal err

        #For title
        title = ev_type.upper()+ " | " + pp_type.upper() + " | " + integration_method + " | BNDS " + bounds + " | MTRX " + curr_idx + " : " + true_values_title
        fig.suptitle(title, fontsize = 15)

        # Output IMG files
        ev_and_pp_type = ev_type + "_" + pp_type
        subdir = "../output/" + ev_and_pp_type + "_" + self.get_date_str() + "/signal_err_graph/signal_err_graph" + "_" + mtrx_type + "/"
        filename = "signal_err_graph" + "_" + ev_and_pp_type
        os.makedirs(subdir, exist_ok = True)
        fig.savefig(subdir + filename + "_mtrx_" + curr_idx + ".png", dpi = 300)
        plt.close('all')
