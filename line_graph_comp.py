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
    def get_bad_matrices(self):
        return self._mtrx_dict

    def get_bad_matrices_exist(self):
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
    def set_bad_matrices_exist(self, bool):
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
        alphas, betas, deltas, gammas = guesses
        alpha, beta, delta, gamma  = true_vals
        total_err_steps = []
        rel_delta_err = []
        rel_gamma_err = []

        if delta != 0 and gamma != 0:
            j = 0
            for i in range(len(deltas)):
                rel_delta_err.append(abs(deltas[i] - delta) / abs(delta))
                rel_gamma_err.append(abs(gammas[i] - gamma) / abs(gamma))
                total_err_steps.append( abs(rel_delta_err[i] + rel_gamma_err[i]) / 2)

            #Select based on final val (i.e. like Tr / Det plane)
            rel_delta_err_final = abs(deltas[-1] - delta) / abs(delta)
            rel_gamma_err_final = abs(gammas[-1] - gamma) / abs(gamma)
            rel_avg_dg_err_final = abs(rel_delta_err_final +  rel_gamma_err_final) / 2

            if rel_avg_dg_err_final  > 1e-3:
                self.plot_graph(cycle_num, total_err_steps)
                self.set_bad_matrices_exist(True)
                return True

        else:
            # print("Delta:", abs(delta), "and Gamma:", abs(gamma) )
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
        if self.get_bad_matrices_exist() != False:
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
