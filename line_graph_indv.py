import matplotlib.pyplot as plt
import numpy as np
import datetime
import time
import math
import os


class LineGraphIndv:
    def __init__(self, ev_type, pp_type):
        self._ev_type = ev_type
        self._pp_type = pp_type

    def get_ev_type(self):
        return self._ev_type

    def get_pp_type(self):
        return self._pp_type

    def plot_graph(self, guesses, true_vals, cycle_num):
        alphas, betas, deltas, gammas = guesses
        alpha, beta, delta, gamma  = true_vals
        total_err_steps = []
        rel_delta_err = []
        rel_gamma_err = []

        if delta != 0 and gamma != 0:
            for i in range(len(deltas)):
                rel_delta_err.append(abs(deltas[i] - delta) / abs(delta))
                rel_gamma_err.append(abs(gammas[i] - gamma) / abs(gamma))
                total_err_steps.append( abs(rel_delta_err[i] + rel_gamma_err[i]) / 2)

            #Select color based on final val (i.e. like Tr / Det plane)
            rel_delta_err_final = abs(deltas[-1] - delta) / abs(delta)
            rel_gamma_err_final = abs(gammas[-1] - gamma) / abs(gamma)
            rel_avg_dg_err_final = abs(rel_delta_err_final +  rel_gamma_err_final) / 2

            if rel_avg_dg_err_final  > 1e-3:
                self.display(true_vals, cycle_num, total_err_steps)
                self.write_to_file("Mtrx. " + str(cycle_num) + ": " + str(true_vals), "Avg. Err. :" + str(rel_avg_dg_err_final))
                #self.set_bad_matrices("MTRX: " + str(cycle_num) + str(true_vals) +"\n" + "AVG. ERR: " + str(rel_avg_dg_err_final))
                # self.plot_fig_comp(cycle_num, total_err_steps)
                return True

        else:
            # print("Delta:", abs(delta), "and Gamma:", abs(gamma) )
            # print("Mtrx.", true_vals, "is not plottable.")
            return False


    def write_to_file(self, str_1, str_2):
        output = str_1 + "\n" + str_2 + "\n" + "\n"
        ev_pp_type = self.get_ev_type() + '_' + self.get_pp_type()
        filename = 'bad_matrices'  + '_' + ev_pp_type + '.txt'
        subdir = "../output/" + ev_pp_type + "_" + self.get_date_str() + "/"
        os.makedirs(subdir, exist_ok = True)
        f = open(subdir + filename, 'a')
        f.write(output)
        f.close()


    def get_date_str(self):
        raw_date = datetime.datetime.now()
        t = time.localtime()
        # current_time = time.strftime("%H.%M", t)
        # date_str = str(raw_date.year) + "." + str(time.strftime('%m')) + "." + str(time.strftime('%d')) + current_time
        date_str = str(raw_date.year) + "." + str(time.strftime('%m')) + "." + str(time.strftime('%d'))
        return date_str


    def display(self, true_vals, cycle_num, total_err_steps):
        ev_type = self.get_ev_type()
        pp_type = self.get_pp_type()
        fig, ax = plt.subplots()
        fig.set_size_inches(7, 5)
        ax.plot(list(range(len(total_err_steps))), total_err_steps, label = "MX" + " "+ str(cycle_num))
        title = ev_type.upper() + " | " + pp_type.upper() + " - MATRIX " + str(cycle_num) + " : " + str(true_vals)
        ax.set_title(label = title, pad = 20)
        ax.set_xlabel("Iterations")
        ax.set_ylabel("Avg. Rel. Errors")
        ax.set_yscale('log')

        #Output IMG files
        output_type = "indv"
        ev_pp_type = self.get_ev_type() + '_' + self.get_pp_type()
        subdir =   "../output/" + ev_pp_type + "_" + self.get_date_str() + "/line_graph" + "_" + output_type + "/"
        filename =  "line_graph" + "_" + ev_pp_type
        os.makedirs(subdir, exist_ok = True)
        fig.savefig(subdir + filename + "_mtrx_" + str(cycle_num) + ".png", dpi = 300)
        plt.close(fig)
