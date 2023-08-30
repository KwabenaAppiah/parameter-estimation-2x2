from line_graph import LineGraph
from non_line_graph import NonLineGraph
from matrix_2x2 import Matrix_2x2
import numpy as np


class LinearNudgingAlg:
    def __init__(self, *args):
        self.prepare_sim(args)

    def prepare_sim(self, args):
        ev_type, pp_type, mu_value, relax_time, case_type = args[0], args[1], args[2], args[3], args[4]
        loop_limit, bound_value = 0, 0
        matrix_estimates = []
        file_import_path = ""

        if len(args) == 6:
            file_import_path = args[5]
            matrix_estimates = self.text_file_to_mtrx_estimates(file_import_path)
            loop_limit = len(matrix_estimates)


        elif len(args) == 7:
            bound_value = args[5]
            loop_limit =  args[6]

        self.init_sim(ev_type, pp_type, mu_value, relax_time, bound_value, matrix_estimates, loop_limit, case_type)

    def init_sim(self, ev_type, pp_type, mu_value, relax_time, bound_value, matrix_estimates, loop_limit, case_type):

        # Initialize display graphs
        threshold_vals = [1e-12, 1e-8, 1e-4, 1e-1]
        line_graph = LineGraph(ev_type, pp_type, [-bound_value, bound_value], loop_limit, case_type)
        non_line_graph = NonLineGraph(ev_type, pp_type, [-bound_value, bound_value], loop_limit, case_type, threshold_vals)
        # threshold = 1e-8 # Values below this threshold are considered bad, while values above are considered good
        threshold = threshold_vals[1] #I.e. 1e-8, Values below this threshold are considered bad, while values above are considered good

        i = 0
        while i < loop_limit:
            print("\n" + "-------- MATRIX:", i, "---------------------------------------------------------------------------------------------" + "\n")
            # # Initialize Matrix
            # A = self.generate_matrix(2, bound_value)
            # print("Matrix type:", self.classify_matrix(A), "\n")

            if matrix_estimates: # check if matrix_estimates is empty
                a11, a12, a21, a22 = matrix_estimates[i][0][0], matrix_estimates[i][0][1], matrix_estimates[i][1][0], matrix_estimates[i][1][1]
                A = np.array([[a11, a12], [a21, a22]])

            else:
                curr_mtrx = Matrix_2x2(-bound_value, bound_value, ev_type, pp_type)
                a11, a12, a21, a22 = curr_mtrx.get_element(0, 0), curr_mtrx.get_element(0, 1), curr_mtrx.get_element(1, 0), curr_mtrx.get_element(1, 1)
                A = np.array([[a11, a12], [a21, a22]])

            trace_A = np.trace(A)
            det_A = np.linalg.det(A)
            true_params = (a11, a12, a21, a22)

            # Simulation parameters
            sim_time = 10    # Stop time
            dt = 0.001       # Time step
            t = np.arange(0, sim_time, dt)   # Evaluation times
            T_R = relax_time # Relaxation period, default = 0.5
            init_param_err = 10

            # Stored update times
            a11_update_times, a12_update_times, a21_update_times, a22_update_times = np.array([0]), np.array([0]), np.array([0]), np.array([0])
            update_times = (a11_update_times, a12_update_times, a21_update_times, a22_update_times)

            # System parameters
            a11, a12, a21, a22 = A[0, 0], A[0, 1], A[1, 0], A[1, 1]
            if(case_type == "main_diagonal"):
                a11_t, a12_t, a21_t, a22_t = a11 + init_param_err, a12, a21, a22 + init_param_err  # a11, a22


            elif(case_type == "anti-diagonal"):
                a11_t, a12_t, a21_t, a22_t = a11, a12 + init_param_err, a21 + init_param_err, a22 # a12, a21


            elif(case_type == "left_column"):
                a11_t, a12_t, a21_t, a22_t = a11 + init_param_err, a12, a21 + init_param_err, a22 # a11, a21


            elif(case_type == "right_column"):
                a11_t, a12_t, a21_t, a22_t = a11, a12 + init_param_err, a21, a22 + init_param_err  # a12, a22


            At = np.array([[a11_t, a12_t], [a21_t, a22_t]])

            # Stored parameter estimates
            a11_estimates, a12_estimates, a21_estimates, a22_estimates  = np.array([a11_t]), np.array([a12_t]), np.array([a21_t]), np.array([a22_t])
            param_estimates = (a11_estimates, a12_estimates, a21_estimates, a22_estimates)

            # Nudging parameters
            mu_11, mu_12, mu_21, mu_22 = mu_value, 0, 0, mu_value
            M = np.array([[mu_11, mu_12], [mu_21, mu_22]])

            # Initialize system
            x_0, y_0, xt_0, yt_0 = 1, 1, 3, 3
            S = np.array([[x_0, y_0, xt_0, yt_0]])  # Each new solution point makes a new row
            U = np.array([[S[0, 2] - S[0, 0], S[0, 3] - S[0, 1]]])  # Signal error terms: U = [u, v] = [xt - x, yt - y]

            # Run the algorithm
            S_lists, U_lists, U_lists_abs_err,  param_estimates, update_times, param_errors, avg_param_errors = self.nuding_algorithm(A, At, t, dt, T_R, param_estimates, update_times, S, M, U, case_type)
            line_graph.organize_data(t, update_times, A, param_estimates, param_errors, avg_param_errors, i, threshold, S_lists, U_lists_abs_err)
            non_line_graph.organize_data(A, param_estimates, avg_param_errors, trace_A, det_A, case_type)
            i = i + 1

        print("\n" + "-------- END SIMULATION ---------------------------------------------------------------------------------------------" + "\n")

        # trace_det_graph.display(ev_type, pp_type, case_type, loop_limit)
        # ev_graph.display(ev_type, pp_type, case_type, loop_limit)
        non_line_graph.display(ev_type, pp_type, case_type, loop_limit)

    def F(self, t, A, At, S, M, U):
        X_dot = np.matmul(A, S[t, 0:2])
        Xt_dot = np.matmul(At, S[t, 2:]) - np.matmul(M, U[t])
        S_dot = np.append(X_dot, Xt_dot)
        return S_dot

    def nuding_algorithm(self, A, At, t, dt, T_R, param_estimates, update_times, S, M, U, case_type):
        a11_update_times, a12_update_times, a21_update_times, a22_update_times  = update_times
        a11_estimates, a12_estimates, a21_estimates, a22_estimates = param_estimates

        a11, a12, a21, a22 = A[0, 0], A[0, 1], A[1, 0], A[1, 1]
        mu_11, mu_12, m_21, mu_22 = M[0, 0], M[0, 1], M[1, 0], M[1, 1]

        a11_err, a12_err = [a11_estimates[0] - a11], [a12_estimates[0] - a12]
        a21_err, a22_err = [a21_estimates[0] - a21], [a22_estimates[0] - a22]

        a11_abs_err, a12_abs_err = [abs(a11_err[0])], [abs(a12_err[0])]
        a21_abs_err, a22_abs_err = [abs(a21_err[0])], [abs(a22_err[0])]

        a11_rel_err, a12_rel_err = [abs(a11_err[0] / a11)],  [abs(a12_err[0] / a12)]
        a21_rel_err, a22_rel_err  =[abs(a21_err[0] / a21)], [abs(a22_err[0] / a22)]


        x_estimates, y_estimates, xt_estimates, yt_estimates = [S[0, 0]], [S[0, 1]], [S[0, 2]], [S[0, 3]]
        U1_indx_j = S1_indx_j = At1_indx_i = At1_indx_j = U2_indx_j = S2_indx_j = At2_indx_i = At2_indx_j = -1
        x_sig_err_estimates, y_sig_err_estimates = [U[0, 0]], [U[0, 1]]
        x_sig_abs_err_estimates, y_sig_abs_err_estimates = [abs(U[0, 0])], [abs(U[0, 1])]
        param_1_update_times, param_2_update_times = [], []
        param_1_estimates, param_2_estimates = [], []
        param_1_err, param_2_err = [], []
        param_1_abs_err, param_2_abs_err = [], []
        param_1_rel_err, param_2_rel_err = [], []
        param_1, param_2 = 0, 0

        if(case_type == "main_diagonal"):
            U1_indx_j, S1_indx_j, At1_indx_i, At1_indx_j = 0, 2, 0, 0
            U2_indx_j, S2_indx_j, At2_indx_i, At2_indx_j = 1, 3, 1, 1
            param_1_estimates, param_2_estimates = a11_estimates, a22_estimates
            param_1_update_times, param_2_update_times = a11_update_times, a22_update_times
            param_1_abs_err, param_2_abs_err = a11_abs_err, a22_abs_err
            param_1_rel_err, param_2_rel_err = a11_rel_err, a22_rel_err
            param_1, param_2 = a11, a22

        elif(case_type == "anti-diagonal"):
            U1_indx_j, S1_indx_j, At1_indx_i, At1_indx_j = 0, 3, 0, 1
            U2_indx_j, S2_indx_j, At2_indx_i, At2_indx_j = 1, 2, 1, 0
            param_1_estimates, param_2_estimates = a12_estimates, a21_estimates
            param_1_update_times, param_2_update_times = a12_update_times, a21_update_times
            param_1_abs_err, param_2_abs_err = a12_abs_err, a21_abs_err
            param_1_rel_err, param_2_rel_err = a12_rel_err, a21_rel_err
            param_1, param_2 = a12, a21

        elif(case_type == "left_column"):
            U1_indx_j, S1_indx_j, At1_indx_i, At1_indx_j = 0, 2, 0, 0
            U2_indx_j, S2_indx_j, At2_indx_i, At2_indx_j = 1, 2, 1, 0
            param_1_estimates, param_2_estimates = a11_estimates, a21_estimates
            param_1_update_times, param_2_update_times = a11_update_times, a21_update_times
            param_1_abs_err, param_2_abs_err = a11_abs_err, a21_abs_err
            param_1_rel_err, param_2_rel_err = a11_rel_err, a21_rel_err
            param_1, param_2 = a11, a21

        elif(case_type == "right_column"):
            U1_indx_j, S1_indx_j, At1_indx_i, At1_indx_j = 0, 3, 0, 1
            U2_indx_j, S2_indx_j, At2_indx_i, At2_indx_j = 1, 3, 1, 1
            param_1_estimates, param_2_estimates = a12_estimates, a22_estimates
            param_1_update_times, param_2_update_times = a12_update_times, a22_update_times
            param_1_abs_err, param_2_abs_err = a12_abs_err, a22_abs_err
            param_1_rel_err, param_2_rel_err = a12_rel_err, a22_rel_err
            param_1, param_2 = a12, a22

        # Run nudging algorithm (main diagonal)
        for i in range(len(t) - 1):
            # If |u(i)| ≥ 0 and xt(i) != 0 and relaxation period passed:

            if abs(U[i, U1_indx_j]) >= 0 and S[i, S1_indx_j] != 0 and t[i] - param_1_update_times[-1] >= T_R:
                # Update the first parameter
                At[At1_indx_i, At1_indx_j] = At[At1_indx_i, At1_indx_j] - mu_11 * U[i, U1_indx_j] / S[i, S1_indx_j]

                # Record new estimate
                new_param_1_estimate = At[At1_indx_i, At1_indx_j]
                param_1_estimates = np.append(param_1_estimates, new_param_1_estimate)

                # Record update time
                param_1_update_times = np.append(param_1_update_times, t[i])

            # If |v(i)| ≥ 0 and yt(i) != 0 and relaxation period passed:
            if abs(U[i, U2_indx_j]) >= 0 and S[i, S2_indx_j] != 0 and t[i] - param_2_update_times[-1] >= T_R:

                # Update the second parameter
                At[At2_indx_i, At2_indx_j] = At[At2_indx_i, At2_indx_j] - mu_22 * U[i, U2_indx_j] / S[i, S2_indx_j]

                # Record new param estimate
                new_param_2_estimate = At[At2_indx_i, At2_indx_j]
                param_2_estimates = np.append(param_2_estimates, new_param_2_estimate)

                # Record update time
                param_2_update_times = np.append(param_2_update_times, t[i])

            # Integrate coupled reference and auxiliary system (forward Euler method)
            new_row = S[i] + self.F(i, A, At, S, M, U) * dt
            S = np.vstack((S, new_row))

            # Create lists for x, y, xt and yt for graphing purposes
            x_estimates.append(S[i, 0])
            y_estimates.append(S[i, 1])
            xt_estimates.append(S[i, 2])
            yt_estimates.append(S[i, 3])

            # Record new state error
            new_err = np.array([S[i + 1, 2] - S[i + 1, 0], S[i + 1, 3] - S[i + 1, 1]])
            U = np.vstack((U, new_err))
            x_sig_err_estimates.append(U[i, 0])
            y_sig_err_estimates.append(U[i, 1])
            x_sig_abs_err_estimates.append(abs(x_sig_err_estimates[-1]))
            y_sig_abs_err_estimates.append(abs(y_sig_err_estimates[-1]))

            #Calculate paramemter error
            param_1_err.append(param_1_estimates[-1] - param_1)
            param_2_err.append(param_2_estimates[-1] - param_2)

            #Calculate absolute paramemter error
            param_1_abs_err.append(abs(param_1_err[-1]))
            param_2_abs_err.append(abs(param_2_err[-1]))

            #Calculate relative paramemter error
            param_1_rel_err.append(abs(param_1_err[-1] / param_1))
            param_2_rel_err.append(abs(param_2_err[-1] / param_2))


        param_errors = param_1_abs_err, param_2_abs_err, param_1_rel_err, param_2_rel_err
        param_estimates = (param_1_estimates, param_2_estimates)
        params = (param_1, param_2)
        update_times = (param_1_update_times, param_2_update_times)

        self.print_alg_output(A, t, S, U, params, param_estimates, param_errors, update_times, case_type)
        S_lists = (x_estimates, y_estimates, xt_estimates, yt_estimates)
        U_lists = (x_sig_err_estimates, y_sig_err_estimates)
        U_lists_abs_err = (x_sig_abs_err_estimates, y_sig_abs_err_estimates)

        avg_abs_param_err = self.get_avg_list(param_1_abs_err, param_2_abs_err)
        avg_rel_param_err = self.get_avg_list(param_1_rel_err, param_2_rel_err)
        avg_param_errors = (avg_abs_param_err, avg_rel_param_err)

        return S_lists, U_lists, U_lists_abs_err, param_estimates, update_times, param_errors, avg_param_errors


    def print_alg_output(self, A, t, S, U, params, param_estimates, param_errors, update_times, case_type):
        param_1_abs_err, param_2_abs_err, param_1_rel_err, param_2_rel_err = param_errors
        param_1_estimates, param_2_estimates = param_estimates
        param_1, param_2 = params

        if(case_type == "main_diagonal"):
            param_1_label, param_2_label = "a11", "a22"

        elif(case_type == "anti-diagonal"):
            param_1_label, param_2_label = "a12", "a21"

        elif(case_type == "left_column"):
            param_1_label, param_2_label = "a11", "a21"

        elif(case_type == "right_column"):
            param_1_label, param_2_label = "a12", "a22"


        print("Matrix:")
        print(A, "\n" + "\n")

        print(f"X = \n {S[:, 0:2]}\n")
        print(f"Xt = \n {S[:, 2:]}\n", "\n" + "\n" )

        # Param 1
        print("........", param_1_label, "........................................", "\n" )
        print("True value =", param_1, "\n")
        print("Final estimate =", param_1_estimates[-1], "\n")
        print("Absolute error =", param_1_abs_err[-1], "\n")
        print("Relative error =", param_1_rel_err[-1], "\n")
        print("Total estimates = \n", param_1_estimates, "\n" + "\n")

        # Param 2
        print("........", param_2_label, ".......................................", "\n"  )
        print("True value =", param_2, "\n")
        print("Final estimate =", param_2_estimates[-1], "\n")
        print("Absolute error =", param_2_abs_err[-1], "\n")
        print("Relative error =", param_2_rel_err[-1], "\n")
        print("Total estimates = \n", param_2_estimates, "\n" + "\n")


    def get_avg_list(self, list_1, list_2):
        avg_list = []
        for i in range(len(list_1)):
            avg_list.append(list_1[i] + list_2[i] / 2)
        return avg_list


    def text_file_to_mtrx_estimates(self, imported_file):
        mtrx_estimates = []
        with open(imported_file, 'r') as file:
            for line in file:
                line = line.strip()  # removes any leading or trailing white space
                if line:  # only processes lines with text (not just white space)
                    str = line.strip()
                    vars = str.split("(")[1].split(")")[0]
                    mtrx_elts = vars.split(',')
                    a11, a12, a21, a22 = mtrx_elts
                    temp_mtrx = np.array([[float(a11), float(a12)], [float(a21), float(a22)]])
                    mtrx_estimates.append(temp_mtrx)
        return mtrx_estimates
