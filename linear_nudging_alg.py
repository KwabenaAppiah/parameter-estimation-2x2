from eigenvalues_graph import EigenvaluesGraph
from trace_det_graph import TraceDetGraph
from scipy.integrate import solve_ivp
# import matplotlib.pyplot as plt
from line_graph import LineGraph
from matrix_2x2 import Matrix2x2
import numpy as np
import time
import math
import sys


class LinearNudgingAlg:
    def __init__(self, *args):
        self.select_simulation(*args)

    def get_trace(self, curr_mtrx):
        return np.trace(curr_mtrx)

    def get_det(self, curr_mtrx):
        return np.linalg.det(curr_mtrx)

    def select_simulation(self, *args):
        self._mtrx_list  = []
        self._curr_loop_pos = 0
        self.event_occurred = False
        self.is_custom_solver_on = False
        integration_method = file_import_path = ""
        loop_limit = 0
        self.set_mu_value(args[2])
        self.set_init_mu_value(args[2])
        self.set_relax_time_value(args[3])
        self.set_init_relax_time_value(args[3])

        self._sim_vars_dict = {
            'ev_type': args[0],
            'pp_type': args[1],
            'integration_method': 'BDF',
            'loop_limit': 0,
            'step_val': 10,
            'file_import_path': "",
            'low_bound': 0,
            'high_bound': 0
        }

        self.set_rule(2)
        u_thold = v_thold = ut_thold = vt_thold = np.inf
        threshold_decay_factor = 10
        self.set_thresholds(np.array([u_thold, v_thold, ut_thold, vt_thold, threshold_decay_factor]))

        if len(args) == 5:
            file_import_path = args[4]
            mtrx_list   = self.text_file_to_mtrx_list(file_import_path)
            self.set_mtrx_list(mtrx_list)
            loop_limit = len(self.get_mtrx_list())


        elif len(args) == 6:
            loop_limit =  args[5]
            bound_value = args[4]
            self.set_sim_vars_dict_elt("low_bound", -bound_value)
            self.set_sim_vars_dict_elt("high_bound", bound_value)

        self.set_sim_vars_dict_elt('loop_limit', loop_limit)
        integration_method = self.get_sim_vars_dict_elt('integration_method')
        step_val           = self.get_sim_vars_dict_elt('step_val')
        ev_type            = self.get_sim_vars_dict_elt('ev_type')
        pp_type            = self.get_sim_vars_dict_elt('pp_type')
        self.set_sim_vars_dict_elt('file_import_path', file_import_path)
        bounds = str([self.get_sim_vars_dict_elt("low_bound"), self.get_sim_vars_dict_elt("high_bound")])
        case_type = "MAIN_DIAG_2x2"
        self.set_trace_det_graph(ev_type, pp_type, integration_method, bounds, loop_limit, case_type)
        self.set_line_graph(ev_type, pp_type, integration_method, bounds,loop_limit, case_type)
        self.set_ev_graph(ev_type, pp_type, integration_method, bounds,loop_limit, case_type)
        self.init_simulation(0, case_type)



    def init_simulation(self, loop_start_pos, case_type):
        loop_limit             = self.get_sim_vars_dict_elt('loop_limit')
        i                      = loop_start_pos
        gui_counter            = i + 1

        if loop_limit - loop_start_pos < 0:
            print("ERROR: Loop limit < Loop start position")
            exit()

        elif abs(loop_limit - loop_start_pos) == 0:
            print("The loop cycle has ended.", '\n')
            exit()

        low_bound              = self.get_sim_vars_dict_elt("low_bound")
        high_bound             = self.get_sim_vars_dict_elt("high_bound")
        integration_method     = self.get_sim_vars_dict_elt('integration_method')
        step_val               = self.get_sim_vars_dict_elt('step_val')
        ev_type                = self.get_sim_vars_dict_elt('ev_type')
        pp_type                = self.get_sim_vars_dict_elt('pp_type')
        thresholds             = self.get_thresholds()
        bounds = str([self.get_sim_vars_dict_elt("low_bound"), self.get_sim_vars_dict_elt("high_bound")])
        u_thold, v_thold, ut_thold, vt_thold, threshold_decay_factor = thresholds
        line_graph = self.get_line_graph()
        trace_det_graph = self.get_trace_det_graph()
        ev_graph = self.get_ev_graph()
        mtrx_list  = self.get_mtrx_list()
        mu_1 = mu_4 = self.get_init_mu_value()
        mu_2 = mu_3 = 0

        print("**********************************************", "START - SIMULATION" ,"***************************************************", '\n')
        while i < loop_limit:
            while i != loop_limit:

                print("**********************************", "CYCLE:", gui_counter, "/" , loop_limit, "*************************************************************************", '\n')
                print(ev_type.upper(), "-", pp_type.upper())

                if mtrx_list: # check if mtrx_list is empty
                    a11, a12, a21, a22 = mtrx_list[i][0][0], mtrx_list[i][0][1], mtrx_list[i][1][0], mtrx_list[i][1][1]
                    curr_mtrx = np.array([[a11, a12], [a21, a22]])

                else:
                    curr_mtrx = Matrix2x2(low_bound, high_bound, ev_type, pp_type)
                    a11, a12, a21, a22 = curr_mtrx.get_element(0, 0), curr_mtrx.get_element(0, 1), curr_mtrx.get_element(1, 0), curr_mtrx.get_element(1, 1)
                    curr_mtrx = np.array([[a11, a12], [a21, a22]]) # Created just for the true_trace / true_det

                true_values = np.array([a11, a12, a21, a22])
                true_trace = self.get_trace(curr_mtrx)
                true_det = self.get_det(curr_mtrx)

                #Reset
                self.reset_solutions() # Reset solutions once each matrix is initialized
                self.set_mu_value(self.get_init_mu_value())
                self.set_relax_time_array(self.get_init_relax_time_array())
                self.reset_signal_err()

                # Main Diagonal Case
                a11_0         = a11 + step_val
                a12_0         = a12
                a21_0         = a21
                a22_0         = a22 + step_val

                initial_guesses = np.array([a11_0, a12_0, a21_0, a22_0])
                updates_to_plot = list(true_values != initial_guesses)
                updates_on      = true_values != initial_guesses
                time_between = np.array([self.get_relax_time_value()] * 4) # creates a 1 x 4 array.
                mus = np.array([mu_1, mu_2, mu_3, mu_4])
                params = np.array([a11, a12, a21, a22, a11_0, a12_0, a21_0, a22_0], dtype = np.float64)
                thresholds = np.array([u_thold, v_thold, ut_thold, vt_thold, threshold_decay_factor])

                sim_time = 100 # Stopping time
                dt = 1 / 2 * (1 / self.get_mu_value())
                t_span = [0, sim_time]
                t = np.arange(0, sim_time, dt)

                # ------------ Initialize system --------------------
                S0 = np.array([1, 1, 3, 3])                 # Initialize [x, y, xt, yt], initial state
                guesses = [[a11_0],[a12_0],[a21_0],[a22_0]]  # Guesses

                try:
                    x_dot, y_dot, xt_dot, yt_dot = self.calc_rhs(S0, mus, params)
                    derivs = [[x_dot, y_dot, xt_dot, yt_dot]]
                    pos_err = [[abs(S0[2] - S0[0]), abs(S0[3] - S0[1])]]  # Position error
                    last_updates     = np.array([0, 0, 0, 0])             # Time of last updates
                    idx_last_updates = np.array([0, 0, 0, 0])             # Index of last updates
                    tfe = []                                              # Used to record the times when model() is called

                    # Package into time_tracker arg for solve_ivp()
                    time_tracker = [last_updates, time_between, tfe]
                    sim_args = (mus, params, thresholds, derivs, guesses, time_tracker, updates_on, pos_err, idx_last_updates)

                    sol_time, time_events, guesses, derivs = self.run_simulation(t_span, S0, t, true_values, integration_method, dt, tfe, sim_args )

                    if(case_type == "MAIN_DIAG_2x2"):
                        param_labels = ["a11", "a22"]
                        print("\n", curr_mtrx, '\n')
                        is_tr_det_graph_plotted = trace_det_graph.organize_data(guesses, true_values, param_labels, true_trace, true_det, case_type)
                        is_ev_graph_plotted = ev_graph.organize_data(curr_mtrx, guesses, true_values, param_labels, case_type)

                    elif(case_type == "OFF_DIAG_2x2"):
                        print("Off Diagonal 2x2 - This option is not available yet.")
                        exit()

                    elif(case_type == "BTM_ROW_2x2"):
                        print("Bottom Row 2x2 - This option is not available yet.")
                        exit()

                    if is_tr_det_graph_plotted == True:
                        static_args = ev_type, pp_type, integration_method, bounds, loop_limit, case_type
                        line_graph.organize_data(static_args, sol_time, time_events, guesses, true_values, param_labels, i, 1e-5, self.get_solutions(), self.get_signal_err())

                        i += 1
                        gui_counter += 1

                    else:
                        print('\n', curr_mtrx, '\n' + '\n', 'UNUSABLE curr_mtrx', "\n")
                        print("CYCLE:", gui_counter, "- SKIP THIS ITERATION (1).")
                        if mtrx_list and len(args) == 5: # I.e. If you're importing from a text file
                            i += 1
                            gui_counter += 1
                        else:
                            break

                except ValueError as e:
                    print(e)

                    print('\n', curr_mtrx, '\n' + '\n', 'UNUSABLE MTRX', "\n")
                    print("CYCLE:", gui_counter, "- SKIP THIS ITERATION (2).")
                    break


        print("**********************************************","END - SIMULATION" ,"***************************************************", '\n')
        line_graph.display_avg_param_err_comp("Relative")
        trace_det_graph.display(ev_type, pp_type, integration_method, loop_limit)
        ev_graph.display(ev_type, pp_type, integration_method, loop_limit)
        # plt.show()


    def stop_solve_ivp(self, t, S, *args):
        is_sol_NaN_or_INF = False

        i = 0
        while i < len(S):
            if math.isnan(S[i]):
                print("NaN!")
                is_sol_NaN_or_INF = True

            elif math.isinf(S[i]):
                print("INF!")
                is_sol_NaN_or_INF = True

            i += 1

        return is_sol_NaN_or_INF


    def solve_ivp_bdf(self, f, t_span, S0, h, *sim_args):
        t0, t_end = t_span
        t_values = np.arange(t0, t_end + h, h)  # Generate t values with step size h
        s_values = np.zeros((len(t_values), len(S0)))  # Initialize array to store y values
        s_values[0] = S0  # Set initial value

        for i in range(len(t_values) - 1):
            tn = t_values[i]
            Sn = s_values[i]
            s_values[i + 1] = np.array(Sn) + h * np.array(f(tn, Sn, *sim_args))

        return s_values, t_values



    #2
    def run_simulation(self, t_span, S0, t, true_values, integration_method, dt, tfe, sim_args):
        mus, params, thresholds, derivs, guesses, time_tracker, updates_on, err, idx_last_updates = sim_args

        start = time.time()

        if self.get_is_custom_solver_on() == True:
            # Custom solver is set, built-in solver is not.
            sol = self.solve_ivp_bdf(self.model, t_span, S0, dt, *sim_args)
            y_values, time_events = sol
            sol_time =  time_events # Ideally this var should not be needed for the custom version, but it's being populated here in order to switch easily from the cutom and built-in options

        else:
            # Built-in solver is set, custom solver is not.
            sol = solve_ivp(self.model, t_span = t_span, y0 = S0, method = integration_method, t_eval = t, args = sim_args, events = self.stop_solve_ivp)
            y_values, sol_time, time_events = sol.y_events, sol.t, sol.t_events


        # ---------- Handle output ------------------------
        guesses      = np.array(guesses, dtype = np.float64)
        derivs       = np.array(derivs)
        num_iter     = round(tfe[-1] / dt)
        t_sol        = np.linspace(0, dt * num_iter, num = num_iter)
        f_eval       = np.searchsorted(tfe, t_sol)
        guesses      = guesses[:, f_eval]
        derivs       = derivs[f_eval,:]
        a11s, a12s, a21s, a22s = guesses
        a11,  a12,  a21,  a22  = true_values

        i = 0
        while i < len(y_values[0]):
            x, y, xt, yt = y_values[0][i][0], y_values[0][i][1], y_values[0][i][2], y_values[0][i][3]
            self.set_solutions(x, y, xt, yt)
            self.set_signal_err(x, y, xt, yt)
            i = i + 1

        return sol_time, time_events, guesses, derivs


    #3b
    def model(self, t, S, *sim_args):
        mus, params, thresholds, derivs, guesses, time_tracker, updates_on, err, idx_last_updates = sim_args

        '''
        Function called by odeint to return the time derivative of the system S at time t.

        INPUT

            mus              : NumPy array of shape (4) with values [mu_1, mu_2, mu_3, mu_4]
                               Note a full nudging matrix is not advised (typically choose
                               diagonal or off diagonal elements depending on which params
                               you want to recover)

            params            : NumPy array of shape (8) with values [a11, a12, a21, a22, a11_t, a12_t, a21_t, a22_t]
                               where a11, a12, a21, a22 are the true parameter values and
                               a11_t, a12_t, a21_t, a22_t are the current guesses.

            thresholds       : NumPy array of shape(5) with values [u_thold,
                               ut_thold, v_thold, vt_thold, d] which are used to decide whether to
                               update parameters. Note these are not being used right now, updates
                               are being made any time the relaxation period has lapsed

            derivs           : List of time derivatives at all function evaluations

            guesses          : List of parameter guesses made at all function evals

            time_tracker     : List containing last_updates (the timestamp of the last parameter updates)
                               time_between (relaxation time)
                               tfe (times model() is called in simulation)

            updates_on       : Boolean array specifying which parameters should be updated

            err              : Local position errors (for get_tholds() -- not in use right now)

            idx_last_updates : Index of last updates (for get_tholds() -- not in use right now)
        '''

        a11, a12, a21, a22, a11_t, a12_t, a21_t, a22_t = params
        last_updates, time_between, tfe = time_tracker
        self.set_relax_time_array(time_between)
        curr_time = t
        time_tracker[-1].append(curr_time)      # Record time of function call in tfe

        u_now = abs(S[2] - S[0])       # Current x error |xt - x|
        v_now = abs(S[3] - S[1])       # Current y error |yt - y|

        err.append([u_now, v_now])              # Take the transpose of the np.array 'err'
        pos_err          = np.array(err).T
        time_since       = curr_time - last_updates
        avg_rel_err      = self.get_elt_avg_rel_err(a11, a22, a11_t, a22_t)
        init_relax_value = self.get_init_relax_time_value()


        # --- MAX ADJUSTMENTS  -------------------------------------------------
        if 1e-1 < avg_rel_err and avg_rel_err < math.inf:
            self.set_relax_time_value(init_relax_value * 1/10) # 1/4, ..., 1/10, ... or less resluts in an infinite loop
            self.set_mu_value(self.get_init_mu_value() * 10)   # 4, ..., 10, ... or more resluts in an infinite loop

        elif 1e-3 < avg_rel_err and avg_rel_err <= 1e-1:
            self.set_relax_time_value(init_relax_value * 1/5) # 1/3, 1/5, ... or less resluts in a stall
            self.set_mu_value(self.get_init_mu_value() * 5)   # 3, 5, ... or more results in a stall

        elif 1e-5 < avg_rel_err and avg_rel_err <= 1e-3:
            self.set_relax_time_value(init_relax_value * 1/2) # 1/2 = Max
            self.set_mu_value(self.get_init_mu_value() * 2) # 2 = Max

        elif 1e-8 < avg_rel_err and avg_rel_err <= 1e-5:
            self.set_relax_time_value(init_relax_value * 2)
            self.set_mu_value(self.get_init_mu_value() * 1/2)

        elif 0 <= avg_rel_err and avg_rel_err <= 1e-8:
            self.set_relax_time_value(init_relax_value * 4)
            self.set_mu_value(self.get_init_mu_value() * 1/4)

        mus[0] = mus[3] = self.get_mu_value()
        # # ----------------------------------------------------------------------

        time_threshold = time_since > self.get_relax_time_array() # Check if relaxation time has lapsed for parameters
        stop_update = 1e-10  # Threshold on relative position error for stopping updates

        # Just for sinks - add condition later
        if(self.get_sim_vars_dict_elt('pp_type') == "sink" or self.get_sim_vars_dict_elt('pp_type') == "sp_sink"):
            if u_now  < stop_update and v_now  < stop_update:
                updates_on[0] = updates_on[1] = updates_on[2] = updates_on[3] = False
        else:
            if u_now / abs(S[0]) < stop_update and v_now / abs(S[1]) < stop_update:
                updates_on[0] = updates_on[1] = updates_on[2] = updates_on[3] = False

        # Update params where relaxation period has lapsed and stop criteria not met
        to_update  = np.logical_and(time_threshold, updates_on)
        update_idx = [i for i, x in enumerate(to_update) if x]

        self.update_params(update_idx, S, mus, params, guesses, thresholds, self.get_rule(), time_tracker, pos_err, idx_last_updates)
        no_update_idx = [i for i, x in enumerate(to_update) if not x]

        for i in no_update_idx:
            guesses[i].append(params[i + 4])

        x_dot, y_dot, xt_dot, yt_dot = self.calc_rhs(S, mus, params)
        St = [x_dot, y_dot, xt_dot, yt_dot]
        derivs.append(St)


        return St


    def update_params(self, which, S, mus, params, guesses, thresholds, rule, time_tracker, pos_err, idx_last_updates):
        '''
        Updates parameter guesses [if thresholds are met]
        '''
        curr_time = time_tracker[-1][-1]

        # Extract input
        x, y, xt, yt = S
        mu_1, mu_2, mu_3, mu_4 = mus
        a11, a12, a21, a22, a11_t, a12_t, a21_t, a22_t = params
        u_thold, v_thold, ut_thold, vt_thold, d = thresholds

        # Calculate derivatives
        dx_dt, dy_dt, dxt_dt, dyt_dt = self.calc_rhs(S, mus, params)
        d_vals = [dx_dt, dy_dt, dxt_dt, dyt_dt]
        d_names = ['dx_dt', 'dy_dt', 'dxt_dt', 'dyt_dt']

        if self.is_nan_inf_or_overflow(d_vals, d_names):
            print('Exit function: update_params().')
            print(S)
            return

        else:
            derivs_now  = np.array([dx_dt, dy_dt, dxt_dt, dyt_dt])
            u      = xt - x
            v      = yt - y
            u_dot = dxt_dt - dx_dt
            v_dot = dyt_dt - dy_dt
            errors = np.array([u, v, u_dot, v_dot])

            if abs(u) / abs(x) <= u_thold and abs(v) / abs(y) <= v_thold:

                for i in which:
                    new_guess = self.update_formula(rule, S, errors, derivs_now, mus, params)[i]
                    params[i + 4] = new_guess
                    guesses[i].append(new_guess)
                    time_tracker[0][i] = curr_time
                    idx_last_updates[i] = pos_err.shape[1] - 1

            else:
                for i in update_idx:
                    parm_idx = i + 4
                    guesses[i].append(parm_idx)


    def get_elt_avg_rel_err(self, a_true, b_true, a_guess, b_guess):
        a_rel_err = abs(a_true - a_guess ) / abs(a_true)
        b_rel_err = abs(b_true - b_guess) / abs(b_true)
        return abs(a_rel_err + b_rel_err) / 2



    def calc_rhs(self, S, mus, params):
        '''
        Returns derivatives of x, y, x, yt
        '''
        x, y, xt, yt = S
        mu_1, mu_2, mu_3, mu_4 = mus
        a11, a12, a21, a22, a11_t, a12_t, a21_t, a22_t = params

        try:
            x_dot = a11 * x + a12 * y
            y_dot = a21 * x + a22 * y
            xt_dot = a11_t * xt + a12_t * yt - mu_1 * (xt - x) - mu_2 * (yt - y)
            yt_dot = a21_t  * xt + a22_t * yt - mu_3 * (xt - x) - mu_4 * (yt - y)

        except BaseException:
            print("\n" + "WARNING: " + "Skipping current matrix: "  + str(curr_mtrx) + "." , '\n')
            pass

        return x_dot, y_dot, xt_dot, yt_dot



    def is_nan_or_inf(self, value, var_name):
        bool = False
        if math.isinf(value):
            bool = True

        elif math.isnan(value):
            bool = True

        return bool


    def is_over_or_under_flow(self, value, var_name):
        bool = False
        if abs(value) > sys.float_info.max:
            print("\n"+ "ERROR: " + var_name + ": Overflow may occur")
            bool = True

        elif abs(value) < sys.float_info.min:
            print("\n"+  "ERROR: " + var_name + ": Underflow may occur")
            bool = True

        return bool


    def is_nan_inf_or_overflow(self, vals, names):
        i = 0
        bool = False
        while i < len(vals):
            if self.is_nan_or_inf(vals[i], names[i]):
                bool = True

            elif self.is_over_or_under_flow(vals[i], names[i]):
                bool = True
            i = i + 1
        return bool


    def update_formula(self, rule, S, errors, derivs_now, mus, params):
        '''
        Returns parameter update formulas.
        '''
        x, y, xt, yt = S
        u, v, u_dot, v_dot = errors
        mu_1, mu_2, mu_3, mu_4 = mus
        a11, a12, a21, a22, a11_t, a12_t, a21_t, a22_t = params

        switcher = {
            'constant'  : [a11_t, a12_t, a21_t, a22_t], # Why is this necessary?

            'exact':      [a11_t - mu_1 * (u / xt) - mu_2 * (v / xt) - (u_dot / xt),
                           a12_t - mu_1 * (u / yt) - mu_2 * (v / yt) - (u_dot / yt),
                           a21_t - mu_3 * (u / xt) - mu_4 * (v / xt) - (v_dot / xt),
                           a22_t - mu_3 * (u / yt) - mu_4 * (v / yt) - (v_dot / yt)],

            'drop_deriv': [a11_t - mu_1 * (u / xt) - mu_2 * (v / xt),
                           a12_t - mu_1 * (u / yt) - mu_2 * (v / yt),
                           a21_t - mu_3 * (u / xt) - mu_4 * (v / xt),
                           a22_t - mu_3 * (u / yt) - mu_4 * (v / yt)]

        }
        return switcher.get(rule, "Update rule not recognized")


    def text_file_to_mtrx_list(self, imported_file):
        mtrx_list = []
        with open(imported_file, 'r') as file:
            for line in file:
                line = line.strip()  # removes any leading or trailing white space
                if line:  # only processes lines with text (not just white space)
                    str = line.strip()
                    vars = str.split("[")[1].split("]")[0]
                    mtrx_elts = vars.split()
                    a11, a12, a21, a22 = mtrx_elts
                    temp_mtrx = np.array([[float(a11), float(a12)], [float(a21), float(a22)]])
                    mtrx_list.append(temp_mtrx)
        return mtrx_list


    def rule_string(self, which):
        '''
        Converts rule index to string for retrieving update formula from dictionary
        '''
        switcher = {0: 'constant', 1: 'exact', 2: 'drop_deriv'}
        return switcher.get(which)


    #-------- SETTERS -------------------------------------------------------------------

    def get_curr_loop_pos(self):
        return self._curr_loop_pos

    def get_ev_graph(self):
        return self._ev_graph

    def get_init_mu_value(self):
        return self.init_mu_value

    def get_init_relax_time_array(self):
        return self._init_relax_time_array

    def get_init_relax_time_value(self):
        return self._init_relax_time_array[0]

    def get_is_custom_solver_on(self):
        return self.is_custom_solver_on

    def get_line_graph(self):
        return self._line_graph

    def get_mu_value(self):
        return self._mu_value

    def get_relax_time_array(self):
        return self._relax_time_array

    def get_relax_time_value(self):
        return self._relax_time_array[0]

    def get_rule(self):
        return self._rule

    def get_signal_err(self):
        return self._x_signal_err_list, self._y_signal_err_list

    def get_sim_vars_dict(self):
        return self._sim_vars_dict

    def get_solutions(self):

        return self._x_list, self._y_list, self._xt_list, self._yt_list


    def get_sim_vars_dict_elt(self, key):
        return self._sim_vars_dict[key]

    def get_thresholds(self):
        return self._thresholds

    def get_trace_det_graph(self):
        return self._trace_det_graph


    def get_mtrx_list(self):
        return self._mtrx_list

    def get_mtrx_list_elt(self, index):
        return self._mtrx_list[index]

    #-------- SETTERS -------------------------------------------------------------------
    def set_mtrx_list(self, list):
        self._mtrx_list = list

    def set_mtrx_list_elt(self, index, value):
        self._mtrx_list[index] = value


    def set_mu_value(self, value):
        self._mu_value = value

    def set_init_mu_value(self, value):
        self.init_mu_value = value

    def set_relax_time_value(self, value):
        self._relax_time_array  = np.array([value] * 4)

    def set_relax_time_array(self, np_array):
        self._relax_time_array = np_array

    def set_init_relax_time_value(self, value):
        self._init_relax_time_array  = np.array([value] * 4) # I.e. A (1 x 4) numpy array

    def set_is_custom_solver_on(self, bool):
        self.is_custom_solver_on = bool

    def set_rule(self, rule_to_use):
        self._rule = self.rule_string(rule_to_use)

    def set_solutions(self, x, y, xt, yt):
        self._x_list.append(x)
        self._y_list.append(y)
        self._xt_list.append(xt)
        self._yt_list.append(yt)

    def set_signal_err(self, x, y, xt, yt):

        if(self.get_sim_vars_dict_elt('pp_type') == "sink" or self.get_sim_vars_dict_elt('pp_type') == "sp_sink"):
            self._x_signal_err_list.append(abs(xt - x))
            self._y_signal_err_list.append(abs(yt - y))
        
        else:
            self._x_signal_err_list.append(abs(xt - x) / x)
            self._y_signal_err_list.append(abs(yt - y) / y)


    def set_curr_loop_pos(self, value):
        self._curr_loop_pos = value

    def set_thresholds(self, arr):
        self._thresholds = arr

    def set_sim_vars_dict_elt(self, key, value):
        self._sim_vars_dict[key] = value

    def set_line_graph(self, ev_type, pp_type, integration_method, bounds, loop_limit, case_type):
        self._line_graph = LineGraph(ev_type, pp_type, integration_method, bounds, loop_limit, case_type)

    def set_trace_det_graph(self, ev_type, pp_type, integration_method, bounds, loop_limit, case_type):
        self._trace_det_graph = TraceDetGraph(ev_type, pp_type, integration_method, bounds, loop_limit, case_type)

    def set_ev_graph(self, ev_type, pp_type, integration_method, bounds, loop_limit, case_type):
        self._ev_graph = EigenvaluesGraph(ev_type, pp_type, integration_method, bounds, loop_limit, case_type)

    #--------------------------------------------------------------

    def reset_signal_err(self):
        self._x_signal_err_list = []
        self._y_signal_err_list = []

    def reset_solutions(self):
        self._x_list = []
        self._y_list = []
        self._xt_list = []
        self._yt_list = []
