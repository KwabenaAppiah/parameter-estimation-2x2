import math
import random
import sys
import time
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

# nudge, project-specific
from tr_det_graph import TrDetGraph
from line_graph import LineGraph
from matrix_2x2 import Matrix2x2
# from sol_graph import SolGraph

class LinearNudgingAlg:
    def __init__(self, *args):
        ev_type, pp_type, mu_val, relax_time = args[0], args[1], args[2], args[3]
        mtrx_list = []
        self.set_mu_value(mu_val)
        self.set_init_mu_value(mu_val)

        self.set_relax_time_value(relax_time)
        self.set_init_relax_time_value(relax_time)

        # has_5_args = has_6_args = False
        has_5_args = False
        if len(args) == 5:
            file_import = args[4]
            mtrx_list = self.text_file_to_mtrx_list(file_import)
            loop_limit = len(mtrx_list)
            has_5_args = True

        else:
        # if len(args) == 6:
            bound_val = args[4]
            loop_limit = args[5]
            low_bnd, high_bnd = -bound_val, bound_val

        self.set_rule(2)
        # ------------ Algorithm parameters ----------------
        # Nudging parameters
        # NOTE: Full nudging matrix not advised. Diagonal entries seem to work best but
        #       could play around with off diagonal entries. I typically use multiple of
        #       identity matrix
        mu_1 = mu_4 = self.get_mu_value()
        mu_2 = mu_3 = 0

        # Position, velocity thresholds for updates
        # NOTE: Right now we are not using any thresholds, they are all set to infinity with no decay
        u_thold    = np.inf
        v_thold    = np.inf
        ut_thold   = np.inf
        vt_thold   = np.inf
        d          = 10         # Threshold decay factor
        step_val   = 10         # The step value for a22 and a21
        tr_det_graph = TrDetGraph(ev_type, pp_type, loop_limit, mu_val)
        line_graph = LineGraph(ev_type, pp_type)

        print("**********************************************", "START - SIMULATION" ,"***************************************************", '\n')
        i = 0
        gui_counter = i + 1

        while i < loop_limit:
            while i != loop_limit:

                print("**********************************", "CYCLE:", gui_counter, "/" , loop_limit, "*************************************************************************", '\n')
                print(ev_type.upper(), "-", pp_type.upper())

                if has_5_args == True:
                    a11, a12, a21, a22 = mtrx_list[i][0][0], mtrx_list[i][0][1], mtrx_list[i][1][0], mtrx_list[i][1][1]
                    mtrx = np.array([[a11, a12], [a21, a22]])
                    # if np.any(np.isnan(mtrx)) == True and np.any(np.isinf(mtrx)) == True:
                    #     print("Pre-skip!")
                    #     break
                else:
                    mtrx = Matrix2x2(low_bnd, high_bnd, ev_type, pp_type)
                    # if np.any(np.isnan(mtrx)) == False and np.any(np.isinf(mtrx)) == False:
                    a11, a12, a21, a22 = mtrx.get_element(0, 0), mtrx.get_element(0, 1), mtrx.get_element(1, 0), mtrx.get_element(1, 1)

                true_vals = np.array([a11, a12, a21, a22])

                #Reset
                self.reset_solutions() # Reset solutions once each matrix is initialized
                self.set_is_mtrx_usable(True)
                self.set_mu_value(self.get_init_mu_value())
                self.set_relax_time_array(self.get_init_relax_time_array())
                self.reset_signal_err()
                # print("init mu value", self.get_init_mu_value())

                # NOTE: Code is set up to only update parameters with guesses different from
                #       the true value. This needs to be toggled within the script now but could
                #       be implemented from the command line if you wanted
                a11_0         = a11
                a12_0         = a12
                a21_0         = a21 + step_val
                a22_0         = a22 + step_val
                initial_guesses = np.array([a11_0, a12_0, a21_0, a22_0])
                updates_to_plot = list(true_vals != initial_guesses)
                updates_on      = true_vals != initial_guesses

                # ------------ Algorithm parameters ---------------- (Moved out of loop)
                # relaxation time period (time to wait between updates)
                # time_between = np.array([relax_time] * 4) # creates a 1 x 4 array.
                time_between = self.get_relax_time_array()
                # Position, velocity thresholds for updates  --- (Moved out of loop)

                # Package mus, parms, thresholds for use in solve_ivp()
                mus = np.array([mu_1, mu_2, mu_3, mu_4])
                parms = np.array([a11, a12, a21, a22, a11_0, a12_0, a21_0, a22_0], dtype = np.float64)
                thresholds = np.array([u_thold, v_thold, ut_thold, vt_thold, d])

                # ------------ Simulation parameters ----------------
                sim_time = 100 # Stopping time
                self.set_dt(1/2 * (1 / self.get_mu_value()))
                t_span = [0, sim_time]
                t = np.arange(0, sim_time, self.get_dt())

                # ------------ Initialize system --------------------
                S0 = np.array([1, 1, 3, 3])                 # Initialize [x, y, xt, yt]
                guesses = [[a11_0],[a12_0],[a21_0],[a22_0]]  # Guesses

                try:
                    r1, r2, r3, r4 = self.calc_rhs(S0, mus, parms)
                    derivs = [[r1, r2, r3, r4]]
                    err = [[abs(S0[2] - S0[0]), abs(S0[3] - S0[1])]]   # Position error
                    last_updates     = np.array([0,0,0,0])             # Time of last updates
                    self.set_idx_last_updates(np.array([0,0,0,0]))     # Index of last updates
                    self.set_tfe([])                                   # To record times when model() is called

                    # Package into time_tracker arg for solve_ivp()
                    time_tracker = [last_updates, time_between, self.get_tfe()]
                    _args = (mus, parms, thresholds, derivs, guesses, time_tracker, updates_on, err, self.get_idx_last_updates())
                    sol, guesses, derivs = self.run_simulation(t_span, S0, t, true_vals, _args)
                    is_tr_det_graph_plotted = tr_det_graph.organize_data(guesses, true_vals)

                    if self.is_mtrx_usable() == True and is_tr_det_graph_plotted == True:
                        line_graph.init(guesses, true_vals, sol.t, i, 1e-5, self.get_solutions(), self.get_signal_err())
                        i += 1
                        gui_counter += 1
                    # if is_tr_det_graph_plotted  != False:
                        # sol_graph.init(sol, guesses, true_vals, i, 1e-5)

                    else:
                        print('\n', mtrx, '\n' + '\n', 'UNUSABLE MTRX', "\n")
                        print("CYCLE:", gui_counter, "- SKIP THIS ITERATION (1).")
                        if len(args) == 5: #i.e. If you're importing from a text file
                            i += 1
                            gui_counter += 1
                        else:
                            break

                except ValueError:
                    print('\n', mtrx, '\n' + '\n', 'UNUSABLE MTRX', "\n")
                    print("CYCLE:", gui_counter, "- SKIP THIS ITERATION (2).")
                    break

        print("**********************************************","END - SIMULATION" ,"***************************************************", '\n')
        line_graph.display_avg_rel_err_comp()
        tr_det_graph.display(ev_type, pp_type, loop_limit)
        plt.show()


    #2
    def run_simulation(self, t_span, S0, t, true_vals, args):
        mus, parms, thresholds, derivs, guesses, time_tracker, updates_on, err, idx_temp = args
        self.set_idx_last_updates(idx_temp)

        # ---------- Run simulation ------------------------
        start        = time.time()
        sol          = solve_ivp(self.model, t_span = t_span, y0 = S0, method ='BDF', t_eval = t, args = args)

        # ---------- Handle output ------------------------
        guesses      = np.array(guesses, dtype = np.float64)
        derivs       = np.array(derivs)

        # Reshape guesses and derivs arrays to match solution
        # This is necessary because solve_ivp() calls model() minimally to improve computational efficiency
        # but returns sol, an interpolation of the solution at all time points of t_eval. However, guesses
        # and derivs are only recorded when model() is called (which is why we recorded tfe :) )
        dt = self.get_dt()
        num_iter     = round(self.get_tfe_elt(-1) / dt)
        t_sol        = np.linspace(0, dt * num_iter, num = num_iter)
        f_eval       = np.searchsorted(self.get_tfe(), t_sol)
        guesses      = guesses[:, f_eval]
        derivs       = derivs[f_eval,:]
        a11s, a12s, a21s, a22s = guesses
        a11,  a12,  a21,  a22  = true_vals
        return sol, guesses, derivs



    #3
    def model(self, t, S, mus, parms, thresholds, derivs, guesses, time_tracker, updates_on, err, idx_last_updates):
        '''
        function called by odeint to return the time derivative of
        the system S at time t.

        INPUT

            mus              : NumPy array of shape (4) with values [mu_1, mu_2, mu_3, mu_4]
                               Note a full nudging matrix is not advised (typically choose
                               diagonal or off diagonal elements depending on which params
                               you want to recover)

            parms            : NumPy array of shape (8) with values [a11, a12, a21, a22, a11_t, a12_t, a21_t, a22_t]
                               where a11, a12, a21, a22 are the true parameter values and
                               a11_t, a12_t, a21_t, a22_t are the current guesses.

            thresholds       : NumPy array of shape(5) with values [u_thold,
                               ut_thold, v_thold, vt_thold, d] which are used to decide whether to
                               update parameters. Note these are not being used right now, updates
                               are being made any time the relaxation period has lapsed

            derivs           : list of time derivatives at all function evaluations

            guesses          : list of parameter guesses made at all function evals

            time_tracker     : list containing last_updates (the timestamp of the last parameter updates)
                               time_between (relaxation time)
                               tfe (times model() is called in simulation)

            updates_on       : boolean array specifying which parameters should be updated

            err              : local position errors (for get_tholds() -- not in use right now)

            idx_last_updates : index of last updates (for get_tholds() -- not in use right now)
        '''

        # Unpack args
        a11, a12, a21, a22, a11_t, a12_t, a21_t, a22_t = parms
        last_updates, time_between, tfe_temp = time_tracker
        self.set_relax_time_array(time_between)
        self.set_tfe(tfe_temp)
        curr_time = t
        time_tracker[-1].append(curr_time)      # Record time of function call in tfe
        u_now = abs(S[2] - S[0])                # Current x error |xt - x|
        v_now = abs(S[3] - S[1])                # Current y error |yt - y|
        err.append([u_now, v_now])
        pos_err = np.array(err).T
        time_since     = curr_time - last_updates
        # time_threshold = time_since > time_between # Check if relaxation time has lapsed for parameters
        avg_rel_err = self.get_elt_avg_rel_err(a21, a22, a21_t, a22_t)
        init_rlx_val = self.get_init_relax_time_value()
        if avg_rel_err > 1e-1:
            # self.set_mu_value(self.get_mu_value() * 20)
            self.set_relax_time_value(init_rlx_val * .1)
            # print("xl")
        elif avg_rel_err > 1e-3:
            # self.set_mu_value(self.get_mu_value() * 15)
            self.set_relax_time_value(init_rlx_val * .2)
            # print("l")
        elif avg_rel_err > 1e-5:
            # self.set_mu_value(self.get_mu_value() * 5)
            self.set_relax_time_value(init_rlx_val * 1.5)
            # print("m")
        elif avg_rel_err > 1e-8:
            # self.set_mu_value(self.get_mu_value() + self.get_mu_value())
            self.set_relax_time_value(init_rlx_val * 2)
            # print("s")

        # print("rta:", self.get_relax_time_array())
        time_threshold = time_since > self.get_relax_time_array() # Check if relaxation time has lapsed for parameters
        stop_update = 1e-10  # Threshold on relative position error for stopping updates

        if u_now / abs(S[0]) < stop_update and v_now / abs(S[1]) < stop_update:
          updates_on[0] = updates_on[1] = updates_on[2] = updates_on[3] = False

        # Update parms where relaxation period has lapsed and stop criteria not met
        to_update  = np.logical_and(time_threshold, updates_on)
        update_idx = [i for i, x in enumerate(to_update) if x]

        ##############################
        # avg_rel_err = self.get_elt_avg_rel_err(a21, a22, a21_t, a22_t)
        # if avg_rel_err >= 1e-5:
        #     self.set_mu_value(2.11e+08)

        # elif avg_rel_err > 1e-5 and self.get_mu_value() == 2.11e+08:
        #     while avg_rel_err > 1e-5:
        #         self.set_mu_value(self.get_mu_value() - 1)
                # print("still bad", self.get_mu_value())
        # else:
        #     print("good", self.get_mu_value())

        # mus[0] = mus[3] = self.get_mu_value()

        ##############################
        # avg_rel_err = self.get_elt_avg_rel_err(a21, a22, a21_t, a22_t)
        # if avg_rel_err >= 1e-5:
        #     self.set_mu_value(2.11e+08)
        # mus[0] = mus[3] = self.get_mu_value()
        ##############################


        self.update_parms(update_idx, S, mus, parms, guesses, thresholds, self.get_rule(), time_tracker, pos_err, self.get_idx_last_updates())
        no_update_idx = [i for i, x in enumerate(to_update) if not x]

        for i in no_update_idx:
            guesses[i].append(parms[i + 4])

        #New version
        r1, r2, r3, r4 = self.calc_rhs(S, mus, parms)
        St = [r1, r2, r3, r4]
        derivs.append(St)
        return St

    #4
    def update_parms(self, which, S, mus, parms, guesses, thresholds, rule, time_tracker, pos_err, ilast):
        '''
        Updates parameter guesses [if thresholds are met]
        '''
        curr_time = time_tracker[-1][-1]

        # Extract input
        x, y, xt, yt = S
        mu_1, mu_2, mu_3, mu_4 = mus

        a11, a12, a21, a22, a11_t, a12_t, a21_t, a22_t = parms
        u_thold, v_thold, ut_thold, vt_thold, d = thresholds

        # Calculate derivatives
        dx_dt, dy_dt, dxt_dt, dyt_dt = self.calc_rhs(S, mus, parms)

        if self.is_derivs_nan_or_inf(dx_dt, dy_dt, dxt_dt, dyt_dt) == True:
            print('Exit function: update_parms.')
            self.set_is_mtrx_usable(False)
            return

        else:
            derivs_now  = np.array([dx_dt, dy_dt, dxt_dt, dyt_dt])

            # Calculate error
            u      = xt - x
            v      = yt - y
            ut = dxt_dt - dx_dt
            vt = dyt_dt - dy_dt
            errors = np.array([u, v, ut, vt])

            if abs(u) / abs(x) <= u_thold and abs(v) / abs(y) <= v_thold:

                for i in which:
                    new_guess = self.update_formula(rule, S, errors, derivs_now, mus, parms)[i]
                    parms[i + 4] = new_guess
                    guesses[i].append(new_guess)
                    time_tracker[0][i] = curr_time
                    self.set_idx_last_updates_elt(i, pos_err.shape[1] - 1)

            else:
                for i in update_idx:
                    parm_idx = i + 4
                    guesses[i].append(parm_idx)

    def get_elt_avg_rel_err(self, a_true, b_true, a_guess, b_guess):
        a_rel_err = abs(a_true - a_guess ) / abs(a_true)
        b_rel_err = abs(b_true - b_guess) / abs(b_true)
        return abs(a_rel_err + b_rel_err) / 2

    #5
    def calc_rhs(self, S, mus, parms):
        '''
        Returns derivatives of x, y, xt, yt
        '''
        # Extract input
        x, y, xt, yt = S
        mu_1, mu_2, mu_3, mu_4 = mus
        a11, a12, a21, a22, a11_t, a12_t, a21_t, a22_t = parms
        r1 = a21 * x + a11 * y
        r2 = a12 * x + a22 * y
        r3 = a21_t * xt + a11_t * yt - mu_1 * (xt - x) - mu_2 * (yt - y)
        r4 = a12_t  * xt + a22_t * yt - mu_3 * (xt - x) - mu_4 * (yt - y)

        self.set_solutions(x, y, xt, yt)
        self.set_signal_err(x, y, xt, yt)
        return r1, r2, r3, r4

    #5.1
    def is_nan_or_inf(self, val):
        if math.isinf(val) == True:
            return True

        elif math.isnan(val) == True:
            return True

        else:
            return False

    #6
    def is_derivs_nan_or_inf(self, r1, r2, r3, r4):
        if self.is_nan_or_inf(r1) == False and self.is_nan_or_inf(r2) == False and self.is_nan_or_inf(r3) == False and self.is_nan_or_inf(r4) == False:
            return False

        else:
            return True


    #7
    def update_formula(self, rule, S, errors, derivs_now, mus, parms):
        '''
        Returns parameter update formulas
        '''
        # Handle input
        x, y, xt, yt = S
        u, v, ut, vt = errors
        dx_dt, dy_dt, dxt_dt, dyt_dt = derivs_now
        mu_1, mu_2, mu_3, mu_4 = mus
        a11, a12, a21, a22, a11_t, a12_t, a21_t, a22_t = parms
        # print(x)

        switcher = {
            'constant'  : [a11_t, a12_t, a21_t, a22_t],

            'exact'     : [(a21_t * xt + a11_t * yt - a21 * x - mu_1 * u - mu_2 * v - ut) / y,
                           (a12_t * xt  + a22_t * yt - a22 * y - mu_3 * u - mu_4 * v - vt) / x,
                           (a21_t * xt + a11_t * yt - a11 * y - mu_1 * u - mu_2 * v - ut) / x,
                           (a12_t * xt  + a22_t * yt - a12 * x  - mu_3 * u - mu_4 * v - vt) / y],

            'drop_deriv': [(a21_t * xt + a11_t * yt - a21 * x - mu_1 * u - mu_2 * v) / y,
                           (a12_t * xt  + a22_t * yt - a22 * y - mu_3 * v - mu_4 * v) / x,
                           (a21_t * xt + a11_t * yt - a11 * y - mu_1 * u - mu_2 * v) / x,
                           (a12_t * xt  + a22_t * yt - a12 * x  - mu_3 * v - mu_4 * v) / y]
        }
        # self.set_solutions(x, y, xt, yt)
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

    def get_tholds(self, pos_err, ilast):
        '''
        Determines update threshold on position error using log linear fit of data
        since last update

        NOTE: No thresholding is being used in current version of algorithm, but this
              is in place in case it needs to be used in future iterations.
        '''
        log_pos_err = np.log(pos_err[:, ilast:])
        x_cords = np.arange(log_pos_err.shape[1])
        y_coords = log_pos_err.T
        p = np.polyfit(x_cords, y_coords, 1)
        u_thold = np.exp(p[0, 0] * x_cords[-1] + p[1, 0])
        v_thold = np.exp(p[0, 1] * x_cords[-1] + p[1, 1])

        return u_thold, v_thold


    # GETTERS
    def get_mu_value(self):
        return self._mu_value

    def get_relax_time_array(self):
        return self._relax_time_array

    def get_relax_time_value(self):
        return self._relax_time_array[0]

    def get_init_relax_time_array(self):
        return self._init_relax_time_array

    def get_init_relax_time_value(self):
        return self._init_relax_time_array[0]

    def get_dt(self):
        return self._dt

    def get_idx_last_updates(self):
        return self._idx_last_updates

    def get_tfe(self):
        return self._tfe

    def get_tfe_elt(self, i):
        return self._tfe[i]

    def get_rule(self):
        return self._rule

    def get_init_mu_value(self):
        return self.init_mu_value

    def is_mtrx_usable(self):
        return self._is_mtrx_usable

    def get_solutions(self):
        return self._x_list, self._y_list, self._xt_list, self._yt_list

    def get_signal_err(self):
        return self._x_signal_err_list, self._y_signal_err_list


    # SETTERS
    def set_mu_value(self, val):
        self._mu_value = val

    def set_init_mu_value(self, val):
        self.init_mu_value = val

    def set_relax_time_value(self, val):
        self._relax_time_array  = np.array([val] * 4)

    def set_relax_time_array(self, np_array):
        self._relax_time_array = np_array

    def set_init_relax_time_value(self, val):
        self._init_relax_time_array  = np.array([val] * 4) # I.e. A (1 x 4) numpy array


    def set_dt(self, val):
        self._dt = val

    def set_idx_last_updates(self, list):
        self._idx_last_updates = list

    def set_tfe(self, val):
        self._tfe = val

    def set_idx_last_updates_elt(self, i, val):
        self._idx_last_updates[i] = val

    def set_rule(self, rule_to_use):
        # rules = {0: 'constant', 1: 'exact', 2: 'drop_deriv'}
        self._rule = self.rule_string(rule_to_use)

    def reset_solutions(self):
        self._x_list = []
        self._y_list = []
        self._xt_list = []
        self._yt_list = []

    def set_is_mtrx_usable(self, bool):
        self._is_mtrx_usable = bool


    # for display
    def set_solutions(self, x, y, xt, yt):
        self._x_list.append(x)
        self._y_list.append(y)
        self._xt_list.append(xt)
        self._yt_list.append(yt)

    def set_signal_err(self, x, y, xt, yt):
        self._x_signal_err_list.append(abs(xt - x) / x)
        self._y_signal_err_list.append(abs(yt - y) / y)

    def reset_signal_err(self):
        self._x_signal_err_list = []
        self._y_signal_err_list = []
