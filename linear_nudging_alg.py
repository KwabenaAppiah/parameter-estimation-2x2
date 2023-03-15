import math
import random
import sys
import time

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

# lorenz nudge, project-specific
from tr_det_graph import TrDetGraph
from line_graph_indv import LineGraphIndv
from line_graph_comp import LineGraphComp
from matrix_2x2 import Matrix2x2

class LinearNudgingAlg:
    def __init__(self, ev_type, pp_type, mu_val, relax_time, bnd_val, loop_limit):
        low_bnd, high_bnd = -bnd_val, bnd_val
        self.set_rule(2)

        # ------------ Algorithm parameters ----------------
        # Nudging parameters
        # NOTE: Full nudging matrix not advised. Diagonal entries seem to work best but
        #       could play around with off diagonal entries. I typically use multiple of
        #       identity matrix
        mu_1 = mu_4 = mu_val
        mu_2 = mu_3 = 0

        # Position, velocity thresholds for updates
        # NOTE: Right now we are not using any thresholds, they are all set to infinity
        #       with no decay
        u_thold    = np.inf
        v_thold    = np.inf
        ut_thold   = np.inf
        vt_thold   = np.inf
        d          = 10         # Threshold decay factor
        step_val   = 10         # The step value for a22 and a21
        tr_det_graph = TrDetGraph(ev_type, pp_type, loop_limit)
        line_graph_indv = LineGraphIndv(ev_type, pp_type)
        line_graph_comp = LineGraphComp(ev_type, pp_type)


        print("**********************************************", "START - SIMULATION" ,"***************************************************", '\n')
        i = 0
        gui_counter = i + 1
        while i < loop_limit:
            while i != loop_limit:
                mtrx = Matrix2x2(low_bnd, high_bnd, ev_type, pp_type)

                if mtrx.get_element(1, 0) == 0 and mtrx.get_element(1, 1) == 0:
                     print("CYCLE:", gui_counter, " - SKIP ITERATION", 1)
                     break

                else:
                    print("**********************************", "CYCLE:", gui_counter, "/" , loop_limit, "*************************************************************************", '\n')
                    print(ev_type.upper(), "-", pp_type.upper())
                    a11, a12, a21, a22 = mtrx.get_element(0, 0), mtrx.get_element(0, 1), mtrx.get_element(1, 0), mtrx.get_element(1, 1)
                    true_vals = np.array([a11, a12, a21, a22])

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

                    # relax_timeation period (time to wait between updates)
                    time_between = np.array([relax_time] * 4)

                    # Position, velocity thresholds for updates  --- (Moved out of loop)

                    # Package mus, parms, thresholds for use in solve_ivp()
                    mus = np.array([mu_1, mu_2, mu_3, mu_4])
                    parms = np.array([a11, a12, a21, a22, a11_0, a12_0, a21_0, a22_0], dtype = np.float64)
                    thresholds = np.array([u_thold, v_thold, ut_thold, vt_thold, d])

                    # ------------ Simulation parameters ----------------
                    sim_time = 100 # Stopping time
                    # dt = 0.0001  # Timestep
                    self.set_dt(0.0001)
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
                        args = (mus, parms, thresholds, derivs, guesses, time_tracker, updates_on, err, self.get_idx_last_updates())
                        sol, guesses, derivs = self.run_simulation(t_span, S0, t, true_vals, args)
                        # line_graph_indv_flag = line_graph_indv.plot_graph(guesses, true_vals, i)
                        # line_graph_comp_flag = line_graph_comp.organize_data(guesses, true_vals, i)
                        line_graph_indv.plot_graph(guesses, true_vals, i)
                        line_graph_comp.organize_data(guesses, true_vals, i)
                        tr_det_graph_plotted = tr_det_graph.organize_data(guesses, true_vals)

                        if tr_det_graph_plotted  != False:
                            i += 1
                            gui_counter += 1
                        else:
                            print('\n', mtrx, '\n' + '\n', 'UNUSABLE MTRX', "\n")
                            print("CYCLE:", gui_counter, "- SKIP THIS ITERATION.")
                            break

                    except ValueError:
                        print('\n', mtrx, '\n' + '\n', 'UNUSABLE MTRX', "\n")
                        print("CYCLE:", gui_counter, "- SKIP THIS ITERATION.")
                        break

        print("**********************************************","END - SIMULATION" ,"***************************************************", '\n')
        line_graph_comp.display()
        tr_det_graph.display(ev_type, pp_type, loop_limit)
        plt.show()

    # GETTERS
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


    # SETTERS
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


    def calc_rhs(self, S, mus, parms):
        '''
        Returns derivatives of x, y, xt, yt
        '''
        # Extract input
        x, y, xt, yt = S
        mu_1, mu_2, mu_3, mu_4 = mus
        a11, a12, a21, a22, a11_t, a12_t, a21_t, a22_t = parms
        # Replaced: alpha = a11, beta = a12, delta = a21, gamma = a22

        r1 = a21 * x + a11 * y
        r2 = a12 * x + a22 * y
        r3 = a21_t * xt + a11_t * yt - mu_1 * (xt - x) - mu_2 * (yt - y)
        r4 = a12_t  * xt + a22_t * yt - mu_3 * (xt - x) - mu_4 * (yt - y)

        return r1, r2, r3, r4


    def rule_string(self, which):
        '''
        Converts rule index to string for retrieving update formula from dictionary
        '''
        switcher = {0: 'constant', 1: 'exact', 2: 'drop_deriv'}
        return switcher.get(which)


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

        return switcher.get(rule, "Update rule not recognized")



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

        if self.nan_or_inf_derivs(dx_dt, dy_dt, dxt_dt, dyt_dt ) != False:
            print('Exit function: update_parms.')
            return

        else:
            derivs_now  = np.array([dx_dt, dy_dt, dxt_dt, dyt_dt])
            # Calculate error
            u      = xt - x
            v      = yt - y
            ut     = dxt_dt - dx_dt
            vt     = dyt_dt - dy_dt
            errors = np.array([u, v, ut, vt])

            if abs(u) / abs(x) <= u_thold and abs(v) / abs(y) <= v_thold:

                for i in which:
                    new = self.update_formula(rule, S, errors, derivs_now, mus, parms)[i]
                    parms[i + 4] = new
                    guesses[i].append(new)
                    time_tracker[0][i] = curr_time
                    self.set_idx_last_updates_elt(i, pos_err.shape[1] - 1)

            else:
                for i in update_idx:
                    parm_idx = i + 4
                    guesses[i].append(parm_idx)


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
        self.set_tfe(tfe_temp)

        curr_time = t
        time_tracker[-1].append(curr_time)  # Record time of function call in tfe
        u_now = abs(S[2] - S[0])               # Current x error |xt - x|
        v_now = abs(S[3] - S[1])               # Current y error |yt - y|
        err.append([u_now, v_now])
        pos_err = np.array(err).T
        time_since     = curr_time - last_updates
        time_threshold = time_since > time_between # Check if relaxation time has lapsed for parameters

        stop_update = 1e-10  # Threshold on relative position error for stopping updates

        if u_now / abs(S[0]) < stop_update and v_now / abs(S[1]) < stop_update:
          updates_on[0] = updates_on[1] = updates_on[2] = updates_on[3] = False

        # Update parms where relaxation period has lapsed and stop criteria not met
        to_update  = np.logical_and(time_threshold, updates_on)

        update_idx = [i for i, x in enumerate(to_update) if x]
        # update_parms(update_idx, S, mus, parms, guesses, thresholds, rule, time_tracker, pos_err, self.idx_last_updates)
        #self.update_parms(update_idx, S, mus, parms, guesses, thresholds, 'drop_deriv', time_tracker, pos_err, self.get_idx_last_updates())
        self.update_parms(update_idx, S, mus, parms, guesses, thresholds, self.get_rule(), time_tracker, pos_err, self.get_idx_last_updates())
        no_update_idx = [i for i, x in enumerate(to_update) if not x]

        for i in no_update_idx:
            guesses[i].append(parms[i + 4])

        #New version
        r1, r2, r3, r4 = self.calc_rhs(S, mus, parms)
        St = [r1, r2, r3, r4]
        derivs.append(St)
        return St

        # St = calc_rhs(S, mus, parms)
        # derivs.append(St)
        # return St

    def nan_or_inf_derivs(self, r1, r2, r3, r4):
        if math.isinf(r1) != False and math.isinf(r2) != False and math.isinf(r3) != False and math.isinf(r4) != False:
            return True

        elif math.isnan(r1) != False and math.isnan(r2) != False and math.isnan(r3) != False and math.isnan(r4) != False:
            return True

        return False


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

        #print("Runtime: {:.4f} seconds.\n".format(time.time() - start))
        #print("Final a11 abs. error: {:.4e}\n".format(abs(a11 - a11s)[-1]))
        #print("Final a12 abs. error: {:.4e}\n".format(abs(a12 - a12s)[-1]))
        #print("Final a21 abs. error: {:.4e}\n".format(abs(a21 - a21s)[-1]))
        #print("Final a22 abs. error: {:.4e}\n".format(abs(a22 - a22s)[-1]))
        return sol, guesses, derivs
