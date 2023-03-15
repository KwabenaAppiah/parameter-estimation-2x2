import math
import random

import numpy as np


class Matrix2x2:
    def __init__(self, low_bnd, high_bnd, ev_type, pp_type):
        self._ev_type = ev_type
        self._pp_type = pp_type
        self._matrix = []
        self.get_matrix(low_bnd, high_bnd, ev_type, pp_type)

    def __str__(self):
        return str(self._matrix)
    
    # For generating matrices populated with float values
    # def get_random_matrix(self, low_bnd, high_bnd):
    #  self._matrix = np.random.uniform(low_bnd, high_bnd, (2, 2))
    #  return self._matrix
    
    # For generating matrices populated with integer values
    def get_random_matrix(self, low_bnd, high_bnd):
        is_mtrx_ready = False

        while is_mtrx_ready == False:
            self._matrix = np.random.randint(low_bnd, high_bnd, (2, 2))
            # if self._matrix[1, 0] != 0 and self._matrix[1, 1] != 0: # i.e. if gamma or delta equals zero, find another mtrx.
            if self.get_element(1, 0) != 0 and self.get_element(1, 1) != 0:  # i.e. if gamma or delta equals zero, find another mtrx.
                is_mtrx_ready = True
                break

        return self._matrix



    def get_matrix(self, low_bnd, high_bnd, ev_type, pp_type):
        ########## Real distinct eigenvalues: T^2 - 4D > 0 Below the parabola ###########
        if ev_type == "rde":
            mtrx = self.get_matrix_rde(low_bnd, high_bnd, ev_type, pp_type)

        ########## Real distinct eigenvalues: T^2 - 4D = 0 Below the parabola ###########
        elif ev_type == "re":
            mtrx = self.get_matrix_re(low_bnd, high_bnd, ev_type, pp_type)

        ########## complex eigenvalues: T^2 - 4D < 0 Below the parabola ###########
        elif ev_type == "ce":
            mtrx = self.get_matrix_ce(low_bnd, high_bnd, ev_type, pp_type)

        else:
            print("ERROR:", ev_type, "is not a valid input.")
            quit()

        return mtrx

    # SETTERS

    def set_ev_type(self, ev_type):
        self._ev_type = ev_type

    def set_pp_type(self, pp_type):
        self._pp_type = pp_type

    #### GETTERS #####################################################

    def get_element(self, i, j):
        return self._matrix[i, j]

    def get_ev_type(self):
        return self._ev_type

    def get_pp_type(self):
        return self._pp_type

    def get_eignevals(self, mtrx):
        x, y = np.linalg.eig(mtrx)
        return x[0], x[1]

    def get_trace(self, mtrx):
        return np.trace(mtrx)

    def get_det(self, mtrx):
        return np.linalg.det(mtrx)

    def get_matrix_rde(self, low_bnd, high_bnd, ev_type, pp_type):
        mtrx = self.get_random_matrix(low_bnd, high_bnd)
        tr, det = self.get_trace(mtrx), self.get_det(mtrx)
        ev_1, ev_2 = self.get_eignevals(mtrx)
        ev_type_2 = "T\N{SUPERSCRIPT TWO} - 4D > 0 (Below the parabola)"
        rule = False

        # ***************** D < 0: ev_1 < 0 < ev_2: SADDLE ******************************************************
        if ev_type == "rde" and pp_type == "saddle":
            while rule != True:
                if (tr**2) - 4 * det > 0 and det < 0 and ev_1 < 0 and 0 < ev_2:
                    rule = True
                    break
                else:
                    mtrx = self.get_random_matrix(low_bnd, high_bnd)
                    ev_1, ev_2 = self.get_eignevals(mtrx)
                    tr, det = self.get_trace(mtrx), self.get_det(mtrx)

        # ***************** D > 0 and T > 0: 0 < ev_1 < ev_2: SOURCE *****************************************
        elif ev_type == "rde" and pp_type == "source":
            while rule != True:
                if (tr**2) - 4 * det > 0 and det > 0 and tr > 0 and 0 < ev_1 and ev_1 < ev_2:
                    rule = True
                    break

                else:
                    mtrx = self.get_random_matrix(low_bnd, high_bnd)
                    ev_1, ev_2 = self.get_eignevals(mtrx)
                    tr, det = self.get_trace(mtrx), self.get_det(mtrx)

        # ***************** D > 0 and T < 0: ev_1 < ev_2 < 0: SINK *********************************************
        elif ev_type == "rde" and pp_type == "sink":
            while rule != True:
                if (tr**2) - 4 * det > 0 and det > 0 and tr < 0 and ev_1 < ev_2 and ev_2 < 0:
                    rule = True
                    break

                else:
                    mtrx = self.get_random_matrix(low_bnd, high_bnd)
                    ev_1, ev_2 = self.get_eignevals(mtrx)
                    tr, det = self.get_trace(mtrx), self.get_det(mtrx)

        else:
            print("ERROR:", pp_type, "is not a valid input.")
            quit()

        return mtrx
        # **************************************************************************************************************************

    def get_matrix_re(self, low_bnd, high_bnd, ev_type, pp_type):
        mtrx = self.get_random_matrix(low_bnd, high_bnd)
        tr, det = self.get_trace(mtrx), self.get_det(mtrx)
        ev_1, ev_2 = self.get_eignevals(mtrx)
        ev_type_2 = "T\N{SUPERSCRIPT TWO} - 4D = 0"
        rule = False

        # ***************** T < 0: ev_1 < 0: SINK *****************************************************
        if ev_type == "re" and pp_type == "sink":
            while rule != True:
                if (tr**2) - 4 * det == 0 and tr < 0 and ev_1 < 0:
                    rule = True
                    break

                else:
                    mtrx = self.get_random_matrix(low_bnd, high_bnd)
                    ev_1, ev_2 = self.get_eignevals(mtrx)
                    tr, det = self.get_trace(mtrx), self.get_det(mtrx)

        # ***************** T > 0: ev_1 > 0: SOURCE *****************************************************
        elif ev_type == "re" and pp_type == "source":
            while rule != True:
                if (tr**2) - 4 * det == 0 and tr > 0 and ev_1 > 0:
                    rule = True
                    break

                else:
                    mtrx = self.get_random_matrix(low_bnd, high_bnd)
                    ev_1, ev_2 = self.get_eignevals(mtrx)
                    tr, det = self.get_trace(mtrx), self.get_det(mtrx)

        else:
            print("ERROR:", pp_type, "is not a valid input.")
            quit()

        return mtrx
        # *****************************************************************************************************************

    def get_matrix_ce(self, low_bnd, high_bnd, ev_type, pp_type):
        mtrx = self.get_random_matrix(low_bnd, high_bnd)
        tr, det = self.get_trace(mtrx), self.get_det(mtrx)
        ev_1, ev_2 = self.get_eignevals(mtrx)
        ev_type_2 = "T\N{SUPERSCRIPT TWO} - 4D < 0 (Above the parabola)"
        rule = False

        # ***************** T < 0 (Left half of the plane): SPIRAL SINK ******************************************************
        if ev_type == "ce" and pp_type == "sp_sink":
            while rule != True:
                if (tr**2) - 4 * det < 0 and tr < 0:
                    rule = True
                    break

                else:
                    mtrx = self.get_random_matrix(low_bnd, high_bnd)
                    ev_1, ev_2 = self.get_eignevals(mtrx)
                    tr, det = self.get_trace(mtrx), self.get_det(mtrx)

        # ***************** T > 0 (Right half of the plane): SPIRAL SOURCE ******************************************************
        elif ev_type == "ce" and pp_type == "sp_source":
            while rule != True:
                if (tr**2) - 4 * det < 0 and tr > 0:
                    rule = True
                    break

                else:
                    mtrx = self.get_random_matrix(low_bnd, high_bnd)
                    ev_1, ev_2 = self.get_eignevals(mtrx)
                    tr, det = self.get_trace(mtrx), self.get_det(mtrx)

        # ***************** T = 0 (D - Axis): CENTER *****************************************************************************
        elif ev_type == "ce" and pp_type == "center":
            while rule != True:
                if (tr**2) - 4 * det < 0 and tr == 0:
                    rule = True
                    break

                else:
                    mtrx = self.get_random_matrix(low_bnd, high_bnd)
                    ev_1, ev_2 = self.get_eignevals(mtrx)
                    tr, det = self.get_trace(mtrx), self.get_det(mtrx)

        else:
            print("ERROR:", pp_type, "is not a valid input.")
            quit()

        return mtrx
        # **************************************************************************************************************************

    def print_ev_pp_dvdr(self, ev_type, pp_type):
        print("**************************", ev_type.upper(), "-", pp_type.upper(),"*****************************************", "\n")
