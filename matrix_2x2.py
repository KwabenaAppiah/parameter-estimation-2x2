import math
import random
import numpy as np


class Matrix2x2:
    def __init__(self, low_bnd, high_bnd, ev_type, pp_type):
        self._ev_type = ev_type
        self._pp_type = pp_type
        self._matrix = []
        self._mtrx_sample_space  = []
        self.init_matrix(low_bnd, high_bnd, ev_type, pp_type)

    def __str__(self):
        return str(self._matrix)

    def init_matrix(self, low_bnd, high_bnd, ev_type, pp_type):
        ########## Real distinct eigenvalues: T^2 - 4D > 0 Below the parabola ###########
        if ev_type == "rde":
            self.set_matrix_rde(low_bnd, high_bnd, ev_type, pp_type)

        ########## Real distinct eigenvalues: T^2 - 4D = 0 Below the parabola ###########
        elif ev_type == "re":
            self.set_matrix_re(low_bnd, high_bnd, ev_type, pp_type)

        ########## complex eigenvalues: T^2 - 4D < 0 Below the parabola ###########
        elif ev_type == "ce":
            self.set_matrix_ce(low_bnd, high_bnd, ev_type, pp_type)

        else:
            print("ERROR:", ev_type, "is not a valid input.")
            quit()


    def get_matrix(self):
        return self._matrix

    def get_mtrx_sample_space(self):
        return self._mtrx_sample_space

    def get_mtrx_sample_space_elt(self, index):
        return self._mtrx_sample_space[index]

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

    def get_each_elt(self):
        return self._matrix[0: 0], self._matrix[0: 1], self._matrix[1: 0], self._matrix[1: 1]

    def set_matrix(self, mtrx):
        self._matrix = mtrx

    def set_mtrx_sample_space(self, new_mtrx):
        self._mtrx_sample_space.append(new_mtrx)

    #### SETTERS #####################################################
    def set_ev_type(self, ev_type):
        self._ev_type = ev_type

    def set_pp_type(self, pp_type):
        self._pp_type = pp_type

    # For generating matrices populated with float values
    def set_random_matrix(self, low_bnd, high_bnd):
        self.set_matrix(np.random.uniform(low_bnd, high_bnd, (2, 2)))


    def set_matrix_rde(self, low_bnd, high_bnd, ev_type, pp_type):
        self.set_random_matrix(low_bnd, high_bnd)
        mtrx = self.get_matrix()
        tr, det = self.get_trace(mtrx), self.get_det(mtrx)
        ev_1, ev_2 = self.get_eignevals(mtrx)
        has_mtrx_been_found = False

        if pp_type != "saddle" and pp_type != "source" and pp_type != "sink":
            print("ERROR:", pp_type, "is not a valid input.")
            quit()

        while has_mtrx_been_found != True:
            # ***************** D < 0: ev_1 < 0 < ev_2: SADDLE ******************************************************
            if pp_type == "saddle" and (tr**2 - 4*det) > 0 and det < 0 and ev_1 < 0 < ev_2:
                has_mtrx_been_found = True
                self.set_matrix(mtrx)
                break
            # ***************** D > 0 and T > 0: 0 < ev_1 < ev_2: SOURCE *****************************************
            elif pp_type == "source" and (tr**2 - 4*det) > 0 and det > 0 and tr > 0 and 0 < ev_1 < ev_2:
                has_mtrx_been_found = True
                self.set_matrix(mtrx)
                break

            # ***************** D > 0 and T < 0: ev_1 < ev_2 < 0: SINK *********************************************
            elif pp_type == "sink" and (tr**2 - 4*det) > 0 and det > 0 and tr < 0 and ev_1 < ev_2 < 0:
                has_mtrx_been_found = True
                self.set_matrix(mtrx)
                break

            else:
                self.set_random_matrix(low_bnd, high_bnd)
                mtrx = self.get_matrix()
                tr, det = self.get_trace(mtrx), self.get_det(mtrx)
                ev_1, ev_2 = self.get_eignevals(mtrx)



    def set_matrix_re(self, low_bnd, high_bnd, ev_type, pp_type):

        if pp_type != "sink" and pp_type != "source":
            print("ERROR:", pp_type, "is not a valid input.")
            quit()

        sample_space = np.arange(low_bnd, high_bnd + 1, 0.5)
        
        for a11 in sample_space:
            for a12 in sample_space:
                for a21 in sample_space:
                    for a22 in sample_space:
                        temp_mtrx = np.array([[a11, a12], [a21, a22]])
                        tr = a11 + a22
                        det = (a11 * a22) - (a12 * a21)

                        if pp_type == "sink" and (tr**2 - 4*det) == 0 and tr < 0:
                            self.set_mtrx_sample_space(temp_mtrx)

                        elif pp_type == "source" and (tr**2 - 4*det) == 0 and tr > 0:
                            self.set_mtrx_sample_space(temp_mtrx)

        random_index = random.randint(0, len(self.get_mtrx_sample_space()))
        mtrx = self.get_mtrx_sample_space_elt(random_index)
        self.set_matrix(mtrx)


    def set_matrix_ce(self, low_bnd, high_bnd, ev_type, pp_type):
        self.set_random_matrix(low_bnd, high_bnd)
        mtrx = self.get_matrix()
        tr, det = self.get_trace(mtrx), self.get_det(mtrx)
        ev_1, ev_2 = self.get_eignevals(mtrx)
        rule = False

        if pp_type == "center":
            self.set_matrix_ce_center(low_bnd, high_bnd, ev_type, pp_type)

        elif pp_type == "sp_sink" or pp_type == "sp_source":
            while rule != True:
                # ***************** T < 0 (Left half of the plane): SPIRAL SINK ******************************************************
                if pp_type == "sp_sink" and (tr**2) - 4 * det < 0 and tr < 0:
                    rule = True
                    self.set_matrix(mtrx)
                    break

                # ***************** T > 0 (Right half of the plane): SPIRAL SOURCE ******************************************************
                elif pp_type == "sp_source" and (tr**2) - 4 * det < 0 and tr > 0:
                    rule = True
                    self.set_matrix(mtrx)
                    break

                # ***************** T = 0 (D - Axis): CENTER *****************************************************************************
                else:
                    self.set_random_matrix(low_bnd, high_bnd)
                    mtrx = self.get_matrix()
                    tr, det = self.get_trace(mtrx), self.get_det(mtrx)
                    ev_1, ev_2 = self.get_eignevals(mtrx)

        else:
            print("ERROR:", pp_type, "is not a valid input.")
            quit()


    def set_matrix_ce_center(self, low_bnd, high_bnd, ev_type, pp_type):
        # sample_space = np.linspace(low_bnd, high_bnd + 1, 1000000)
        sample_space = np.arange(low_bnd, high_bnd + 1, 0.5)
        for a11 in sample_space:
            for a12 in sample_space:
                for a21 in sample_space:
                    for a22 in sample_space:
                        #temp_mtrx = np.array([[a11, a12], [a21, a22]])
                        tr = a11 + a22
                        det = (a11 * a22) - (a12 * a21)
                        if (tr ** 2) - 4 * det < 0 and tr == 0:
                            self.set_mtrx_sample_space(temp_mtrx)

        random_index = random.randint(0, len(self.get_mtrx_sample_space()))
        mtrx = self.get_mtrx_sample_space_elt(random_index)
        self.set_matrix(mtrx)


    def set_matrix_ce_center(self, low_bnd, high_bnd, ev_type, pp_type):
        # sample_space = np.linspace(low_bnd, high_bnd + 1, 1000000)
        sample_space = np.arange(low_bnd, high_bnd + 1, 0.5)
        for a11 in sample_space:
            for a12 in sample_space:
                for a21 in sample_space:
                    for a22 in sample_space:
                        tr = a11 + a22
                        det = (a11 * a22) - (a12 * a21)
                        if (tr**2 - 4*det) < 0 and tr == 0:
                            self.set_mtrx_sample_space(temp_mtrx)

        random_index = random.randint(0, len(self.get_mtrx_sample_space()))
        mtrx = self.get_mtrx_sample_space_elt(random_index)
        self.set_matrix(mtrx)


    # def set_matrix_ce_center(self, low_bnd, high_bnd, ev_type, pp_type):
    #     sample_space = np.arange(low_bnd, high_bnd + 1, 0.1)
    #
    #     A = np.meshgrid(sample_space, sample_space, sample_space, sample_space)
    #     a11, a12, a21, a22 = A
    #     tr = a11 + a22
    #     det = a11 * a22 - a12 * a21
    #
    #     valid_indices = (tr**2 - 4*det < 0) & (tr == 0)
    #     valid_mtrx = [(a11[i], a12[i], a21[i], a22[i]) for i in np.where(valid_indices)[0]]
    #     print(valid_mtrx)
    #
    #     mtrx = random.choice(valid_mtrx)
    #     self.set_matrix(mtrx)


    def print_ev_pp_dvdr(self, ev_type, pp_type):
        print("**************************", ev_type.upper(), "-", pp_type.upper(),"*****************************************", "\n")
