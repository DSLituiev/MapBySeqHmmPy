__author__ = 'dmytro lituiev'

import numpy as np
import numpy.random as rn
from population import *

def calcCumStateChange(dx_v, mut_v, tind):
    "multiply the state-derivative and mutation vectors:"
    DxMatr = dx_v[:,np.newaxis] * mut_v
    "sum up the states changes in the order they occur:"
    cumdx = np.cumsum( DxMatr.ravel()[tind] );
    return cumdx

class sim_pool_seq:

    def __init__(self, population_size, reads_expected, n_sampling_points, 
                 t_causative = None,  t_max = 1, max_rec_points_rate = 20, false_positives = 0):
        self.t_causative = t_causative
        self.hidstates = population(population_size)
        self.reads_expected = reads_expected
        self.n_sampling_points = n_sampling_points
        self.selection = (not isinstance(t_causative,type(None)))
        self.t_max = t_max
        self.max_rec_points = max_rec_points_rate * self.t_max
        self.false_positives = false_positives

    def _state_deriv_matrix_(self):
        dx_v = np.ones(self.max_rec_points)
        dx_v[1::2] = -1
        return dx_v

    def _generate_rec_times_(self):
        propensity = 1
        waiting_times = - np.log(rn.rand(self.max_rec_points, self.hidstates.N))
        t_matr = 1/ propensity * np.cumsum( waiting_times, axis= 0 )

        tind = np.argsort(t_matr.ravel())
        self.tchr = t_matr.ravel()[tind]
        self.tchr = self.tchr[self.tchr <= self.t_max]
        return t_matr, tind

    def generate_path(self):
        t_matr, tind  = self._generate_rec_times_()
        dx_v = self._state_deriv_matrix_()

        if self.selection:
            self.i_causative = np.argmin(np.abs(self.tchr - self.t_causative))
            self.t_causative = self.tchr[self.i_causative]
            "indt0"
            ic_causative = np.empty(t_matr.shape[1], dtype = int)
            for ii in range(0, t_matr.shape[1]):
                ic_causative[ii] = np.nonzero(t_matr[:,ii] > self.t_causative)[0][0]
            
            first_mut_v = -dx_v[ic_causative]
            "include the false positives"
            first_mut_v[ 1:(self.false_positives) ] = - first_mut_v[ 1:(self.false_positives) ];
            cumdx = calcCumStateChange(dx_v, first_mut_v, tind);
            "offset cumdx(indCausativeSNP)"
            self.kChr = self.hidstates.N - self.false_positives - cumdx[self.i_causative] + cumdx;



