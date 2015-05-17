# -*- coding: utf-8 -*-
"""
Created on Fri Nov  7 19:55:36 2014
@author: Dmitri
"""

import unittest
from population import *
from hmm_cont import *
from scipy.stats import binom
import numpy as np


class hmm_cont_test(unittest.TestCase):
    def setUp(self):
        N = 5;
        pop = population(N);

        q = np.array([1, 2, 2, 1, 3, 4, 1, 3]);
        r = np.array([4, 6, 3, 8, 6, 7, 5, 9]);
        x = np.array([4, 5, 8, 17, 30, 32, 44, 53]);
        # f = q/r;
        def binomial_emission(q, r):
            return binom.pmf(q, r, pop.f_vect)

        emission_handle = binomial_emission

        self.HMM = hmm_cont(pop, emission_handle(q, r), 0.01 * x)
        print('setup completed')

    def test_calc_transition(self):
        self.HMM.calcT()
        stat_prob_propagation = np.dot(self.HMM.hidstates.Pstat, self.HMM.Transition)
        probab_error = np.linalg.norm(stat_prob_propagation - self.HMM.hidstates.Pstat, np.inf)
        assert (probab_error < 1e-12)
        print('`calc_transition` works well on stationary distribution!')
        print('error: %2.2e' % probab_error)

    def test_symmetry(self):
        x = np.arange(0,21)
        M = len(x)
        
        r = [21,] *M
        q = [10,] *M
        pop = population(5)
        p_bin = binom.pmf(q, r, pop.f_vect)
        
        a = hmm_cont(pop, p_bin, x)
        
        (xL, xkL) = a.getLikelihoodOfAModel(pop.Psel)
        
        diffr = (xL - xL[::-1])
        assert all( diffr < 1e-12 )
        print( 'max difference: %2.2e' % (max(diffr)) )
        
if __name__ == "__main__":
    unittest.main(exit=False)