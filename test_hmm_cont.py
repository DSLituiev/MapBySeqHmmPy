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
        x = np.array([ 4 ,5, 8, 17, 30, 32, 44, 53]);
        # f = q/r;
        def binomial_emission(q, r):
            return binom.pmf(q, r, pop.f_vect)
        
        emission_handle = binomial_emission
    
        self.HMM = hmm_cont(pop, emission_handle(q,r), 0.01*x )    

    
    def test_calc_transition(self):
        self.HMM.calcT()
        stat_prob_propagation = np.dot(self.HMM.pop.Pstat, self.HMM.Transition)
        probab_error = np.linalg.norm(stat_prob_propagation - self.HMM.pop.Pstat)
        assert(probab_error < 1e-12)
        print('`calc_transition` works well on stationary distribution!')


if __name__ == "__main__":
    unittest.main() 