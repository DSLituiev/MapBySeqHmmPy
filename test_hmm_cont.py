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
    absTolerance = 1e-2
    relTolerance = 1e-12
    
    def setUp(self):
        N = 5;
        pop = population(N);

        q = np.array([1, 3, 3, 2, 1, 2, 4, 1]);
        r = np.array([4, 9, 6, 3, 8, 6, 7, 5]);
        x = np.array([4, 5, 8, 17, 30, 32, 44, 53]);
        q = q[:4]
        r = r[:4]
        x = x[:4]
        # f = q/r;
        def binomial_emission(q, r):
            return binom.pmf(q, r, pop.f_vect)

        emission_handle = binomial_emission

        self.HMM = hmm_cont(pop, emission_handle(q, r), 0.01 * x)
        print('setup completed')

    def test_calc_transition_propagation(self):
        print('----------------------------------------------------------------------')
        self.HMM.calcT()
        stat_prob_propagation = np.dot(self.HMM.hidstates.Pstat, self.HMM.Transition)
        probab_error = np.linalg.norm(stat_prob_propagation - self.HMM.hidstates.Pstat, np.inf)
        assert (probab_error < self.relTolerance)
        print('`calc_transition` works well on stationary distribution!')
        print('error: %2.2e' % probab_error)
        
    def test_manual_propagation(self):
        print('----------------------------------------------------------------------')
        m = np.matrix
        mdiag = lambda x: m(np.diag(x))
        P_ = m(np.diag(self.HMM.hidstates.Pstat))
        self.HMM.calcT()
        
        out = [0,]* 4
        out[0] = m(np.ones(self.HMM.Np)) * \
        P_ * \
        mdiag(self.HMM.E[:,0]) * m( self.HMM.Transition[0] ) * \
        mdiag(self.HMM.E[:,1]) * m( self.HMM.Transition[1] ) * \
        mdiag(self.HMM.E[:,2]) * m( self.HMM.Transition[2] ) * \
        m(self.HMM.E[:,3] ).T
        
        out[1] = m(np.ones(self.HMM.Np)) * \
        mdiag(self.HMM.E[:,0] ) * m(self.HMM.Transition[0]).T * \
        P_ * \
        mdiag(self.HMM.E[:,1] ) * m( self.HMM.Transition[1] ) * \
        mdiag(self.HMM.E[:,2] ) * m( self.HMM.Transition[2] ) * \
        m(self.HMM.E[:,3] ).T
        
        out[2] = m(np.ones(self.HMM.Np)) * \
        mdiag(self.HMM.E[:,0] ) * m(self.HMM.Transition[0]).T * \
        mdiag(self.HMM.E[:,1] ) * m( self.HMM.Transition[1] ).T * \
        P_ * \
        mdiag(self.HMM.E[:,2] ) * m( self.HMM.Transition[2] ) * \
        m(self.HMM.E[:,3] ).T
       
        out[3] = m(np.ones(self.HMM.Np)) * \
        mdiag(self.HMM.E[:,0] )  * m(self.HMM.Transition[0] ).T * \
        mdiag(self.HMM.E[:,1] ) * m( self.HMM.Transition[1] ).T * \
        mdiag(self.HMM.E[:,2] ) * m( self.HMM.Transition[2] ).T * \
        P_ * m(self.HMM.E[:,3] ).T
        
#        self.HMM.cumMatr()
        
#        print( 'diff Beta: ', \
#        m(np.ones(self.HMM.Np)) * mdiag(self.HMM.E[:,0] ) * m(self.HMM.Transition[0]).T - \
#        m( 10**self.HMM.logBeta[:,1] ) )
        
#        print( 'diff Alpha [M]: ', \
#        m(self.HMM.E[:,3] ).T - m( 10**self.HMM.logAlpha[:,-1] ).T \
#        )
        
#        print( 'diff Alpha: ', \
#        mdiag(self.HMM.E[:,2] ) * m( self.HMM.Transition[2] ) * m(self.HMM.E[:,3] ).T - \
#        m( 10**self.HMM.logAlpha[:,-2] ).T )
        
        mean_lh =  np.mean(out)
        print('manual calculation [mean]: %f' % mean_lh )
        
        for n in range(0,4):
            print( '%u: %f' % (n, np.log10(out[n])) )
        assert ( out[1] - out[0] ) / mean_lh < self.relTolerance and \
        ( out[3] - out[0] ) / mean_lh < self.relTolerance , \
        "the stationary distribution is too different at different points"
        
        return
    
    def test_calc_stationary(self):
        print('----------------------------------------------------------------------')
        (x_p_stat, xk_p_stat) = self.HMM.get_model_likelihood(self.HMM.hidstates.Pstat)
        
        print('number of loci / points %u' %  len(x_p_stat))
        print('stationary model log-lh:')
        print(x_p_stat)
        print('stationary model lh:')
        print(10** x_p_stat)
        assert (abs(x_p_stat + np.mean(x_p_stat)) / np.mean(x_p_stat) < self.relTolerance ).all(), \
        "the stationary distribution is too different at different points"

    def test_symmetry(self):
        print('----------------------------------------------------------------------')
        x = np.arange(0,4)
        M = len(x)
        
        r = [21,] *M
        q = [10,] *M
        pop = population(5)
        p_bin = binom.pmf(q, r, pop.f_vect)
        
        a = hmm_cont(pop, p_bin, x)
        (xL, xkL) = a.get_model_likelihood(pop.Psel)
        
        print('selection model lh:')
        print( xL )
        
        diffr = (xL - xL[::-1])
        assert all( diffr / np.mean(xL) < self.relTolerance )
        print( 'max difference: %2.2e' % (max(diffr)) )
        
if __name__ == "__main__":
    unittest.main(exit=False)