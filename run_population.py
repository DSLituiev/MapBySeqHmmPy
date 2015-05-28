# -*- coding: utf-8 -*-
"""
Created on Tue Nov  4 10:31:05 2014
@author: Dmytro Lituiev
"""
from population import *
from hmm_cont import *
from scipy.stats import binom
import numpy as np

N = 5;

pop = population(N);

assert np.linalg.norm(np.dot(pop.Pstat, pop.Q)) < 1e-12, \
    "stationary vector must be in the null-space of the matrix Q!"

q = np.array([1, 2, 2, 1, 3, 4, 1, 3]);
r = np.array([4, 6, 3, 8, 6, 7, 5, 9]);
x = np.array([4, 5, 8, 17, 30, 32, 44, 53]);
f = q / r;


def binomial_emission(q, r):
    return binom.pmf(q, r, pop.f_vect)


emission_handle = binomial_emission

HMM = hmm_cont(pop, emission_handle(q, r), 0.01 * x)


def test_calc_transition(HMM):
    HMM.calcT()
    stat_prob_propagation = np.dot(HMM.hidstates.Pstat, HMM.Transition)
    probab_error = np.linalg.norm(stat_prob_propagation - HMM.hidstates.Pstat)
    assert (probab_error < 1e-12)
    print('calc_transition works well!')


test_calc_transition(HMM)

HMM.cumMatr()

HMM._runFB_()

HMM.xk_P_flat

(x_p_stat, xk_p_stat) = HMM.getLikelihoodOfAModel(pop.Pstat)

assert (abs(x_p_stat + np.mean(x_p_stat)) < 5e-2).all(), \
"the stationary distribution is distorted too much"

x_p_flat, junk = HMM.getLikelihoodOfAModel(1 / pop.Np)


