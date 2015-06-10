# -*- coding: utf-8 -*-
"""
Created on Sun Jun  7 14:18:25 2015

@author: dima
"""
import numpy as np
from scipy.stats import binom

def binomial_emission(pop, q, r):
    return binom.pmf(q, r, pop.f_vect)


def binomial_uniform_emission(pop, q, r, membership=None):
    p_bin = binom.pmf(q, r, pop.f_vect)
    membership = np.array(membership)
    if membership is None:
        return p_bin
    else:
        assert (membership.shape == q.shape), "dimension mismatch"

        out_emission = membership * p_bin + (1 - membership) * (np.sum(p_bin, axis=0) / pop.Np)
        assert (out_emission >= 0).all(), "negative values in the emission matrix!"
        # assert (out_emission <= 1).all()
        # assert ( abs(np.sum(out_emission , axis = 0)-1) < 1e-20).all()
        return out_emission