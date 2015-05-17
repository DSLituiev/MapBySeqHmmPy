# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 15:20:27 2015

@author: dima
"""
import numpy as np
import numpy.random as rn
from scipy.stats import binom
from hmm_cont import *
from population import *
from simulator import *

from unmix_repeats import *

import matplotlib.pyplot as plt
###############################################################################

s = sim_pool_seq(population_size=4, reads_expected=30, n_sampling_points=20,
                 t_causative=0.5)
s.generate_path()
###############################################################################
VERBOSE = False
x = np.arange(0,21)
x0 = 3
x1 = 16
dx0 = 0.1
dx1 = 0.05
x = np.concatenate( (x[0:x0] , x0 + np.arange(0, 1, dx0), 
                     x[x0+1:x1], x1 + np.arange(0, 1, dx1), x[x1+1:]) )

dx, log10_dx = min_distances(x)
membership, mixtObj, iTrue, out_mode_number, _ = unmix(log10_dx)
    
M = len(x)
II = 500
ixL = np.empty((II, M))

for ii in range(0, II):
#    r = [21,] *M
#    q = [10,] *M
    r = rn.poisson(20, (M,) )
    q = rn.binomial(r, 0.25 )
    
    if VERBOSE:
        print('iteration %u' % ii)
        print('r: ' , end = '')
        print(r)
        
        print('q: ' , end = '')
        print(q)
    
    pop = population(5)
    p_bin = binom.pmf(q, r, pop.f_vect)
    
    a = hmm_cont(pop, p_bin, x)
    
    (xL, xkL) = a.getLikelihoodOfAModel(pop.Psel)
    ixL[ii,:] = xL


import matplotlib.pyplot as plt

###############################################################################

plt.figure()
plt.plot(x, len(x) * [np.mean(ixL[:,1:-2][:,dx[1:-2] == 1].ravel() ),] , 'g-', 
    linewidth = 0.5, label = 'mean of all')
plt.plot(x, np.mean(ixL, 0), 'r-', label = 'mean')
plt.plot(x, np.min(ixL, 0), 'k-', linewidth = 0.5, label = 'min')
plt.plot(x, np.max(ixL, 0), 'k-', linewidth = 0.5, label = 'max')
plt.plot(x, xL, 'b.-', linewidth = 0.4, markersize = 1, label = 'sample')

plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1),
          ncol=3, fancybox=True, shadow=True)
plt.show()

###############################################################################
