# -*- coding: utf-8 -*-
"""
Created on Mon Jan 19 15:36:51 2015

@author: dima
"""

import numpy as np
import numpy.random as nr
import matplotlib.pyplot as plt
from genetic_map import *
# runfile('/home/dima/data/repos/MapBySeqHmmPy/src/read_genetic_map.py', wdir=r'/home/dima/data/repos/MapBySeqHmmPy/src')
file_path = 'reference/Athaliana_Phys_Rec_Corresp_Singer.csv'
gm = GeneticMap(file_path)

gm.visualise_map()
x = np.linspace(0, max(gm._phys_pos_[0]), 1e2)

plt.plot(x, gm(x, 1), 'rx-')
plt.show()

gm([34534, 123124], 2)

chrs = 4
L = gm._last_x_(1) + 4

M = 500
x = np.array(nr.random_integers(1, L, M))
x.sort(axis=0)

plt.figure()
gm.visualise_map(chrs)
plt.plot(x, gm(x, chrs), 'rx-')
plt.show()