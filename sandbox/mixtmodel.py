# -*- coding: utf-8 -*-
"""
Created on Wed Nov 12 17:04:02 2014

@author: Dmitri
"""

import numpy as np
from sklearn import mixture


obs = np.concatenate((np.random.randn(100, 1),
                       2 + np.random.randn(300, 1)))
          
g = mixture.GMM(n_components=2)             
g.fit(obs)
np.round(g.weights_, 2)
np.round(g.means_, 2)

g.predict_proba(obs)

pr = g.predict_proba(obs)
np.sum(pr , axis=1)


g.bic(obs)

import matplotlib.pyplot as plt

plt.plot(pr[:,1], pr[:,2], 'ro')

plt.show()