# -*- coding: utf-8 -*-
"""
Created on Sun May 17 22:22:46 2015

@author: dima
"""

m = np.matrix

T = m([[1,2],[3,4]])
e = m([1,3]).T
E = m(np.diag(np.array(e).ravel()))

p = m([1,2])
P = m(np.diag( np.array(p).ravel()))

#######################
Alpha = (E * T.T)
Beta =  (E * T)
# Beta = np.array(e) * np.array(T)
Alpha = 
np.array(e) * np.array(T)

one_.T * (Alpha**4) * P * e

one_.T * P * (Beta**4 ) * e

one_.T * P * (Beta**4  ) * E * one_

one_.T * Alpha **2 * P * Beta **2 * e

one_.T * Alpha  * P * Beta **3 * e

one_.T * Alpha **3 * P * Beta * e
