# -*- coding: utf-8 -*-
"""
Created on Wed Nov 12 19:40:08 2014

@author: D S Lituiev
"""

import numpy as np
from sklearn import mixture
from hmm_cont import *
###############################################################################
def min_distances(x):
    dx_raw = np.diff(x[:])
    np.concatenate(([np.Inf], dx_raw[:-1]))
    dx = np.min([np.concatenate(([np.Inf], dx_raw)),
                 np.concatenate((dx_raw, [np.Inf] ))], axis=0)
    log10_dx = np.log10(dx)
    return dx,log10_dx
###############################################################################
def unmix(log_x, max_mode_number=3):
    """
    a wrapper for a Gaussian mixture model;
    returns BIC optimal number of mixtures
    """
    mixtBICs = np.zeros(max_mode_number)
    mixt = [None] * max_mode_number

    for n_modes in range(0, max_mode_number):
        mixt[n_modes] = mixture.GMM(n_components=n_modes + 1)
        mixt[n_modes].fit(log_x)
        mixtBICs[n_modes] = mixt[n_modes].bic(log_x)

    " find the BIC-optimal number of the mixture components ( >=2 ) "
    out_mode_number = np.argmin(mixtBICs)
    " use the BIC-optimal GMM parameters "
    mixtObj = mixt[out_mode_number]

    del (mixt)

    iTrue = np.argmax(mixtObj.means_);
    "  P[x| component_true ] "
    piTrue = mixtObj.weights_[iTrue];

    membership_funct = lambda x: mixtObj.predict_proba(x)[:, iTrue]
    Membership = membership_funct(log_x)
    return (Membership, mixtObj, iTrue, out_mode_number+1, membership_funct)


###############################################################################
class hmm_cont_unmixer(hmm_cont):
    def unmix_repeats():
        firstParam = 'dx'
        max_mode_number = 2

        x = float(getattr(obj, firstParam));
        """ combine the log-scaled feature vectors into one array"""
        log_x = np.log10(x(~np.isnan(x)));

        (Membership, mixtObj, iTrue) = unmix(log_x, max_mode_number)