# -*- coding: utf-8 -*-
"""
Created on Tue Nov  4 10:58:02 2014

@author: Dmitri
"""

import numpy as np

import scipy.misc as sm
import scipy.linalg as la
import warnings

def logsumexp10(x, *args, **kwargs):
    return sm.logsumexp( x * np.log(10), *args, **kwargs) / np.log(10)
calcMarginal = logsumexp10
#class SmallNegativeValuesWarning(Warning, ValueError):
#    pass
class ErrorAllNaNs(Exception):
    def __init__(self):
        self.value = 'All entries in the vector are NaN!'
    def __str__(self):
        return self.value

class hmm_cont:
    
    def __init__(self, population_obj, emission_matrix, dist_genetic = None):
        self.E = emission_matrix;
        self.hidstates = population_obj;
        self.M = self.E.shape[1]
        self.Np = self.E.shape[0]
        
        assert self.hidstates.Np == self.Np,  "population size mismatch: internal %u, external %u" % (self.Np, self.hidstates.Np)
        
        if (isinstance(dist_genetic, (np.ndarray, np.generic) )) :
            self.t = dist_genetic
            assert(len(self.t) == self.M )
        else:
            self.t = np.array(dist_genetic) # self.T_E_productrray([]);

    def calcT(self, linkage_loosening=1.0):
        self.linkage_loosening = linkage_loosening

        warnings.filterwarnings('once', '.*small.*', UserWarning)

        assert len(self.t) > 0, 'the node (recombination) distances "Transition" have not been set!'
        assert (isinstance(self.linkage_loosening, float) or \
                isinstance(self.linkage_loosening, int)), \
            'linkage_loosening factor must be a scalar'

        """ calculate Transition matrices """
        self.Transition = np.zeros((self.M - 1, self.hidstates.Np, self.hidstates.Np));
        delta_t = np.diff(self.t, 1) * self.linkage_loosening;
        """ calculate """
        for ii in range(0, self.M - 1):
            if (delta_t[ii] != 0) and not np.isinf(delta_t[ii]):
                self.Transition[ii, :, :] = la.expm(self.hidstates.Q * delta_t[ii]);
            elif (delta_t[ii] == 0):
                self.Transition[ii, :, :] = np.identity(self.hidstates.Np);
            elif np.isinf(delta_t[ii]):
                self.Transition[ii, :, :] = np.tile(self.hidstates.Pstat, self.hidstates.Np);

        """sanity check : negative values"""
        if any(self.Transition.ravel() < 0):
            assert all(abs(self.Transition.ravel()[
                self.Transition.ravel() < 0.0]) < 1e-20), 'Transition matrix has negative values!'
            warnings.warn('\n'.join(['Transition matrix has very small negative values (possibly numeric error)!',
                                     'Replacing negative numbers by zeros']))

            self.Transition.ravel()[self.Transition.ravel < 0] = 0;

    def crossMatr(self, linkage_loosening=1.0):
        assert isinstance(self.E, (np.ndarray, np.generic)), "provide emission matrices first!"

        if not hasattr(self, 'Transition') or \
                not (isinstance(self.Transition, (np.ndarray, np.generic))):
            """define transition matrix T first!"""
            self.calcT(linkage_loosening)

        resh_E = np.transpose(self.E[:, 0:self.M - 1][:, np.newaxis], (2, 0, 1))
        resh_E.shape

        self.TxE_B  = np.multiply(resh_E, np.transpose(self.Transition, (0, 2, 1)) );
        self.TxE_A = np.multiply(resh_E, self.Transition);

        if not (self.TxE_A == 0)[:, :, 0].all() and \
        (self.TxE_A == 0)[:, :, :-1].all():
            print( 'non-zero entries: ')
            print( np.where(self.TxE_A[:, :, 0] != 0) )
            tmp = self.TxE_A[:, :, 0]
            print( tmp[np.where(tmp != 0)] )
            # print( np.where((self.TxE_A == 0)[:, :, 0].all) )
            
            import matplotlib.pyplot as plt
            plt.figure(num = 1)
            plt.pcolor(self.TxE_A[:, :, 0])
           
            assert len(self.t)>0, 'the node (recombination) distances "Transition" have not been set!'
            assert (isinstance(self.linkageLoosening, float) or isinstance(self.linkageLoosening, int)), 'linkageLoosening factor must be a scalar'
            """ calculate Transition matrices """
            self.Transition = np.zeros( (self.M-1, self.hidstates.Np, self.hidstates.Np) );
            delta_t = np.diff(self.t, 1);
            """ calculate """
            for ii in range(0, self.M-1 ):
                if (delta_t[ii] !=0) and not np.isinf(delta_t[ii]) :
                    self.Transition[ii, :,:] = la.expm(self.hidstates.Q * delta_t[ii]);
                elif  (delta_t[ii]==0) :
                    self.Transition[ii, :,:] = np.identity(self.hidstates.Np);
                elif np.isinf(delta_t[ii]):
                    self.Transition[ii, :,: ] = np.tile(self.hidstates.Pstat, self.hidstates.Np);

            """sanity check : negative values"""
            if any(self.Transition.ravel()<0):
                assert all(abs(self.Transition.ravel()[self.Transition.ravel<0] ) < 1e-25) , 'Transition matrix has negative values!'
                warnings.warn( '\n'.join(['Transition matrix has very small negative values (possibly numeric error)!',
               'Replacing negative numbers by zeros']))

                self.Transition.ravel()[self.Transition.ravel<0] = 0;

    def cumMatr(self, linkage_loosening=1.0):
        assert isinstance(self.Np, int), "please define the number of states (`Np`) first"

        if not hasattr(self, 'TxE_A') or \
        not isinstance(self.TxE_A, (np.ndarray, np.generic)):
            """calculate the Transition-Emission product matrix first"""
            self.crossMatr();
        
        self.logAlpha = -float('inf') * np.ones( ( self.Np, self.M) );
        self.logBeta  = -float('inf') * np.ones( ( self.Np, self.M) );
        
        """forward"""
        aCumul = self.E[ :, -1]; # len(aCumul) == self.Np   """number of states"""
        self.logAlpha[:, self.M-1] = np.log10(aCumul);
        
        scalePrev = 0;
        scaleA = np.zeros(self.M);
        
        for m in range(self.M - 2, -1, -1):
            order_of_magnitude = 10 ** (scaleA[m + 1] - scalePrev) 
            aCumul = np.dot(self.TxE_A[m, :, :], aCumul ) * order_of_magnitude 
            self.logAlpha[:, m] = np.log10(aCumul) - scaleA[m + 1];
            scalePrev = scaleA[m + 1];
            scaleA[m] = - max( self.logAlpha[:, m] );
        
        """backward"""
        self.logBeta = -float('inf') * np.ones(( self.Np, self.M));
        # self.Transition[ 0].T * self.E[:,0] - self.TxE_A[0]
        bCumul =  np.ones(self.Np) # self.E[:, 0] # 
        self.logBeta[:, 0] =  np.log10(bCumul);

        scalePrev = 0;
        scaleB = np.zeros( self.M )
        
        for m in range(1, self.M):
            order_of_magnitude = 10 ** (scaleB[m - 1] - scalePrev)
            bCumul = np.dot(bCumul * order_of_magnitude, self.TxE_B[m - 1, :, :]) 
            self.logBeta[:, m] = np.log10(bCumul) - scaleB[m - 1]
            scalePrev = scaleB[m - 1];
            scaleB[m] = - max(self.logBeta[:, m]);
        
        return

    def _runFB_(self):
            if not hasattr(self, 'T_E_product') or not hasattr(self, 'logAlpha' ) or not hasattr(self, 'logBeta' ) :
                self.crossMatr();
                self.cumMatr();
            
            self.xk_P_flat = self.logAlpha + self.logBeta 
            # log10(self.Pz);
            
            # self.x_P_flat = calcMarginal(self.xkPout, axis = 1);
            if np.isnan(self.xk_P_flat ).any():
                warnings.warn('some entries in the probability matrix are NaN!')
            
            if np.isnan(self.xk_P_flat).all():
                raise ErrorAllNaNs
                # exception('all entries in the probability vector are NaN!')
    
    def getLikelihoodOfAModel(self, model_P_z):
        # assert hasattr(self, 'Pz') ,'no distribution over the hidden states has been submitted!'
        assert (isinstance(model_P_z, float)) or (len(model_P_z) == self.hidstates.Np) , \
        'a wrong distribution over the hidden states z_m is submitted!'
        
        if not hasattr(self, 'xk_P_flat'):
            self._runFB_()
        
        xk_P_out = self.xk_P_flat + np.log10(model_P_z)[np.newaxis].T
        x_P_out = calcMarginal(xk_P_out, axis = 0)
        return (x_P_out, xk_P_out)
 
