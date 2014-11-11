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
        self.pop = population_obj;
        self.M = self.E.shape[1]
        self.Np = self.E.shape[0]
        
        assert self.pop.Np == self.Np,  "population size mismatch: internal %u, external %u" % (self.Np, self.pop.Np)
        
        if (isinstance(dist_genetic, (np.ndarray, np.generic) )) :
            self.t = dist_genetic
            assert(len(self.t) == self.M )
        else:
            self.t = np.T_E_productrray([]);
    
    def calcT(self, linkageLoosening = 1.0):
            self.linkageLoosening = linkageLoosening
            
            warnings.filterwarnings('once', '.*small.*', UserWarning)            
            
            assert len(self.t)>0, 'the node (recombination) distances "Transition" have not been set!'
            assert (isinstance(self.linkageLoosening, float) or isinstance(self.linkageLoosening, int)), 'linkageLoosening factor must be a scalar'
            """ calculate Transition matrices """
            self.Transition = np.zeros( (self.M-1, self.pop.Np, self.pop.Np) );
            delta_t = np.diff(self.t, 1);
            """ calculate """
            for ii in range(0, self.M-1 ):
                if (delta_t[ii] !=0) and not np.isinf(delta_t[ii]) :
                    self.Transition[ii, :,:] = la.expm(self.pop.Q * delta_t[ii]);
                elif  (delta_t[ii]==0) :
                    self.Transition[ii, :,:] = np.identity(self.pop.Np);
                elif np.isinf(delta_t[ii]):
                    self.Transition[ii, :,: ] = np.tile(self.pop.Pstat, self.pop.Np);

            """sanity check : negative values"""
            if any(self.Transition.ravel()<0):
                assert all(abs(self.Transition.ravel()[self.Transition.ravel<0] ) < 1e-25) , 'Transition matrix has negative values!'
                warnings.warn( '\n'.join(['Transition matrix has very small negative values (possibly numeric error)!',
               'Replacing negative numbers by zeros']))

                self.Transition.ravel()[self.Transition.ravel<0] = 0;
    
    def crossMatr(self):
            assert isinstance(self.E, (np.ndarray, np.generic) ),  "provide emission matrices first!"
          
            if not (isinstance(self.Transition, (np.ndarray, np.generic) )):
                """define transition matrix T first!"""
                self.calcT() 
                
            resh_E = np.transpose(self.E[:, 0:self.M - 1][:,np.newaxis], (2,0,1))
            resh_E.shape
    
            self.T_E_product = np.transpose( np.multiply( resh_E, self.Transition ), (0,2,1));
            
            assert (self.T_E_product == 0)[:,:,0].all()
            
    def cumMatr(self):
        assert isinstance(self.Np, int), "please difine the number of states (`Np`) first"
        
        if not hasattr(self, 'T_E_product') or not isinstance(self.T_E_product, (np.ndarray, np.generic)):
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
            aCumul = np.dot( aCumul, self.T_E_product[m, :, :]) * 10**(scaleA[m+1] - scalePrev) 
            self.logAlpha[:, m]= np.log10(aCumul) - scaleA[m+1];
            scalePrev = scaleA[m + 1];
            scaleA[m] = - max( self.logAlpha[:, m] );
        
        """backward"""
        bCumul = np.ones(self.Np);
        self.logBeta[ :, 0 ] = 0;
        
        scalePrev = 0;
        scaleB = np.zeros( self.M )
        
        for m in range(1, self.M):
            bCumul = np.dot( bCumul,  self.T_E_product[m-1, :,:]) *  10**(scaleB[m-1] - scalePrev) 
            self.logBeta[:, m] = np.log10(bCumul) - scaleB[m-1];
            scalePrev = scaleB[m-1];
            scaleB[m] = - max(self.logBeta[:, m] );
        
        return

    def runFBinternal(self):
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
        assert (isinstance(model_P_z, float)) or (len(model_P_z) == self.pop.Np) , \
        'a wrong distribution over the hidden states z_m is submitted!'
        
        xk_P_out = self.xk_P_flat + np.log10(model_P_z)[np.newaxis].T
        x_P_out = calcMarginal(xk_P_out, axis = 0)
        return (x_P_out, xk_P_out)
 