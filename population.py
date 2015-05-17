import numpy as np
from scipy.special import gammaln
from scipy.misc import comb


class OutOfRangeError(Exception):
    pass


class population:
    """a class to store population properties, incl. infinitisemal transition
    martix Q for the MapBySeqHmm project"""

    def __init__(self, N_in, crossType="backcrossed"):

        if not isinstance(N_in, int) or not (N_in > 1):
            raise OutOfRangeError("the `N_in` should be more than 1")

        self.N = N_in;
        self.Np = N_in + 1;
        self.k_vect = np.arange(0, N_in + 1)[np.newaxis].T;

        crossDict = {'homozygous': False, 'selfed': False, 'self': False, 's': False, \
                     'heterozygous': True, 'backcrossed': True, 'back': True, 'b': True}

        if ( crossType in crossDict ) and (not crossDict[crossType] ):
            self.f_vect = self.k_vect / self.N;
        else:
            """ heterozygous population from a back-cross """
            self.f_vect = self.k_vect / (2 * self.N);

        self.calcQmatrix();
        self.stationaryDistribution();
        self.Pflat = [1 / self.Np] * self.Np;

        self.Psel = [0 for ii in range(self.Np)]
        self.Psel[self.Np - 1] = 1

    def stationaryDistribution(self):
        self.Pstat = [0] * self.Np;

        if self.N <= 20:
            for ii in range(0, self.Np):
                self.Pstat[ii] = 2 ** (-self.N) * comb(self.N, ii);
        else:
            """ for larger numbers """
            for ii in range(0, self.Np):
                self.Pstat[ii] = 2 ** (-self.N) * round(
                    np.exp(gammaln(self.N + 1) - gammaln(ii + 1) - gammaln(self.N - ii + 1)));
        self.Pstat /= sum(self.Pstat)

    def calcQmatrix(self):
        """ infinitisemal transition matrix
        Calculates infinitesimal transition matrix Q
        for continuous Ehrenfest process

        ===== Input: =====
        N - number of plants/molecules
        """

        """ short k-vector"""
        k_vect_sh = np.arange(0, self.N)
        """ initialize """
        Q_dimensions = (self.Np, self.Np)
        self.Q = np.zeros(Q_dimensions);
        """ over diag: increment """
        # subs = np.vstack( ( np.arange(0, self.N-1), np.arange(1, self.N) ))
        inds = np.arange(1, self.Np ** 2, self.N + 2)
        self.Q.ravel()[inds] = self.N - k_vect_sh;
        """ under diag: decrement """
        inds = np.arange(self.Np, self.Np ** 2, self.N + 2)
        self.Q.ravel()[inds] = k_vect_sh + 1;
        """ staying same """
        inds = np.arange(0, self.Np ** 2, self.N + 2)
        self.Q.ravel()[inds] = - np.sum(self.Q, axis=1);
        # Q is done!