# -*- coding: utf-8 -*-
"""
reader for chromosome genetic map

Created on Mon Jan 19 13:59:50 2015

@author: dima
"""

import sys
import numpy as np
from scipy.interpolate import PchipInterpolator
import matplotlib.pyplot as plt


class GeneticMap:
    "the genetic distances are stored in Morgans (not cM, i.e. Morgans = cM/100)"

    def __init__(self, file_path='reference/Athaliana_Phys_Rec_Corresp_Singer.csv'):
        self._count_lines_(file_path, tot_chromosome_num=5)
        self._initialize_arrays_()
        self._read_table_(file_path)
        "assert that the map points are non-decreasing"
        for cc in range(len(self._phys_pos_)):
            assert (np.diff(self._phys_pos_[cc]) > 0).all(), "decreasing values found!"

        print('initialization of the genetic map is completed', file=sys.stderr)

    def _count_lines_(self, file_path, tot_chromosome_num=5):
        # num_lines = int(check_output(['wc', '-l', file_path]).decode().split(' ')[0]) - 1
        self._num_lines_ = []

        from subprocess import check_output, Popen, PIPE

        for cc in range(tot_chromosome_num):
            ps = Popen(['cut', '-d', ';', '-f', '2', file_path], stdout=PIPE)
            self._num_lines_.append( \
                int(check_output(('grep', '-c', '%u' % (cc + 1) ), stdin=ps.stdout)))
            ps.wait()
        # print(num_lines)
        return self._num_lines_

    def _initialize_arrays_(self):
        "the first entry must be zero in each array: allocate a spare entry for it"
        self._phys_pos_ = []
        self._gen_pos_ = []
        for nn in self._num_lines_:
            # print(nn)
            self._phys_pos_.append(np.zeros((1 + nn,), dtype=int))
            self._gen_pos_.append(np.zeros((1 + nn,)))

    def _read_table_(self, file_path):
        " "
        ii = 1
        with open(file_path) as fi:
            header = next(fi)
            chromosome = -1
            for line in fi:
                lspl = line.split(';')
                if not (chromosome == int(lspl[1])):
                    chromosome = int(lspl[1])
                    ii = 1
                self._phys_pos_[chromosome - 1][ii] = (int(lspl[3]) + int(lspl[4])) // 2
                self._gen_pos_[chromosome - 1][ii] = 0.01 * float(lspl[5])
                ii += 1
        return

    def visualise_map(self, cc=1):
        plt.plot(self._phys_pos_[cc - 1], self._gen_pos_[cc - 1], 'o-')
        plt.show()

    def interpolate(self, phys_x, chromosome):
        " note that the interpolation scheme does not guarantee non-decreasing output"

        NUM_LAST_POINTS = 5
        assert (np.diff(phys_x) > 0).all(), "decreasing input values found!"
        "using  Piecewise Cubic Hermite Interpolating Polynomial"
        ipf = PchipInterpolator(self._phys_pos_[chromosome - 1], \
                                self._gen_pos_[chromosome - 1], extrapolate=False)
        gen_x = ipf(phys_x)

        nan_inds, = np.where(gen_x != gen_x)
        # print("extrapolating outer points:", file=sys.stderr)
        # print(nan_inds, file=sys.stderr)

        "last segment"
        inds = self._last_inds_(chromosome, phys_x)

        g_last = self._fit_last_points_(chromosome, 2)
        if inds.any():
            w_last = self._weight_last_(chromosome)(phys_x[inds])

            gen_x[inds] = g_last(phys_x[inds]) * (1 - w_last) + \
                          gen_x[inds] * w_last

        "hanging segment"
        inds = self._hanging_inds_(chromosome, phys_x)
        if inds.any():
            w_h = self._weight_hanging_(chromosome)(phys_x[inds])
            g_h = self._fit_last_points_(chromosome, num_last_points=NUM_LAST_POINTS)

            gen_x[inds] = g_h(phys_x[inds]) * (1 - w_h) + \
                          g_last(phys_x[inds]) * w_h
        "==="

        # for nn in range(0, NUM_LAST_POINTS):
        #            if (gen_x[nan_inds[0]-nn] >= gen_x[nan_inds[0]]):
        #                print("resetting positin # %u" % (nan_inds[0] - nn), file=sys.stderr)

        "assert non-decreasing"
        if not (np.diff(gen_x) > 0).all():
            print("decreasing output values found!", file=sys.stderr)
            nz = [i for i, e in enumerate(np.diff(gen_x) > 0) if not e]
            print(nz, file=sys.stderr)
            self.visualise_map(chromosome)
            plt.plot(phys_x, gen_x, 'gx')
            plt.show()
        assert (np.diff(gen_x) > 0).all(), "decreasing output values found!"

        return gen_x

    def _last_x_(self, chromosome):
        return self._phys_pos_[chromosome - 1][-1]

    def _last_inds_(self, chromosome, phys_x):
        li = phys_x > self._phys_pos_[chromosome - 1][-2]
        li = li & ( phys_x < self._phys_pos_[chromosome - 1][-1])
        return li

    def _hanging_inds_(self, chromosome, phys_x):
        hi = phys_x > self._phys_pos_[chromosome - 1][-1]
        return hi

    def _scale_(self, chromosome):
        return self._phys_pos_[chromosome - 1][-1] - \
               self._phys_pos_[chromosome - 1][-2]

    def _weight_(self, scale, x):
        return np.exp2(- (x / scale / 2) ** 2)

    def _weight_last_(self, chromosome):
        return lambda phys_x: \
            ( (self._last_x_(chromosome) - phys_x) / self._scale_(chromosome))

    def _weight_hanging_(self, chromosome):
        return lambda x: self._weight_(self._scale_(chromosome), x - self._last_x_(chromosome))

    def _fit_last_points_(self, chromosome, num_last_points=5):
        fit = np.polyfit(self._phys_pos_[chromosome - 1][-num_last_points:], \
                         self._gen_pos_[chromosome - 1][-num_last_points:], 1)
        line = np.poly1d(fit)
        return line

    __call__ = interpolate
        