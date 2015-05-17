# -*- coding: utf-8 -*-
"""
Created on Tue Jan 20 20:37:51 2015

@author: dima
"""

import sys

sys.path.append('../../snpdb/')
from snp_tables import *
from genetic_map import *
import weakref

from unmix_repeats import *
from population import *
from hmm_cont import *
from scipy.stats import binom
import numpy as np


def close(series, ele):
    series = dilate(series, ele)
    series = erode(series, ele)
    return series


def dilate(series, ele):
    new_series = np.copy(series)
    for n in range(len(series)):
        new_series[n] = max(series[max(0, n - ele):min(len(series), n + ele + 1)])
    return new_series


def erode(series, ele):
    new_series = np.copy(series)
    for n in range(len(series)):
        new_series[n] = min(series[max(0, n - ele):min(len(series), n + ele + 1)])
    return new_series


# f = q/r;
def binomial_emission(pop, q, r):
    return binom.pmf(q, r, pop.f_vect)


def binomial_uniform_emission(pop, q, r, membership=None):
    p_bin = binom.pmf(q, r, pop.f_vect)
    if membership is None:
        return p_bin
    else:
        assert (membership.shape == q.shape), "dimension mismatch"

        out_emission = membership * p_bin + (1 - membership) * (np.sum(p_bin, axis=0) / pop.Np)
        assert (out_emission >= 0).all(), "negative values in the emission matrix!"
        # assert (out_emission <= 1).all()
        # assert ( abs(np.sum(out_emission , axis = 0)-1) < 1e-20).all()
        return out_emission


###############################################################################
class ChromosomeLength:
    def __init__(self, path):
        self.L_dict = {}
        self.L = []
        with open(path) as fi:
            for line in fi:
                cols = line.split('\t')
                chr_length = int(cols[1])
                self.L_dict[cols[0].split(' ')[0]] = chr_length
                self.L.append(chr_length)

    def __getitem__(self, chromosome):
        return self.L[chromosome - 1]

    def __repr__(self):
        st = 'chromosome lengths:\n'
        for nn, ll in enumerate(self.L):
            st += '%u\t%u\n' % (nn + 1, ll)
        return st

    __str__ = __repr__


###############################################################################
def check_str_array_field(data, field):
    for s in data.__array_interface__['descr']:
        if s[0] == field:
            return True
    return False


###############################################################################
class SnpReadsDbChr():
    def __init__(self, parent, chromosome_name, \
                 condition_):
        self.membership_set = False
        self.chromosome_name = chromosome_name
        self.parent = weakref.ref(parent)  # <= garbage-collector safe parent tracker
        print(condition_)
        self._read_input_table_(condition_)
        self._set_length_()

    def _set_length_(self):
        if hasattr(self.parent(), 'chromosome_lengths'):
            self.length = self.parent().chromosome_lengths[self.chromosome_name]
        else:
            self.length = max(self.data['x'])

    def prepare(self):
        if not self.membership_set:
            print("membership not set!")
            self.unmix()
        self.emission_array = binomial_uniform_emission(self.parent().pop, \
                                                        self.data['q'], self.data['r'], self.data['membership'])

        "genetic distance for given data points : `gen_x`"
        self.gen_x = self.parent().gen_map(self.data['x'], self.chromosome_name)
        self.HMM = hmm_cont(self.parent().pop, self.emission_array, self.gen_x)

    def _read_input_table_(self, condition_=''):
        self.membership_set = False
        if self.column_exists('membership'):
            res = self.parent().db_table.selectPositionArray(self.chromosome_name, \
                                                             pos='*', columns='pos, altCount, totCount, membership',
                                                             condition_=condition_)
        else:
            res = self.parent().db_table.selectPositionArray(self.chromosome_name, \
                                                             columns='pos, altCount, totCount', condition_=condition_)
            # self.x = res[:,0]
        #        self.q = res[:,1]
        #        self.r = res[:,2]
        #        self.f = self.q/self.r
        print('shape: ', res.shape)
        self.data = np.empty(res.shape[0],
                             dtype=[('x', 'i4'), ('q', 'i4'), ('r', 'i4'), ('dx', 'i4'), ('membership', 'f4')])
        self.data['x'] = res[:, 0]
        self.data['q'] = res[:, 1]
        self.data['r'] = res[:, 2]

        if res.shape[1] > 3 and not ((res[:, 3] == 0).all()):
            self.data['membership'] = res[:, 3]
            self.membership_set = True
        assert (np.diff(self.data['x']) > 0).all()

    def get_columns(self):
        res = self.parent().db_table.run_query('PRAGMA table_info(%s)' % self._get_table_name_())
        return res

    def column_exists(self, column_name):
        columns = self.get_columns()
        for rr in columns:
            if (rr[1] == column_name):
                return True
        return False

    def _get_table_name_(self):
        return self.parent().db_table.get_name(self.chromosome_name)

    def add_membership(self):
        self.unmix()
        self.parent().db_table.add_membership(self.chromosome_name, \
                                              self.data['x'], self.data['membership'])
        print('membership data has been added to the table, chromosome %u' % \
              self.chromosome_name, file=sys.stderr)


    def run_query(query):
        res = self.parent().db_table.run_query(query)
        return res

    def run(self, plot=False):
        self.p_x_sel = self.HMM._runFB_(prior='Psel')
        self.p_x_stat = self.HMM._runFB_(prior='Pstat')
        if plot:
            self.plot_both()

    def plot_both(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        fig.suptitle('probability')
        ax.plot(self['x'] * 1e-6, self.p_x_sel, 'r*-', label='selection')
        ax.plot(self['x'] * 1e-6, self.p_x_stat, 'g-', label='no selection')
        plt.show()
        return fig

    def unmix(self):
        print('unmixing', file=sys.stderr)
        dx_raw = np.diff(self.data['x'])
        np.concatenate(([np.Inf], dx_raw[:-1]))
        self.data['dx'] = np.min([np.concatenate(([np.Inf], dx_raw)),
                                  np.concatenate((dx_raw, [np.Inf] ))], axis=0)
        log10_dx = np.log10(self['dx'])
        self.data['membership'], self.mixt_obj, iTrue, out_mode_number, _ = unmix(log10_dx)
        self.membership_set = True

    def __getitem__(self, name):
        return self.data[name]

    def find_dense_snps_grid(self, scale=1, shift=0):
        "characteristic length scale of the Poisson distribution"
        mu = scale * int(self.length / self.data.shape[0])
        x_grid = np.arange(0 - shift * mu, self.length, mu)
        snps_on_grid = np.empty(x_grid.shape, dtype=int)

        for ii in range(0, len(x_grid)):
            snps_on_grid[ii] = sum(abs(self['x'] - x_grid[ii]) < mu)

        return (x_grid, snps_on_grid)

    def find_dense_snps(self, scale=1, shift=0):
        "characteristic length scale of the Poisson distribution"
        mu = scale * int(self.length / self.data.shape[0])
        x_grid = np.copy(self['x'])
        snps_on_grid = np.empty(x_grid.shape, dtype=int)

        for ii in range(0, len(x_grid)):
            snps_on_grid[ii] = sum(abs(self['x'] - x_grid[ii]) < mu)

        return (snps_on_grid)

    def threshold_dense_snps(self, scale=1, shift=0, threshold=3, inv_subscale=10, grid=True):
        x = [None] * inv_subscale;
        y = [None] * inv_subscale
        if grid:
            for nn in range(inv_subscale):
                x[nn], y[nn] = self.find_dense_snps_grid(scale, nn / inv_subscale)
            x = np.fliplr(np.array(x).T).ravel()
            y = np.fliplr(np.array(y).T).ravel()
        else:
            y = np.array(self.find_dense_snps(scale, 0))
            x = np.copy(self['x'])
            inv_subscale = 1

        self.crowded_region = close(y > threshold * scale, 1 * inv_subscale)
        return (self.crowded_region, y, x)

    def plot_snp_density(self, scale=1, shift=0, threshold=3, inv_subscale=10, grid=True):
        self.crowded_region, y, x = self.threshold_dense_snps(scale=1, shift=0, threshold=3, inv_subscale=10, grid=True)

        import matplotlib.pyplot as plt

        y_max = np.percentile(y, 99);
        fig = plt.figure()
        fig.suptitle('snp density')
        plt.plot(x * 1e-6, y, 'r.-', label = 'density')
        plt.plot([0, self.length * 1e-6], np.array([1, 1]) * threshold * scale, 'g-', label = 'threshold')
        plt.plot(x * 1e-6, y_max * self.crowded_region, 'b-', label = 'crowded region label')
        plt.legend()        
        plt.ylim((0, np.ceil(1.1*y_max))) 
        plt.show()

        return self.crowded_region, x


###############################################################################
class SnpReadsData():
    def __init__(self, file_path, experiment_name, pop_size, emmission_type=None, \
                 lengths_path='reference/TAIR10-chr-counts.dat'):
        self.pop = population(pop_size);
        self.db_table = GenomeVariantDbReader(file_path, experiment_name, subtable='mtCounts')
        self.snp_reads_on_chrs = {}
        self.gen_map = GeneticMap()
        self._set_chromosome_lengths_(lengths_path)

    def _set_chromosome_lengths_(self, path=''):
        if len(path) == 0:
            return
        else:
            self.chromosome_lengths = ChromosomeLength(path)

    def read_chromosome(self, chromosome_name, \
                        condition_='(totCount >=10) AND (altCount >= 7) AND (CAST(altCount AS FLOAT) / CAST(totCount AS FLOAT) < 0.8)'):

        self.snp_reads_on_chrs[chromosome_name] = SnpReadsDbChr(self, chromosome_name, condition_)
        self[chromosome_name].prepare()

    def check_data_tables(self):
        tables = self.db_table.showTables()
        print(tables)

        # def run(self):

    #        for self[]

    def __getitem__(self, chromosome_name):
        """
        call `s = SnpReadsData(...); s[chromosome_name]`
        to access child `SnpReadsDbChr` objects containing data for `chromosome_name`
        """
        if chromosome_name in self.snp_reads_on_chrs:
            return self.snp_reads_on_chrs[chromosome_name]
        else:
            try:
                self.read_chromosome(chromosome_name)
                return self.snp_reads_on_chrs[chromosome_name]
            except:
                chromosome_name = chromosome_name if isinstance(chromosome_name, str) else '%u' % chromosome_name
                raise KeyError('chromosome %s is not found in the data.' % chromosome_name)
