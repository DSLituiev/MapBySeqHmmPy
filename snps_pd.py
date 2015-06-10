# -*- coding: utf-8 -*-
"""
Created on Sat Jun  6 17:06:01 2015

@author: dima
"""

import sys
sys.path.append('../../snpdb/')
from snp_tables import *
from population import *
from unmix_repeats import *
from hmm_cont import *
from chromosome_length import ChromosomeLength
from genetic_map import *
from snp_density import *

from binomial_emission import *
import pandas as pd

import os.path

class snps_pd(pd.DataFrame):
    membership_set = False
    @property
    def _constructor(self):
        return snps_pd
    
    def __init__(self, file_path, experiment_name = None, pop_size = None, \
            conditions = '', lengths_path='reference/TAIR10-chr-counts.dat'):
        if type(file_path) == snps_pd:
            self = file_path
            return
        
        if not type(file_path) == str or \
            not os.path.isfile(file_path) and \
            not type(experiment_name) == str and \
            not (type(pop_size) == int):
            super(snps_pd, self).__init__(data=file_path, index=experiment_name, columns=pop_size, dtype=None,
                 copy=False)
            return
        
        columns ='chr, pos, altCount, totCount'
        db_table = GenomeVariantDbReader(file_path, experiment_name, subtable='mtCounts')
        db_table.countRows()
        super(snps_pd, self).__init__( pd.read_sql_query( 'SELECT %s ' % columns + \
            db_table._from_where_(None, conditions),  
            db_table.connection) )
        self.db_table = db_table
        self.pop = population(pop_size);

        self.snp_reads_on_chrs = {}
        self.gen_map = GeneticMap()
        self._set_chromosome_lengths_(lengths_path)
        
        self.rename(columns={'pos': 'x', 'totCount':'r', 'altCount': 'q'}, inplace=True)
        self['Mb'] = self['x']*1e-6
        self.set_index(['chr'], inplace=True)
        return
#        print( '==========')
#        for chromosome, chr_data in self.groupby(level=0):
#            print( '%s', '%u entries' % (chromosome, chr_data.shape), sep = '\t', file=sys.stderr)
#        print( '==========')
        
    def _set_chromosome_lengths_(self, path=''):
        if len(path) == 0:
            return
        else:
            self.chromosome_lengths = ChromosomeLength(path)

    def unmix(self):        
        for chromosome, chr_data in self.groupby(level=0):
            print('unmixing %s' % chromosome, file=sys.stderr)
            dx_raw = np.diff(chr_data['x'])
            self.loc[chromosome, 'dx'] = np.min([np.concatenate(([np.Inf], dx_raw)),
                                  np.concatenate((dx_raw, [np.Inf] ))], axis=0)
            self.loc[chromosome,'log10_dx'] = np.log10(self.loc[chromosome, 'dx'])
            self.loc[chromosome,'membership'], self.mixt_obj, iTrue, out_mode_number, _ = \
                unmix( self.loc[chromosome,'log10_dx'] )
        
        self.membership_set = True


    def prepare(self):
        if not self.membership_set:
            print("membership not set!", file=sys.stderr)
            self.unmix()
        self.HMM = {} 
        for chromosome, chr_data in self.groupby(level=0):
            "genetic distance for given data points : `gen_x`"
            self.loc[chromosome, 'gen_x'] = self.gen_map(chr_data['x'], chromosome)
            if not np.isnan(self.loc[chromosome, 'gen_x'].iloc[0]):
                print('preparing emission and transition for %s' % chromosome, file=sys.stderr)
                "skip contigs with no genetic map (which is set to `NaN`)"
                emission_array = binomial_uniform_emission(self.pop, \
                    chr_data['q'], chr_data['r'], chr_data['membership'])
                self.HMM[chromosome] = hmm_cont(self.pop, emission_array, self.loc[chromosome, 'gen_x'])
        
    def run(self, plot=False, linkage_loosening = 1):
        for chromosome, chr_data in self.groupby(level=0):
            if not chromosome in self.HMM:
                continue
            print('calculating %s' % chromosome, file=sys.stderr)
            self.loc[chromosome,'p_slct'], _ = self.HMM[chromosome].get_model_likelihood(prior='Psel', linkage_loosening=linkage_loosening, fresh=True)
            self.loc[chromosome,'p_stat'], _ = self.HMM[chromosome].get_model_likelihood(prior='Pstat', linkage_loosening=linkage_loosening, fresh=True)
            self.loc[chromosome, 'lod'] = self.loc[chromosome, 'p_slct'] - self.loc[chromosome, 'p_stat']
        
        if plot:
            self.plot_both()
    
    def plot_chr(self, y = ['p_stat', 'p_slct'], peaks = False):
#        ax.plot(self['Mb'] , self['p_slct'], 'r.-', label='selection')
#        ax.plot(self['Mb'] , self['p_stat'], 'g-', label='no selection')
        
        def find_peaks(y, tolerance = None):
            ddy = np.diff(y, 2)
            signs = np.diff(np.sign(np.diff(y) )) !=0
            if tolerance is None:
                tolerance = np.median(np.abs(ddy[signs])) / 8
            return np.concatenate(([False,], \
            signs * (ddy < -tolerance), [False,]) )
        
        def plot_peaks(x_, y_, xlab_ = None, ax = plt.gca()):
            if type(xlab_) == type(None):
                xlab_ = x_
            "find peaks"
            ind = find_peaks(y_)
            ax.plot(x_[ind] , y_[ind], 'o', markersize = 8, markerfacecolor = 'none',\
            markeredgecolor = 'r')
            for n in range(len(ind)):
                if ind[n]:
                    print('peak: %u' % n, file=sys.stderr)
                    x, y, xl = x_[n] , y_[n], xlab_[n] 
                    ax.text(x, y, "{:,}".format(xl), style='italic')
            print('number of peaks: %u' % sum(ind), file=sys.stderr)
            return
        
        fig = plt.figure()
        fig.suptitle( (', '.join(y) if type(y) is list else y) + ' for all contigs' )
        nsubpl = len(self.HMM)
        nn = [1]
        ax = []
        for chromosome, chr_data in self.groupby(level=0):
            if not chromosome in self.HMM:
                continue
            ax.append( fig.add_subplot(nsubpl, 1, nn[0]) )
            nn[0] += 1
            chr_data.plot(ax = ax[-1], x = 'Mb', y = y, 
               label = y)
            if peaks:
                plot_peaks(np.array(chr_data['Mb']),
                           np.array(chr_data[y[-1] if type([1,2]) is list else y]),
                           np.array(chr_data['x']), 
                           ax = ax[-1])
            
        ax[-1].set_xlabel('position, Mb', fontsize=12)
        for aa in ax[:-1]:
            aa.legend('', frameon = False)
        fig.show()
        print('===================================' , file=sys.stderr)
        return fig, ax