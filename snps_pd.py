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

class snps_pd():
    membership_set = False
    @property
    def _constructor(self):
        return snps_pd
    
    def __init__(self, file_path, experiment_name = None, pop_size = None, \
            conditions = '', lengths_path='reference/TAIR10-chr-counts.dat'):
        columns ='chr, pos, altCount, totCount'
        db_table = GenomeVariantDbReader(file_path, experiment_name, subtable='mtCounts')
#        db_table.count_rows()
        self.data = pd.DataFrame( pd.read_sql_query( 'SELECT %s ' % columns + \
            db_table._from_where_(None, conditions),  
            db_table.connection) )
        self.db_table = db_table
        self.pop = population(pop_size);

        self.snp_reads_on_chrs = {}
        self.gen_map = GeneticMap()
        self._set_chromosome_lengths_(lengths_path)
        
        self.data.rename(columns={'pos': 'x', 'totCount':'r', 'altCount': 'q'}, inplace=True)
        self.data['Mb'] = self.data['x']*1e-6
        self.data.set_index(['chr'], inplace=True)
        "if index is unsorted, assignment of a new column in `unmix` will fail"        
        self.sort_index(inplace = True)
        return
#        print( '==========')
#        for chromosome, chr_data in self.groupby(level=0):
#            print( '%s', '%u entries' % (chromosome, chr_data.shape), sep = '\t', file=sys.stderr)
#        print( '==========')
    def __getitem__(self, name):
        return self.data[name]
    def __getattr__(self, name):
        return self.data.__getattribute__(name)
    
    def _set_chromosome_lengths_(self, path=''):
        if len(path) == 0:
            return
        else:
            self.chromosome_lengths = ChromosomeLength(path)

    def unmix(self):        
        for chromosome, chr_data in self.data.groupby(level=0):
            print('unmixing %s' % chromosome, file=sys.stderr)
            dx_raw = np.diff(chr_data['x'])
            
            self.data.loc[chromosome, 'dx'] = np.min([np.concatenate(([np.Inf], dx_raw)),  np.concatenate((dx_raw, [np.Inf] ))], axis=0)
            self.data.loc[chromosome,'log10_dx'] = np.log10(self.data.loc[chromosome, 'dx'])
            self.data.loc[chromosome,'membership'], self.mixt_obj, iTrue, out_mode_number, _ = \
                unmix( self.data.loc[chromosome,'log10_dx'] )
        
        self.membership_set = True
        
    @property
    def contigs(self):
        return self.db_table.contigs

    def prepare(self):
        if not self.membership_set:
            print("membership not set!", file=sys.stderr)
            self.unmix()
        self.HMM = {} 
        for chromosome, chr_data in self.data.groupby(level=0):
            "genetic distance for given data points : `gen_x`"
            self.data.loc[chromosome, 'gen_x'] = self.gen_map(chr_data['x'], chromosome)
            if   type(self.data.loc[chromosome, 'gen_x']) is pd.core.series.Series and \
            not np.isnan(self.data.loc[chromosome, 'gen_x'].iloc[0]):
                print('preparing emission and transition for %s' % chromosome, file=sys.stderr)
                "skip contigs with no genetic map (which is set to `NaN`)"
                emission_array = binomial_uniform_emission(self.pop, \
                    chr_data['q'], chr_data['r'], chr_data['membership'])
                self.HMM[chromosome] = hmm_cont(self.pop, emission_array, self.data.loc[chromosome, 'gen_x'])
        
    def run(self, plot=False, linkage_loosening = 1):
        for chromosome, chr_data in self.data.groupby(level=0):
            if not chromosome in self.HMM:
                continue
            print('calculating %s' % chromosome, file=sys.stderr)
            self.data.loc[chromosome,'p_slct'], _ = self.HMM[chromosome].get_model_likelihood(prior='Psel', linkage_loosening=linkage_loosening, fresh=True)
            self.data.loc[chromosome,'p_stat'], _ = self.HMM[chromosome].get_model_likelihood(prior='Pstat', linkage_loosening=linkage_loosening, fresh=True)
            self.data.loc[chromosome, 'lod'] = self.data.loc[chromosome, 'p_slct'] - self.data.loc[chromosome, 'p_stat']
        
        if plot:
            self.plot_both()
    
