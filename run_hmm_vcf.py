# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 23:00:51 2015

@author: dima
"""
"""
ISSUES
======
 - works with sqlite, not vcf immediately
 - not clear if the results are the same as in MATLAB
 - needs refactoring (collect high-level functions)
 
"""

import matplotlib.pyplot as plt
from snps_pd import *
import sys
from plot_contigs import plot_contig_row
 
###############################################################################
###############################################################################
plt.close('all')

dbPath = '../../snpdb/vcf.db'
dbName = 'HL7'
N = 50;
linkage_loosening = 4
###############################################################################
count_conditions = ['(totCount >=10)', ' (altCount >= 7) ',
    '(CAST(altCount AS FLOAT) / CAST(totCount AS FLOAT) < 2.0/3.0)', 
    '(CAST(altCount AS FLOAT) / CAST(totCount AS FLOAT) > 1.0/8)', ]

sr = snps_pd(dbPath, dbName, N, conditions = count_conditions)

sr = sr[ sr['r'] < sr['r'].quantile(0.8)]


sr.unmix()
sr.prepare()
sr.run()

sr.plot_chr(['lod'])

fig, _ = plot_contig_row(sr, x= 'x', kind = 'connected', rescale_x = 1e-6, len_by_chr = sr.chromosome_lengths.dict)
plt.show()

#def show_only_some(x, pos):
#    x_chr = chr_break_points[np.argmax( chr_break_points > x)] - x
#    return x_chr

#ax.xaxis.set_minor_formatter(plt.FuncFormatter(show_only_some))
#
################################################################################
#sr = SnpReadsData(dbPath, dbName, N)
#sr.check_data_tables()
#andjoin = ' AND '.join
#
#cc = 1
#count_conditions = ['(totCount >=10)', ' (altCount >= 7) ',
#    '(CAST(altCount AS FLOAT) / CAST(totCount AS FLOAT) < 2.0/3.0)', 
#    '(CAST(altCount AS FLOAT) / CAST(totCount AS FLOAT) > 1.0/8)', ]
#
#sr.read_chromosome(cc, andjoin(count_conditions) )
#sr[cc].column_exists('membership')
#
#print( '==========')
#print( '%u entries' % sr[cc]['x'].shape )
#print( '==========')
#
#" DIAGNOSTICS"
#sr[cc].data.plot(y = 'log10_dx', kind = 'hist', bins = 20)
#plt.show()
#sr[cc].data.plot(y = 'r', kind = 'hist', bins = 20, label = 'r')
#
################################################################################
#max_r = sr[cc].data['r'].quantile(1-1/16)
#sr[cc].data = sr[cc].data[sr[cc].data['r'] < max_r]
#sr[cc].data.reset_index(drop=True, inplace=True)
#sr[cc].prepare()
##sr[cc].data['membership'].plot(kind='hist')
#
#sr[cc].run(linkage_loosening=linkage_loosening)
#
#sr[cc].plot_both()
#
#sr.db_table.countRows(count_conditions)
#
#sr.db_table.contigs
#
################################################################################
 