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
from plot_contigs import *

###############################################################################
###############################################################################
plt.close('all')

dbPath = '../../snpdb/vcf.db'
dbName = 'ABD159'
N = 50;
linkage_loosening = 4
###############################################################################
count_conditions = ['(totCount >=10)', ' (altCount > 7) ',
    '(CAST(altCount AS FLOAT) / CAST(totCount AS FLOAT) < 2.0/3.0)', 
    '(CAST(altCount AS FLOAT) / CAST(totCount AS FLOAT) > 1.0/8)' ]

sr = snps_pd(dbPath, dbName, N, conditions = count_conditions)

#####################################################################
sr.unmix()

sr.data = sr[ sr['r'] < sr['r'].quantile(0.98)]

sr.prepare()
sr.run()

#fig,_ = plot_contig_col(sr, y = 'r')

fig, _ = plot_contig_row(sr.data, x= 'x', y = 'lod', kind = 'connected', rescale_x = 1e-6, len_by_chr = sr.chromosome_lengths.dict)
plt.show()

sr.to_csv( 'output/' + dbName+'.csv')
