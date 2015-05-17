# -*- coding: utf-8 -*-
"""
Created on Mon Jan  5 16:08:03 2015

@author: dima
"""

from snp_reads_data import *
import matplotlib.pyplot as plt
###############################################################################

###############################################################################
plt.close('all')

dbPath = '../../snpdb/vcf.db'
dbName = 'HL7'
N = 50;
linkage_loosening = 32

###############################################################################

def generate_pos_query_condition(intervals_to_exclude):
    cond = []
    for line in intervals_to_exclude:
        cond.append('NOT ( pos >= %u AND pos <= %u )' % tuple(line))

    return ' AND '.join(cond)


# intervals = np.zeros((2, 4))
intervals = np.array([[123, 5445], [6443, 6890], [12000, 12453], [18000, 20000]])

pos_cond = generate_pos_query_condition(intervals)

###############################################################################
sr = SnpReadsData(dbPath, dbName, N)
sr.check_data_tables()

cc = 2
sr.read_chromosome(cc, '(totCount >=5) AND (altCount >= 3)')
sr[cc].add_membership()
sr[cc].plot_snp_density(cc, threshold=2, grid=False);

sr.read_chromosome(cc,
                   '(totCount >=10) AND (altCount >= 7) AND (CAST(altCount AS FLOAT) / CAST(totCount AS FLOAT) < 0.8)')
sr[cc].column_exists('membership')

###############################################################################
sr[cc].prepare()
sr[cc].data['membership']

"here it takes time:"
p_out = sr[cc].HMM._runFB_(linkage_loosening, prior='Psel')
###############################################################################

fig = plt.figure()
ax = fig.add_subplot(111)
fig.suptitle('probability of selection')
ax.plot(sr[cc]['x'] * 1e-6, p_out)

###############################################################################
###############################################################################
y, _, _ = sr[cc].threshold_dense_snps(1, threshold=2, grid=False)
"set membership of crowded snps to `0`"
sr[cc].data['membership'][y] = 0
sr[cc].prepare()

###############################################################################
p_out = sr[cc].HMM._runFB_(linkage_loosening, prior='Psel', fresh=True)

fig = plt.figure()
ax = fig.add_subplot(111)
fig.suptitle('probability of selection')
ax.plot(sr[cc]['x'] * 1e-6, p_out, 'r*-', label='without crowded')
plt.show()
###############################################################################
###############################################################################

sr[cc].data = sr[cc].data[~y]
sr[cc].prepare()
###############################################################################
p_out = sr[cc].HMM._runFB_(linkage_loosening, prior='Psel', fresh=True)
###############################################################################

plt.figure()
plt.plot(sr[cc]['x'] * 1e-6, p_out, 'r*-')
plt.show()

###############################################################################

fig1 = plt.figure()
fig.suptitle('log dx')
plt.plot(sr[cc]['x'] * 1e-6, np.log10(sr[cc]['dx']))
plt.show()

fig2 = plt.figure()
plt.plot(np.log10(sr[cc]['dx']), sr[cc]['membership'], 'r.')
plt.show()
