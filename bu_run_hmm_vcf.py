# -*- coding: utf-8 -*-
"""
Created on Mon Jan  5 16:08:03 2015

@author: dima
"""
"""
ISSUES
======
 - works with sqlite, not vcf immediately
 - not clear if the results are the same as in MATLAB
 - needs refactoring (collect high-level functions)
 
"""

from snp_reads_data import *
import matplotlib.pyplot as plt
###############################################################################

###############################################################################
plt.close('all')

dbPath = '../../snpdb/vcf.db'
dbName = 'HL7'
N = 50;
linkage_loosening = 16

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
andjoin = ' AND '.join

cc = 1

count_conditions = ['(totCount >=10)', ' (altCount >= 7) ']
sr.read_chromosome(cc, andjoin(count_conditions))

sr.db_table.countRows(count_conditions)

sr[cc].add_membership()
sr[cc].plot_snp_density(cc, threshold=2, grid=False);

print( '==========')
print( '%u entries' % sr[cc]['x'].shape )
print( '==========')

count_conditions = ['(totCount >=10)', ' (altCount >= 7) ',
    '(CAST(altCount AS FLOAT) / CAST(totCount AS FLOAT) < 2.0/3.0)', 
    '(CAST(altCount AS FLOAT) / CAST(totCount AS FLOAT) > 1.0/8)', ]

sr.read_chromosome(cc, andjoin(count_conditions) )
sr[cc].column_exists('membership')

print( '==========')
print( '%u entries' % sr[cc]['x'].shape )
print( '==========')
###############################################################################
sr[cc].prepare()
sr[cc].data['membership'].plot(kind='box')
sr[cc].data['q'].plot(kind='box')

"here it takes time:"
p_out, _ = sr[cc].HMM.get_model_likelihood(prior='Psel', linkage_loosening = linkage_loosening)
###############################################################################

fig_sel = plt.figure()
ax = fig_sel.add_subplot(111)
fig_sel.suptitle('probability of selection')
ax.plot(sr[cc]['Mb'], p_out, '.-')
plt.show()

###############################################################################
###############################################################################
crowded_inds, _, _ = sr[cc].threshold_dense_snps(1, threshold=2, grid=False)
"set membership of crowded snps to `0`"
sr[cc].data['membership'][crowded_inds] = 0
sr[cc].prepare()

sr[cc].run(linkage_loosening=linkage_loosening)

sr[cc].data.plot(x = 'Mb', y = ['p_slct', 'p_stat'], style = ['r.-',  'b.-'])

def plot_both(self):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.suptitle('probability')
    ax.plot(self['Mb'] , self.data['p_slct'], 'r.-', label='selection')
    ax.plot(self['Mb'] , self.data['p_stat'], 'g-', label='no selection')
    
    def peaks(lod, tolerance = 0):
        return np.concatenate(([False,], \
        (0!=np.diff(np.sign(np.diff(lod) ))) * (np.diff(lod, 2) < -tolerance), [False,]) )
    
    def plot_peaks(x_, y_, xlab_ = None, ax = plt.gca()):
        if type(xlab_) == type(None):
            xlab_ = x_
        "find peaks"
        ind = peaks(y_)
        ax.plot(x_[ind] , y_[ind], 'ro-', markersize = 8)
        for n in range(len(ind)):
            if ind[n]:
                print('peak: %u' % n)
                x, y, xl = x_[n] , y_[n], xlab_[n] 
                ax.text(x, y, "{0:n}".format(xl), style='italic')
        print('number of peaks: %u' % sum(ind))
        return
        
    plot_peaks(self.data['Mb'], self.data['p_stat'], self.data['x'], ax = ax)
    plt.show()
    return fig

sr[cc].data.reset_index(drop=True, inplace=True)
plot_both(sr[cc])

###############################################################################
p_sel, _ = sr[cc].HMM.get_model_likelihood(prior='Psel', linkage_loosening=linkage_loosening, fresh=True)
p_stat, _ = sr[cc].HMM.get_model_likelihood(prior='Pstat', linkage_loosening=linkage_loosening, fresh=True)
lod = p_sel - p_stat

fig_sel = plt.figure()
ax = fig_sel.add_subplot(111)
fig_sel.suptitle('probability of selection')
ax.plot(sr[cc]['Mb'], p_sel, 'r.-', label='P[sel] without crowded')
ax.plot(sr[cc]['Mb'], p_stat, 'cx-', label='P[stat] without crowded')
plt.show()

fig_sel = plt.figure()
ax = fig_sel.add_subplot(111)
fig_sel.suptitle('log-odds of selection')
ax.plot(sr[cc]['Mb'][[0,-1]] , np.array([0,0]) , 'k--')
ax.plot(sr[cc]['Mb'] , lod, 'r+-', label='log-oddsm without crowded')
ind = np.concatenate(([False,], \
    (0!=np.diff(np.sign(np.diff(lod) ))) * (np.diff(lod, 2)<-1e-1 ), [False,]) )
#ax.plot(sr[cc]['x']* 1e-6, ind , 'b+-')
ax.plot(sr[cc]['Mb'][ind], lod[ind], 'mx', label='log-oddsm without crowded')
###########
# import locale
# locale.setlocale(locale.LC_ALL, '')
###########3
for n in range(len(lod)):
    if ind[n]:
        x, y = sr[cc]['x'][n] , lod[n]
        ax.text(x* 1e-6, y, "{0:n}".format(x), style='italic')
plt.show()

###############################################################################
###############################################################################
p_flat = sr[cc].HMM.plot_raw_lh(linkage_loosening=linkage_loosening, fresh=True)

###############################################################################

sr[cc].data = sr[cc].data[~crowded_inds]
sr[cc].prepare()
###############################################################################
p_out,_ = sr[cc].HMM.get_model_likelihood(prior='Psel', linkage_loosening=linkage_loosening, fresh=True)
###############################################################################

fig_dx = plt.figure()

fig_dx.suptitle('log dx')
ax = fig_dx.add_subplot(111)
ax.plot(sr[cc]['Mb'], p_out, 'r*-')
fig_dx.show()

###############################################################################

fig1 = plt.figure()
fig1.suptitle('log dx')
plt.plot(sr[cc]['x'] * 1e-6, np.log10(sr[cc]['dx']))
plt.show()

fig2 = plt.figure()
plt.plot(np.log10(sr[cc]['dx']), sr[cc]['membership'], 'r.')
plt.show()
##########################################################
"play with pandas package tools for visualization"
import pandas as pd
disk_engine = sr.db_table.connection

df = pd.read_sql_query('SELECT * FROM "%s__vcf" WHERE (chr == "Chr%u")' % (dbName, cc), disk_engine)

#df['dx'] = 
df2 = pd.read_sql_query('SELECT * FROM "%s__mtCounts" WHERE (chr == "Chr%u")' % (dbName, cc), disk_engine)

df2.plot(x='pos', y = ['totCount', 'refCount',  'altCount'])
df2['f'] = df2['altCount'] / df2['totCount']
df2['Mb'] = df2['pos'] * 1e-6

ax = df2.plot(x='Mb', y = 'f', style = 'b.', label = 'all')
df2[df2['f']< 2/3].plot(x='Mb', y = 'f', style = 'r.', ax = ax )
lines, labels = ax.get_legend_handles_labels()
ax.legend(lines, ['all', '<2/3'], loc='best')  

plt.show()
df2['f'].plot(kind='box')
plt.show()
df2['f'].plot(kind='hist')
plt.show()

print('full \t: %u' % df2.shape[0] )
print( '2/3 \t: %u' % df2[df2['f']< 2/3].shape[0] )
print( '3/4 \t: %u' % df2[df2['f']< 3/4].shape[0] )


df2['totCount'].plot(kind='box')

df2['totCount'][df2['f']< 2/3].plot(kind='box')

result = pd.concat([df, df2], axis=1, join = 'inner')
result[0:10]
