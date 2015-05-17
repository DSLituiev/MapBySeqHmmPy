# -*- coding: utf-8 -*-
"""
Created on Mon Dec 22 14:46:15 2014

@author: dima
"""

import sys

sys.path.append('/home/dima/data/repos/snpdb/')
from snp_tables import *

from unmix_repeats import *

dbPath = '../../snpdb/vcf.db'
dbName = 'HL7__mtCounts'
st = GenomeVariantDbReader(dbPath, dbName)


#res = st.selectPosition()
chromosome = 'Chr1'
pos = '*'
res = st.selectPosition(chromosome, pos)

# print(res)
# print(res[2])

import numpy as np

# res_new = [ [int(line[0][-1]),] + list(line[1:]) for line in res]
res_new = [ list(line[1:]) for line in res]

x = np.array(res_new, dtype=int)
dx, log10_dx = min_distances(x)

Membership, mixtObj, iTrue, out_mode_number, _ = unmix(log10_dx)

import matplotlib.pyplot as plt

plt.figure()
plt.plot(x[:, 0], log10_dx)

plt.figure()
plt.plot(log10_dx, Membership, 'r.')