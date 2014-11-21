# -*- coding: utf-8 -*-
"""
Created on Fri Nov 21 16:26:50 2014

@author: dima
"""

def dummy4():
    a = 1
    b = 2
    c = 'c'
    # d = {'alpha':121, 'beta': 'yes'}
    d = {'Hsp_hseq': 'TCCCACAAATGATCAATTAATGCATT', 'Hsp_hit-to': '23704228', 'Hsp_identity': '26', 'Hsp_query-from': '76', 'Hsp_midline': '||||||||||||||||||||||||||', 'Hsp_num': '1', 'Hsp_qseq': 'TCCCACAAATGATCAATTAATGCATT', 'Hsp_gaps': '0', 'Hsp_positive': '26', 'Hsp_hit-frame': '1', 'Hsp_score': '26', 'Hsp_bit-score': '52.0340900999248', 'Hsp_align-len': '26', 'Hsp_hit-from': '23704203', 'Hsp_query-frame': '1', 'Hsp_query-to': '101', 'Hsp_evalue': '5.09586899086248e-06'}
    
    return (a,b,c,d)

mathched_flag,q,contig,hsp_fields = dummy4()

print(mathched_flag,q,contig)

hsp_fields