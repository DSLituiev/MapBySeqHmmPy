# -*- coding: utf-8 -*-
"""
Created on Sat Jun  6 17:27:47 2015

@author: dima
"""

class ChromosomeLength:
    def __init__(self, path):
        self.dict = {}
        self.L = []
        self.chromosome_order = []
        with open(path) as fi:
            for line in fi:
                cols = line.split('\t')
                chr_length = int(cols[1])
                chr_name = cols[0].split(' ')[0]
                self.dict[chr_name] = chr_length
                self.L.append(chr_length)
                self.chromosome_order.append(chr_name)

    def __getitem__(self, chromosome):
        if type(chromosome) is int:
            return self.L[chromosome - 1]
        else:
            return self.dict[chromosome]

    def __repr__(self):
        st = 'chromosome lengths:\n'
        for nn, cc in enumerate(self.chromosome_order):
            st += '  \t'. join(['%u' % (nn + 1), '%s' % cc, '%u' % self.dict[cc] ]) + '\n'
        return st
        
    def __len__(self):
        return len(self.chromosome_order)
    
    def __floatlist__(self):
        return [float(self.dict[ss]) for ss in self.chromosome_order]

    __str__ = __repr__