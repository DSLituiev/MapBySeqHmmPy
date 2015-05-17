# -*- coding: utf-8 -*-
"""
Created on Fri Nov 21 11:27:49 2014

@author: dima
"""

from itertools import groupby


def fasta_iter(fasta_name):
    """
    By Brent Pedersen
    https://www.biostars.org/p/710/
    given a fasta file. yield tuples of header, sequence
    """
    "first open the file outside "
    # fh = open(fasta_name)    
    with open(fasta_name) as fh:
        # ditch the boolean (x[0]) and just keep the header or sequence since
        # we know they alternate.
        faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
        for header in faiter:
            # drop the ">"
            headerStr = header.__next__()[1:].strip()
            # join all sequence lines to one.
            seq = "".join(s.strip() for s in faiter.__next__())
            yield (headerStr, seq)
            
