# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 14:00:11 2014

@author: dima
"""
import sys
# import os
from runshell import *

ORGANISM = 'Vitis'

#ORGANISM = 'Arabidopsis'
FMT = '"2 scomnames stitle"'
fasta_str = ' "CATTTGTTATATTGGATACAAGCTTTGCTACGATCTACATTTGGGAATGTGAGTCTCTTA" '
cmd_blast = 'echo '+ fasta_str + ' | ' + \
'blastn -db nr -entrez_query "' + ORGANISM + '[Organism]" -evalue 1e-3 -num_alignments 5 -outfmt ' + FMT
out, err = runcmd( cmd_blast  )
out_lines = out.decode('utf-8').split('\n')
print('\n'.join(out_lines))

print(err.decode('utf-8'))

"""
echo  "CATTTGTTATATTGGATACAAGCTTTGCTACGATCTACATTTGGGAATGTGAGTCTCTTA"  | blastn -db nr -entrez_query "Vitis[Organism]" -evalue 1e-3 -num_alignments 5 -outfmt "2 scomnames stitle"

echo  "CATTTGTTATATTGGATACAAGCTTTGCTACGATCTACATTTGGGAATGTGAGTCTCTTA"  | blastn -db TAIR10 -num_alignments 5 -outfmt "5"
"""