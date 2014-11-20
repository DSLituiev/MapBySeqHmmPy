# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 14:00:11 2014

@author: dima
"""
import sys
import os
import xml.etree.ElementTree as ET

import xml.dom.minidom  
#######################
try:
    CURR_DIR = os.path.dirname(os.path.realpath(__file__))
except NameError:
    CURR_DIR = '.'
SNPDB_PATH = os.path.normpath(os.path.join(CURR_DIR, '..'))
print('adding path: %s' % SNPDB_PATH)
sys.path.append(SNPDB_PATH)
from runshell import *
#######################
REFPATH = "/home/dima/scripts/reference"

def add_path_cmd(species_folder):
    return 'export BLASTDB="$BLASTDB:'+ os.path.join(REFPATH, species_folder) + '" '

species_list = ['TAIR10', 'Brassica_rapa', 'Solanum_lycopersicum', 'Vitis_vinifera']

for sl in species_list:
    out, err = runcmd( add_path_cmd(sl)  )

#####################################################################
def run_blast(fasta_str, genome_db_name, FMT = '5', REFPATH = "/home/dima/scripts/reference"):
    db_path = os.path.normpath(os.path.join(REFPATH, genome_db_name, genome_db_name))
    cmd_blast = 'echo ' + fasta_str + ' | ' + \
    'blastn -db ' + db_path  + ' -evalue 1e-3 -num_alignments 5 -outfmt "' + \
    FMT + ' " '
    out, err = runcmd( cmd_blast  )
    return (out, err)
#####################################################################
# ORGANISM = 'Vitis_vinifera'
ORGANISM = 'TAIR10'

fasta_str = ' "CATTTGTTATATTGGATACAAGCTTTGCTACGATCTACATTTGGGAATGTGAGTCTCTTA" '
out, err = runcmd( "echo $BLASTDB" )

out, err = run_blast(fasta_str, ORGANISM)

out_lines = out.decode('utf-8').split('\n')
print('\n'.join(out_lines))

print(err.decode('utf-8'))
################################3

"""
echo  "CATTTGTTATATTGGATACAAGCTTTGCTACGATCTACATTTGGGAATGTGAGTCTCTTA"  | blastn -db nr -entrez_query "Vitis[Organism]" -evalue 1e-3 -num_alignments 5 -outfmt "2 scomnames stitle"

echo  "CATTTGTTATATTGGATACAAGCTTTGCTACGATCTACATTTGGGAATGTGAGTCTCTTA"  | blastn -db TAIR10 -num_alignments 5 -outfmt "5"

echo  "CATTTGTTATATTGGATACAAGCTTTGCTACGATCTACATTTGGGAATGTGAGTCTCTTA"  | blastn -db /home/dima/scripts/reference/TAIR10/TAIR10  -num_alignments 5 -outfmt "5"
"""

################################3
def compare_alignments(out):
    dom = xml.dom.minidom.parseString( out.decode('utf-8') )   
    hsp = dom.getElementsByTagName('BlastOutput_iterations')[0].getElementsByTagName('Iteration')[0].getElementsByTagName('Iteration_hits')[0].getElementsByTagName('Hit')[0].getElementsByTagName('Hit_hsps')[0].getElementsByTagName('Hsp')[0]
    q = hsp.getElementsByTagName('Hsp_qseq')[0]
    qdata = q.childNodes[0].data
    
    h = dom1.getElementsByTagName('Hsp_hseq')[0]
    hdata = h.childNodes[0].data
    
    return hdata[-30] == qdata[-30]


################################3