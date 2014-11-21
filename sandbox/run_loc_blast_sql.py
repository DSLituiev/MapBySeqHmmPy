# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 14:00:11 2014

@author: dima
"""
import sys
import os
import xml.dom.minidom  
# import xml.etree.ElementTree as ET
from fasta_tools import *
###############################################################################
import argparse
parser = argparse.ArgumentParser('''
''')

parser.add_argument("inFile",
                   help="input file (.fasta/.fas) with the gene IDs given in the 'geneID' column")        

parser.add_argument("outFile", nargs='?', type=str, default='',
                    help="output file (.csv); for stdout type '-' ")

parser.add_argument("-r", "--reference_fasta", type=str, default="../reference/TAIR10/TAIR10",
                    help="")     

parser.add_argument("-o", "--species", type=str, default="Brapa_197",
                    help="")     
                    
parser.add_argument("-s", "--csvseparator", type=str, default= r';',
                    help="separator for the output .csv file [;]")
                             
args = parser.parse_args()

ORGANISM = args.species
###############################################################################
if not args.outFile == r'-':
    if not args.outFile == r'':
        sys.stdout = open(args.outFile, 'w')
    else:
        import re
        outFile = re.sub(r'\.fa.+', '-conservation_'+ ORGANISM + '.csv', args.inFile)
        assert(not outFile==args.inFile)
        sys.stdout = open(outFile, 'w')
        print('output file: %s' % outFile,  file=sys.stderr) 

###############################################################################
try:
    CURR_DIR = os.path.dirname(os.path.realpath(__file__))
except NameError:
    CURR_DIR = '.'
SNPDB_PATH = os.path.normpath(os.path.join(CURR_DIR, '..'))
print('adding path: %s' % SNPDB_PATH,  file=sys.stderr) 
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
    
    assert os.path.isfile( db_path ), ('no reference for blast found: %s' % db_path )
    cmd_blast = 'echo ' + fasta_str + ' | ' + \
    'blastn -task blastn-short -db ' + db_path  + ' -evalue 1e-5 -num_alignments 5 -outfmt "' + \
    FMT + ' " ' #  dc-megablast
    
    out, err = runcmd( cmd_blast  )
    return  (out, err)

#####################################################################
def get_xml_field_data(root, field):
    return root.getElementsByTagName(field)[0].childNodes[0].data

def check_hit_bounds(qdata, qfr, qto, LOC_IND):
    if (qfr > LOC_IND):
        raise IndexError
    
    jj = qfr
    while (jj < LOC_IND):
        if qdata[jj - qfr] == '-':
            LOC_IND+=1
            if (qto < LOC_IND):
                raise IndexError
        jj+=1
        
    if (qto < LOC_IND):
        raise IndexError
        
    final_ind = LOC_IND - (qfr-1)
    return final_ind
    
def compare_alignments(xml_string, LOC_IND):
    
    try:
        dom = xml.dom.minidom.parseString( xml_string )   
    except:
        print("empty blast result submitted", file=sys.stderr) 
        return (False, None, None, {})
        
    try:
        hits = dom.getElementsByTagName('BlastOutput_iterations')[0].getElementsByTagName('Iteration')[0].getElementsByTagName('Iteration_hits')[0].getElementsByTagName('Hit')
        
        hsp  = hits[0].getElementsByTagName('Hit_hsps')[0].getElementsByTagName('Hsp')[0]
        contig = get_xml_field_data(hits[0], 'Hit_id')   
        
    except IndexError as err:
        return (False, None, None,  {})

    hsp_fields = {}
    
    for child_n in hsp._get_childNodes():
        if (child_n.nodeType == child_n.ELEMENT_NODE):
            hsp_fields[child_n.tagName] = child_n.childNodes[0].data

    qfr = int(hsp_fields['Hsp_query-from'])
    qto = int(hsp_fields['Hsp_query-to'])
   
    
    try:
        final_ind = check_hit_bounds(hsp_fields['Hsp_qseq'], qfr, qto, LOC_IND)
    except IndexError:
        return (False, '', '', {})
    
    q_allele = hsp_fields['Hsp_qseq'][final_ind]
    mathched_flag = ( q_allele == hsp_fields['Hsp_hseq'][final_ind])
    
    return (mathched_flag, q_allele, contig, hsp_fields)

#####################################################################
def initializeSnpTable(conn, species, chrNumber = 5):
    
    dbases = [''] * chrNumber
    
    with conn:
        curs = conn.cursor()
        
        for cc in range(1,chrNumber+1):
            dbases[cc-1] = 'Conservation_%s_%u'% (species, cc)
            query = 'DROP TABLE IF EXISTS '+ dbases[cc-1]
            print(query, file=sys.stderr)
            curs.execute(query)
        
            query = ' CREATE TABLE ' + dbases[cc-1] + '''(
            pos   INT    PRIMARY KEY,
            genotype TEXT
            conserved TEXT
            )
            '''
            curs.execute(query)

#####################################################################            
# ORGANISM = 'Vitis_vinifera'
# ORGANISM = 'TAIR10'

col_names = ['chromosome', 'position', 'conserved', 'refAllele', 'contig', \
'hit_from', 'hit_to', 'score', 'e_value']

print(args.csvseparator.join(col_names))
    
fiter = fasta_iter(args.inFile)

for ff in fiter:
    headerStr, seq = ff
    (cc, pp) = headerStr.split(':')
    LOC_IND = int(len(seq)//2)+1
    ref_allele = seq[LOC_IND]
    out, err = run_blast(seq, ORGANISM)    
    out_lines = out.split('\n')
    
    """
    out_plain, err = run_blast(seq, ORGANISM, '1')    
    print(out_plain[182:] , file=sys.stderr) 
    print("______________________________________", file=sys.stderr) 
    """
    mathched_flag, q, contig, hsp_fields = compare_alignments(out, LOC_IND)
         
    for kk in hsp_fields:
        print(kk + '\t' + hsp_fields[kk] , file=sys.stderr) 
    print("______________________________________", file=sys.stderr) 
        
    assert (not mathched_flag) or (ref_allele == q), 'reference allele mismatch'
    if len(hsp_fields)>0:
        print(cc, pp, mathched_flag, ref_allele, contig, 
              hsp_fields['Hsp_hit-from'], hsp_fields['Hsp_hit-to'],
              hsp_fields['Hsp_score'], hsp_fields['Hsp_evalue'],
              sep = args.csvseparator)
    else:
         print(cc, pp, mathched_flag, ref_allele, contig, 
              sep = args.csvseparator)
#####################################################################

print('finished successfully',  file=sys.stderr) 
        
if not args.outFile == r'-':
    sys.stdout.close()
    
sys.stdout = sys.__stdout__