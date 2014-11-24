#!/usr/local/bin/python3
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
import sqlite3
###############################################################################
import argparse
parser = argparse.ArgumentParser('''
''')

parser.add_argument("inFile",
                   help="input file (.fasta/.fas) with the gene IDs given in the 'geneID' column")        

parser.add_argument("-d", "--dbPath", default = '/home/dima/scripts/snpdb/vcf.db',
                    help="a path to the database")

parser.add_argument("-r", "--reference_fasta", type=str, default="../reference/TAIR10/TAIR10",
                    help="")     

parser.add_argument("-o", "--species", type=str, default="Brapa_197",
                    help="")     
                    
parser.add_argument("-s", "--csvseparator", type=str, default= r';',
                    help="separator for the output .csv file [;]")
                             
args = parser.parse_args()

ORGANISM = args.species

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
class conservation_tables():
    table_template = 'Conservation_%s_%u'
    colNum = 8
    
    
    def __init__(self, dbPath, species, chrNumber = 5, purge = False):
        self.conn = sqlite3.connect(dbPath)
        self.chrs_tables = [''] * chrNumber
        self.species = species
        if purge:
            self.insert_query = 'INSERT INTO ' + self.table_template + ' VALUES ('+ ','.join(['?']*self.colNum) +')'
        else:
            self.insert_query = 'INSERT OR IGNORE INTO ' + self.table_template + ' VALUES ('+ ','.join(['?']*self.colNum) +')'
        
        with self.conn:
            curs = self.conn.cursor()
            
            for cc in range(1,chrNumber+1):
                self.chrs_tables[cc-1] = self.table_template % (self.species, cc)
                
                if purge:
                    query = 'DROP TABLE IF EXISTS '+ self.chrs_tables[cc-1]
                    print(query, file=sys.stderr)
                    curs.execute(query)

                query = ' CREATE TABLE if not exists ' + self.chrs_tables[cc-1] + '''(
                pos   INT    PRIMARY KEY,            
                conserved INT,
                refAllele TEXT,
                contig TEXT,
                hit_from INT,
                hit_to INT,
                score INT,
                e_value REAL
                )
                '''
                curs.execute(query)
                
    def __del__(self):
        self.conn.close()
        
    def populateSnpTable(self, theIter):
         
        self.conn.isolation_level = None
        
        with self.conn:
            curs = self.conn.cursor()
            try:
                curs.execute("begin")
                for cons_obj in theIter:
                    entry = cons_obj.get_content()
                    assert len(entry) == self.colNum+1, 'number of columns mismatch'      
                    chromosome = int(entry[0][3])
                    query = self.insert_query % (self.species, chromosome)
                    curs.execute(query, entry[1:])
                curs.execute("commit")
            
            except self.conn.Error as ee:
                    print("failed to write to the database!")
                    curs.execute("rollback")
                    raise ee

#####################################################################
class conservation_from_fasta():
    _col_names_ = ['chromosome', 'position', 'conserved', 'refAllele', 'contig', \
            'hit_from', 'hit_to', 'score', 'e_value']
             
    def __init__(self, *args):
        if len(args)>=2:
            ff, ORGANISM = args
        else:
            return
        headerStr, seq = ff
        (self.cc, self.pp) = headerStr.split(':')
        LOC_IND = int(len(seq)//2)+1
        self.ref_allele = seq[LOC_IND]
        out, err = run_blast(seq, ORGANISM)    
        out_lines = out.split('\n')
        """
        out_plain, err = run_blast(seq, ORGANISM, '1')    
        print(out_plain[182:] , file=sys.stderr) 
        print("______________________________________", file=sys.stderr) 
        """
       
        self.mathched_flag, self.q, self.contig, self._hsp_fields_ = compare_alignments(out, LOC_IND)

        for kk in self._hsp_fields_:
            print(kk + '\t' + self._hsp_fields_[kk] , file=sys.stderr) 
        print("______________________________________", file=sys.stderr) 
            
        assert (not self.mathched_flag) or (self.ref_allele == self.q), 'reference allele mismatch'
        
        if len(self._hsp_fields_)>0:
            self.hit_from = self._hsp_fields_['Hsp_hit-from']
            self.hit_to   = self._hsp_fields_['Hsp_hit-to']        
            self.score    = self._hsp_fields_['Hsp_score']
            self.evalue   = self._hsp_fields_['Hsp_evalue']
        else:
            self.hit_from = 0
            self.hit_to = 0
            self.score = 0
            self.evalue = 1.0            
        return
        
    def get_content_order(self):
        return self._col_names_

    def get_content(self):
        return (self.cc, self.pp, int(self.mathched_flag), self.ref_allele, self.contig, 
                self.hit_from, self.hit_to, self.score, self.evalue)

    def print_content(self,  *args, **kwargs):
        print( *(self.get_content() + args), **kwargs)
                
#####################################################################            
class fasta_to_conserv():
    
    def __init__(self, inFile, ORGANISM):
        self.fasta_iter = fasta_iter(inFile)
        self.organism = ORGANISM

    def __iter__(self):
        return self
    
    def __next__(self):
        ff = self.fasta_iter.__next__()
        cons = conservation_from_fasta(ff, self.organism)
        return cons
#####################################################################
cons = conservation_from_fasta()
col_names = cons.get_content_order()

print(args.csvseparator.join(col_names))
    
fasta_iterator = fasta_to_conserv(args.inFile, ORGANISM)

cons_table = conservation_tables(args.dbPath, ORGANISM)
cons_table.populateSnpTable( fasta_iterator)

#####################################################################

print('finished successfully',  file=sys.stderr) 

sys.stdout = sys.__stdout__