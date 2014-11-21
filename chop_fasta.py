import sys
import os
from runshell import *

###############################################################################
import argparse
parser = argparse.ArgumentParser('''
''')

parser.add_argument("inFile",
                   help="input file (.csv) with the gene IDs given in the 'geneID' column")        

parser.add_argument("outFile", nargs='?', type=str, default='',
                    help="output file (.csv); for stdout type '-' ")

parser.add_argument("-n", "--chunk_size", type=int, default=25,
                    help="")                    

parser.add_argument("-r", "--reference_fasta", type=str, default="../reference/TAIR10/TAIR10",
                    help="")     
                    
parser.add_argument("-s", "--csvseparator", type=str, default= r';',
                    help="separator for the output .csv file [;]")
                             
                    
args = parser.parse_args()
###############################################################################

if not args.outFile == r'-':
    if not args.outFile == r'':
        sys.stdout = open(args.outFile, 'w')
    else:
        sys.stdout = open(args.inFile.replace('.csv','.fas'), 'w')
        print('output file: %s' %args.inFile.replace('.csv','.fas'),  file=sys.stderr) 

###############################################################################
header = 'chromosome; position;\n'
snp_table = '1; 20;\n 1; 100;\n 1; 120;\n 1; 340;'.split('\n')

assert os.path.isfile( args.reference_fasta)
in_fas_file  = args.reference_fasta
in_bed_file  = 'test.bed'
out_fas_file = 'test.fa.out'

FLANK = 60
# FASTA_VIEW_CMD = "echo"

def chop_cmd(chrm_range):
    FASTA_VIEW_CMD = 'samtools faidx'
    curr_cmd = [FASTA_VIEW_CMD, in_fas_file, chrm_range]
    out, err = runcmd(' '.join(curr_cmd) )
    return (out, err)

with open(args.inFile) as f:
    header = next(f)
    splitHeader = header.split(args.csvseparator)
    colInds = {}
    
    for ii in range(0,len(splitHeader)):
        colInds[splitHeader[ii].split('/')[0].rstrip()] = ii
        
    # print(header.rstrip(';\n'), '', 'ratioPass', 'expressionPiresPass' , 'ratioAndExpressionPass', sep = args.csvseparator)
    # [elem for i, elem in enumerate(inputlist) if i not in excluded_indices]
    
    for line in f:     

        lineVals = line.split(args.csvseparator)
    
        chromosome = int(lineVals[0])
        chrom_str = 'Chr%u' % chromosome
    
        snp_position = int(lineVals[1])
        region_start = max(1, snp_position - args.chunk_size)
        region_end = snp_position + args.chunk_size
        # chr2:1,000,000-2,000,000
        chrm_range = chrom_str+':'+ '%u' % region_start + '-'+'%u' % region_end
    
        out, err = chop_cmd(chrm_range)
    
        out_lines = out.split('\n')
        out_lines[0] = '>'+chrom_str+':%u' % snp_position
        
        # print(out_lines, file=sys.stderr) 
        tot_size = 2*args.chunk_size + 1
        if len(out_lines[1]) < tot_size :
            out_lines[1] = 'N'* (len(out_lines[1]) - tot_size) + out_lines[1]
            
        out_decoded =  '\n'.join(out_lines) # python3
        print(out_decoded[:-1])

#####################################################################

print('finished successfully',  file=sys.stderr) 
        
if not args.outFile == r'-':
    sys.stdout.close()
    
sys.stdout = sys.__stdout__
