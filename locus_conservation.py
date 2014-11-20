import sys
# import os
import subprocess


header = 'chromosome; position;\n'
snp_table = '1; 20;\n 1; 100;\n 1; 120;\n 1; 340;'.split('\n')
in_fas_file  = './reference/TAIR10-chr.fas'
in_bed_file  = 'test.bed'
out_fas_file = 'test.fa.out'

FLANK = 60
FASTA_VIEW_CMD = 'samtools faidx'
# FASTA_VIEW_CMD = "echo"

def runcmd(x):
   assert (type(x) == str), 'not a string!'
   p = subprocess.Popen( x , shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   return p

for line in snp_table:
    lineVals = line.split(';')
    # print(lineVals)

    chromosome = int(lineVals[0])
    chrom_str = 'Chr%u' % chromosome

    snp_position = int(lineVals[1])
    region_start = max(1, snp_position - FLANK)
    region_end = snp_position + FLANK
    # chr2:1,000,000-2,000,000
    chrm_range = chrom_str+':'+ '%u' % region_start+ '-'+'%u' % region_end
    curr_cmd = [FASTA_VIEW_CMD, in_fas_file, chrm_range]

    p = runcmd(' '.join(curr_cmd) )
    out, err = p.communicate()

    out_lines = out.decode('utf-8').split('\n')
    out_lines[0] = '>'+chrom_str+':%u' % snp_position

    out_decoded =  '\n'.join(out_lines) # python3
    print(out_decoded)


ORGANISM = 'Vitis'

ORGANISM = 'Arabidopsis'
FMT = '"2 scomnames stitle"'
fasta_str = ' "CATTTGTTATATTGGATACAAGCTTTGCTACGATCTACATTTGGGAATGTGAGTCTCTTA" '
cmd_blast = 'echo '+ fasta_str + ' | ' + \
'blastn -remote -db nr -entrez_query "' + ORGANISM + '[Organism]" -evalue 1e-20 -num_alignments 5 -outfmt ' + FMT
p = runcmd( cmd_blast  )
out, err = p.communicate()
out_lines = out.decode('utf-8').split('\n')
print('\n'.join(out_lines))