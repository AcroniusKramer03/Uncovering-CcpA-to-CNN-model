'''
Author: Anne de Jong 

Prepare sequence list from bed or gff 
2024 April 17

Needed
1) my_genomename.fna  # the DNA sequence of the genome
2) my_genomename.bed  # the fragments to be extracted
3) ListOfGenomes.txt is a file with all full path genomes names to be process. Without the .fna or .bed extension
	this tab delimited file has 4 columns. e.g.,
ID	Name	Genome	FilePrefix
ID001	Acinetobacter_baumannii	Acinetobacter_baumannii	/userapps/hb-molgen/ProPr/examples/Acinetobacter_baumannii_ATCC_17978_ASM1542v1_genomic

4) To run the script
module load Biopython/1.81-foss-2022b

python3 "/scratch/hb-molgen/TsIsZiJo/tsjits/scripts/ProPr_Bed2Seq_v2.py"

'''

import os
import sys
import pandas as pd
import numpy as np
import argparse
import re
import random
import time

random.seed(0)

# ---------------------------------------------------------------- parse parameters ---------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Extract sequences based on bed. Bed file should contain at least the following headers; chrom, start, end, strand")
parser.add_argument('-bed', dest='bed', help='bed file', nargs='?', default='/scratch/hb-molgen/TsIsZiJo/Anne/collecttf_meme/Bacillus_subtilis_subsp_subtilis_str_168_ASM904v1.meme_motifs.bed')
parser.add_argument('-fna', dest='fna', help='multi fasta DNA file', default='/scratch/hb-molgen/TsIsZiJo/Anne/collecttf_meme/Bacillus_subtilis_subsp_subtilis_str_168_ASM904v1_genomic.fna')
parser.add_argument('-out', dest='out', help='output filename', default='/scratch/hb-molgen/TsIsZiJo/Anne/collecttf_meme/Bacillus_subtilis_subsp_subtilis_str_168_ASM904v1_genomic.sequences')
parser.add_argument('--version', action='version', version='Anne de Jong, version 1.0, June 2024')

args = parser.parse_args()


#bed_header = ["chrom","start", "end", 'name','score', "strand"]
#my_bed_columns=['chrom','start','end','locus_tag','score','strand', 'Sequence']


''' --------------------------------------  DEFINE FUNCTIONS ---------------------------------------------------'''

def readMultiFasta(fna_file):
	fasta = {}
	key= ''
	seq=''
	with open(fna_file) as f:
		for line in f:
			if line.startswith(">"):
				if (key != ''):   # chk previous key
					fasta[key] = seq.strip().upper()
					seq = ''
				items = re.match(">(.*?)\s", line)
				if items: key=items.group(1)
			else: seq += line.rstrip()
		if (key != ''): fasta[key] = seq.strip().upper()	# add last record	
	return fasta	



def reverse_complement(dna):
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	return ''.join([complement[base] for base in dna[::-1]])



''' ======================================   MAIN ================================================'''

FASTA = readMultiFasta(args.fna)

bed_df = pd.read_csv(args.bed, sep='\t', comment="#")
bed_df = bed_df.dropna()

print()
print("FASTA chroms:", FASTA.keys())
print("BED chroms:",bed_df['chrom'].unique())
print()

# makes all motifs to be padded to the left and right, making all motifs ~71 bases
sequenceSize=71

sequences = []
for index, row in bed_df.iterrows():
    len=row['end']-row['start'] 
    padding=int((sequenceSize-len)/2)
    chrom = row['chrom']
    start = row['start'] - padding
    end = row['end'] + padding
    sequence = FASTA[chrom][start:end]
    if row['strand']=='-': 
        sequence = reverse_complement(sequence)
        print("- strand")
    sequences.append(sequence)

bed_df['sequence'] = sequences
print(bed_df)

# store sequences in output file
if sequences:
	with open(args.out, 'w') as file: file.write("\n".join(sequences))
	print("List of sequences saved to "+args.out)
else:
	print("sequence list is empty")

