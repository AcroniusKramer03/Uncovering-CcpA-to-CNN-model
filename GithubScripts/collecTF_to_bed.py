# -*- coding: utf-8 -*-
"""
Created on Sat Jun 2024

python3 "/scratch/hb-molgen/TsIsZiJo/tsjits/scripts/collecTF_to_bed.py"

@author: tsjit
"""
#----------------------------------MODULES-------------------------------------
import pandas as pd
import argparse

#------------------------------PARSE PARAMETERS--------------------------------
parser = argparse.ArgumentParser(description='test')
parser.add_argument('-collectf', dest='collectf', help='The regulon_intergenic file, output from motifs_to_bed', nargs='?', default='Path/to/bacteria.regulon_intergenic')
parser.add_argument('-outdir', dest='outdir', help='bed_file location + name, code will add .gff', nargs='?', default='Bacillus_subtilis_subsp_subtilis_str_168_ASM904v1.meme_motifs.bed')
parser.add_argument('--version', action='version', version='Acronius Tsjits Kramer, version 1.0, Jun 2024')
args = parser.parse_args()

#----------------------------------MAIN CODE-----------------------------------
# opens the collectf downloaded fasta file, splits it up, extracts valuable information and stores it in a list 
collectf_bed = []

with open(args.collectf, 'r') as f:
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            header_parts = line[1:].split('|')
            locus_tag = header_parts[0]
            chrom = header_parts[1].replace('genome_', '').split()[0]
            start = int(header_parts[2].split('=')[1])
            end = int(header_parts[3].split('=')[1])
            strand = header_parts[4].split('=')[1]
            strand = '-' if strand == '-1' else '+' if strand == '1' else strand
            score = 1
        else:
            sequence = line
            collectf_bed.append([chrom, start, end, score, strand])
  
# turns fasta into dataframe          
collectf_df = pd.DataFrame(collectf_bed, columns=['chrom', 'start', 'end','score', 'strand'])
# turns dataframe into the tab seperated .bed file
collectf_df.to_csv(args.outdir+'collecTF.bed', sep = '\t', index = False)