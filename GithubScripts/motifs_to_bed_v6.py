# -*- coding: utf-8 -*-
"""
Created on Wed May  8 12:35:58 2024

@author: tsjit
Update 5-6-2024 v5,  Anne
Update 7-6-2024 v6, Tsjit

module load Biopython/1.81-foss-2022b
python3 /scratch/hb-molgen/TsIsZiJo/tsjits/scripts/motifs_to_bed_v6.py

"""

#----------------------------------MODULES-------------------------------------
import pandas as pd
import numpy as np
import argparse
import re 
import os

#------------------------------PARSE PARAMETERS--------------------------------
parser = argparse.ArgumentParser(description='test')
parser.add_argument('-sessiondir', dest='sessiondir', help='Full path sessiondir', nargs='?', default='/scratch/hb-molgen/TsIsZiJo/Anne/collecttf_meme')
parser.add_argument('-meme', dest='meme', help='Meme motifs results as fasta file', nargs='?', default='Bacillus_subtilis_subsp_subtilis_str_168_ASM904v1.MEME_results.fna')
parser.add_argument('-inter', dest='inter', help='Intergenic regions table as .gff', nargs='?', default='Bacillus_subtilis_subsp_subtilis_str_168_ASM904v1.intergenic.gff')
parser.add_argument('-out', dest='outputfile', help='Bed result file', nargs='?', default='Bacillus_subtilis_subsp_subtilis_str_168_ASM904v1.meme_motifs.bed')
parser.add_argument('--version', action='version', version='Tsjit Kramer, version 6.0, June 2024')
args = parser.parse_args()


#----------------------------------MAIN CODE-----------------------------------
# Regular expression pattern to extract data
pattern = r'>(.*\d+)_site_1 offset= (\d+)( RC)?\n(.*)'
pattern = r'>(.*\d+)_.* offset= (\d+)( RC)?\n(.*)'

# Read MEME data from the file
with open(args.sessiondir+'/'+args.meme, 'r') as file: MEME = file.read()

# Read Intergenic regions table
intergenic_df = pd.read_csv(args.sessiondir+'/'+args.inter, sep='\t')

# List to store extracted MEME fasta data
extracted_data = []

# Iterate over each line of data
for match in re.finditer(pattern, MEME):
	locus_tag, offset, rc, sequence = match.groups()
	strand = '-' if rc else '+'
	extracted_data.append((locus_tag, offset, strand, sequence))

# Create DataFrame
MEME_df = pd.DataFrame(extracted_data, columns=['locus_tag', 'Offset', 'strand', 'Sequence'])
MEME_df.replace([np.inf, -np.inf], np.nan, inplace=True) # Replace infinite values with NaN

MEME_df = pd.merge(MEME_df, intergenic_df[['locus_tag', 'strand']], on='locus_tag', how='left')
MEME_df.drop(columns=['strand_x'], inplace=True)
MEME_df.rename(columns={'strand_y':'strand'}, inplace=True)

print(MEME_df)

# Iterate over rows of MEME_df
for index, row in MEME_df.iterrows():
    locus_tag = row['locus_tag']
    offset = row['Offset']
    strand = row['strand']
   
    # Check if locus_tag exists in intergenic_df
    if locus_tag in intergenic_df['locus_tag'].values:
        # Find corresponding locus_tag in intergenic_df
        intergenic_row = intergenic_df[intergenic_df['locus_tag'] == locus_tag].iloc[0]
       
        # Convert 'end' to int
        end = int(intergenic_row['end'])
       
        # Convert 'offset' to int
        offset = int(offset)
       
        # Calculate new value based on strand
        if strand == '+':
            new_start = intergenic_row['start'] + offset
            new_end = intergenic_row['start'] + offset + len(row['Sequence'])
        else:
            new_end = end - offset
            new_start = end - offset - len(row['Sequence'])
       
        # Assign new value to new column
        MEME_df.at[index, 'start'] = new_start
        MEME_df.at[index, 'end']   = new_end
    else:
        print(f"Locus tag {locus_tag} not found in intergenic_df.")


# Drop rows with NaN values in 'start' or 'end'
MEME_df.dropna(subset=['start', 'end'], inplace=True)
MEME_df['start'] = MEME_df['start'].astype(int)
MEME_df['end'] = MEME_df['end'].astype(int)
MEME_df['Offset'] = MEME_df['Offset'].astype(int)
MEME_df['score']  = 1
MEME_df['chrom']  = intergenic_df['chrom']

print(MEME_df)

# Saving the output file 
my_bed_columns=['chrom','start','end','score','strand']
MEME_df[my_bed_columns].to_csv(args.sessiondir+'/'+args.outputfile , sep = '\t', index = False)


