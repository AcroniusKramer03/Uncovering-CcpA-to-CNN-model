# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 21:35:07 2024

python /scratch/hb-molgen/TsIsZiJo/tsjits/scripts/intergenic_gff.py


@author: tsjit
"""
#----------------------------------MODULES-------------------------------------
import pandas as pd
import argparse

#------------------------------PARSE PARAMETERS--------------------------------
parser = argparse.ArgumentParser(description='test')
parser.add_argument('-inter', dest='inter', help='The regulon_intergenic file, output from motifs_to_bed', nargs='?', default='Path/to/bacteria.regulon_intergenic')
parser.add_argument('-out', dest='outdir', help='bed_file location + name, code will add .gff', nargs='?', default='Bacillus_subtilis_subsp_subtilis_str_168_ASM904v1.meme_motifs.bed')
parser.add_argument('--version', action='version', version='Acronius Tsjits Kramer, version 1.0, Jun 2024')
args = parser.parse_args()

#----------------------------------MAIN CODE-----------------------------------
# stores the filtered regulon_intergenic file created using the 'regulon_intergenic' script into the variable fasta_file
fasta_file = args.inter

# strips the regulon_intergenic and stores it in a list
intergenic_gff = []
with open(fasta_file, 'r') as f:
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            header = line[1:].split('|')
            locus_tag=header[0]
            chrom = header[1]
            start = int(header[2])
            end = int(header[3])
            strand = header[4]
        else:
            sequence = line 
            intergenic_gff.append([locus_tag, chrom, start, end ,strand])

# converts list to dataframe
GFF_df = pd.DataFrame(intergenic_gff, columns=['locus_tag', 'chrom', 'start', 'end', 'strand'])

# as a quality control prints the dataframe
print(GFF_df)

# converts dataframe into a tab seperated csv file
GFF_df.to_csv(args.outdir+'.gff', sep='\t', index=False)