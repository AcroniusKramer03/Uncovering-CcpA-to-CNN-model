# -*- coding: utf-8 -*-
"""
Created on Thu May 16 10:20:06 2024

@author: tsjit
"""

#----------------------------------MODULES-------------------------------------
import pandas as pd
import argparse
import os

#------------------------------PARSE PARAMETERS--------------------------------

parser = argparse.ArgumentParser(description='Unvocering regulon of genomes')
parser.add_argument('-faa', dest= 'proteins', help= 'Give up the faa file containing all proteins of a organism', nargs='?', default = '')
parser.add_argument('-gff', dest= 'gff', help='The general feature format of the organism of interest, preferably refbank')
parser.add_argument('-genes', dest= 'genes_of_interest', help= 'As complete a list as possible of names of genes known to be regulated by transcription factor, gene names must be in column called Name', nargs='?', default = '')
parser.add_argument('-output', dest= 'outdir', help= "The output directory where created file is stored", nargs='?', default = '')
parser.add_argument('--version', action='version', version='Acronius Tsjits Kramer, version 1, May 2024')

args = parser.parse_args()

#----------------------------------ARGUMENTS-----------------------------------
gff_header = ["chrom","db", "type", "start", "end", "name", "strand", "score", "description"]
gff_descriptions = ['ID', 'Name', 'locus_tag', 'old_locus_tag']	

# -----------------------------DEFINING FUNTIONS-------------------------------
# Function used to write log messages to a file consisting of important messages druring program execution


def gff_add_description_entries(gff):
    def add_description(row, gff_descriptions):
        items = row['description'].split(";")
        for item in items:
            entry = item.split("=")
            if entry and (entry[0] in gff_descriptions):
                row[entry[0]] = entry[1]
        return row
    
    gff = gff.apply(lambda row: add_description(row, gff_descriptions), axis=1)
    return gff

# Function to filter fasta file based on locus tags
def filter_fasta(fasta_file, locus_tags):
    filtered_intergenic = {}
   
    with open(fasta_file, 'r') as f:
        current_header = None
        current_sequence = ''
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header:
#                    print(f"Gene found: {current_header}")
                    if current_header in locus_tags:
                        filtered_intergenic[current_header] = current_sequence             
                current_header = line[1:]  # Remove the '>'
                current_sequence = ''
            else:
                current_sequence += line
        if current_header:
#            print(f"Gene found: {current_header}")
            if current_header in locus_tags:
                filtered_intergenic[current_header] = current_sequence
    return filtered_intergenic

# Function to convert dictionary to fasta format
def convert_to_fasta(intergenic_dict):
    fasta_string = ''
    for locus_tag, sequence in intergenic_dict.items():
        fasta_string += f'>{locus_tag}\n{sequence}\n'
    return fasta_string

# Function to drop rows with NaN values in a specific column
def drop_rows(df, column_name):
    cleaned_df = df.dropna(subset=[column_name])
    return cleaned_df

#----------------------------------MAIN CODE-----------------------------------
##1. Loading in given files and turning them to dataframes
# Loading of csv
genes_df = pd.read_csv(args.genes_of_interest, sep=',', comment='#')


# loading of gff with attributes in description being split into columns of their own
gff_df = pd.read_csv(args.gff, header = None, comment= '#', sep = '\t', names = gff_header)                                  
gff_df = gff_add_description_entries(gff_df)

##2. Filtering .faa file using locus tags of genes in regulon
# Adding locus tag to genes by matching on the column Name
genes_df_merged = pd.merge(genes_df, gff_df[['Name', 'locus_tag']], on='Name', how='left') 
genes_df_merged_dropped = drop_rows(genes_df_merged, 'locus_tag')

# Converging the locus_tags column to a list
locus_tags = ','.join(genes_df_merged_dropped['locus_tag'].tolist())

# Splits the string into a list of locus_tags
locus_tags_list_split = locus_tags.split(',')
print(locus_tags_list_split)

# Cornverts the list into a set of locus_tags so it can be used as input for the filter fasta function
locus_tags_set = set(locus_tags_list_split)
modified_locus_tags_set = {tag.replace('_', '') for tag in locus_tags_set}
print(modified_locus_tags_set)

# Filtering .faa based on locus_tags 
filtered_sequences = filter_fasta(args.proteins, modified_locus_tags_set)

# Turning filtered set back into fasta format
fasta_output = convert_to_fasta(filtered_sequences)

##3. Saving the regulon proteins sequences in a file in the ouput directory
local_output_file = args.outdir+'/'+'filtered_sequences.faa'
with open(local_output_file, 'w') as f:
    f.write(fasta_output)
