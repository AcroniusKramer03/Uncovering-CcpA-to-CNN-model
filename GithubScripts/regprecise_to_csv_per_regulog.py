# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 2024

@author: tsjit
"""

#---------------------------------Modules--------------------------------------

import pandas as pd
import re
import argparse
import os

# ------------------------------FOR CLUSTER USE--------------------------------

parser = argparse.ArgumentParser(description = "Convert regprecise regulog to csv of bacteria and the genes regulated")
parser.add_argument('-regulog', dest='regulog', help='The regulog file in fasta format from the regprecise website without .txt', nargs='?', default = '')
parser.add_argument('-ouput', dest='outdir', help='output folder', nargs='?', default='.')
args = parser.parse_args()

promgramdir = os.path.dirname(os.path.realpath(__file__))
genus = os.path.basename(args.regulog)

# -----------------------------DEFINING FUNTIONS-------------------------------

def drop_rows(df, column_name):
    cleaned_df = df.dropna(subset=[column_name])
    return cleaned_df

#------------------------------------CODE--------------------------------------

# Define file path
file = args.regulog+'.txt'

# Initialize empty lists to store gene and bacteria names
gene_names = []
bacteria_names = []

# Open the text file and read its contents
with open(file) as f:
    file_contents = f.readlines()

    # Regular expressions to match gene names and bacteria names
    gene_pattern = r'\(([^)]+)\)'
    bacteria_pattern = r'\[([^\]]+)\]'

    # Extract gene names and bacteria names using regular expressions
    for line in file_contents:
        gene_match = re.search(gene_pattern, line)
        bacteria_match = re.search(bacteria_pattern, line)

        if gene_match:
            gene_names.append(gene_match.group(1))
        else:
            gene_names.append(None)  # Append None if no gene name is found

        if bacteria_match:
            bacteria_names.append(bacteria_match.group(1))
        else:
            bacteria_names.append(None)  # Append None if no bacteria name is found

# Create a DataFrame
df = pd.DataFrame({
    'Bacteria': bacteria_names,
    'Name': gene_names
})

# Drop the rows with no gene in the Name column
clean_df = drop_rows(df, 'Name')

# transcer the cleaned dataframe of all bacteria to a csv
clean_df.to_csv(args.outdir+'/'+genus+'.csv', index = False)

grouped = clean_df.groupby('Bacteria')

for Bacteria, group_df in grouped:
    file_name = f"{Bacteria}"
    group_df.to_csv(args.outdir+'/'+file_name+'.csv', index=False)