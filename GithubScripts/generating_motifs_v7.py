# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 15:30:24 2024

@author: tsjit
"""

# Modules
from Bio import SeqIO
import pandas as pd
import re
import argparse
import os
import glob
import subprocess
import json

#---------------------------------Arguments------------------------------------
parser = argparse.ArgumentParser(description='Annotation of Genome')
parser.add_argument('-genome', dest='genome', help='Reference sequence without .fna and .gff', nargs='?', default='')
parser.add_argument('-output', dest='outdir', help='output folder', nargs='?', default='.')
parser.add_argument('--version', action='version', version='Acronius Tsjits Kramer, version 7, May 2024')
parser.add_argument('-genes', dest='genes_of_interest', help='A comma seperated values file with the names of the genes of the regulon in a column called *Name*', nargs='?', default='')
args = parser.parse_args()

promgramdir = os.path.dirname(os.path.realpath(__file__))
genomename = os.path.basename(args.genome)

ProTransTerm = 'python3 /userapps/hb-molgen/ProTransTerm/ProTransTerm.py'

gff_header = ["chrom","db", "type", "start", "end", "name", "strand", "score", "description"]
gff_descriptions = ['ID', 'Name', 'locus_tag', 'old_locus_tag']
min_intergenic_len = 15
max_intergenic_len = 600
max_operon_intergenic = 150


# -----------------------------DEFINING FUNTIONS-------------------------------


# Function usefull for executing shell commands from python code and capturing the ouput
def run_cmd(cmd):
	print(cmd)
	result = ''
	result = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
	print("Stdout:", result.stdout)
	#print("Stderr:", result.stderr)


# Function used to write log messages to a file consisting of important messages druring program execution
def write_log(text):
	f = open(args.outdir+'/sessionprogress', "a")
	f.write(text+'<br>')
	f.close
	print(text)

# Function for reading and parsing multi-FASTA files in python and sroting sequence identifiers and sequences in a dictionary
def readMultiFasta(fna_file):
	# Note: string adding to dict is slow, here we use an intermediate string; seq
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

#
def GFF_add_description_entries(GFF):
	for description in gff_descriptions: GFF[description]=""  # add entries to GFF
	for index, row in GFF.iterrows():
		items=row['description'].split(";")
		for item in items:
			entry=item.split("=")
			if entry and (entry[0] in gff_descriptions) : # add description entry if entry name exists
				GFF.loc[index,entry[0]] = entry[1]
	return GFF	

# To generate the reverse compliment of the sequence in the fasta file
def reverse_complement(dna):
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	dna = re.sub('[BDEFHIJKLMNOPQRSUVWXYZ]','G', dna)  # only keep GATC
	return ''.join([complement[base] for base in dna[::-1]])

def get_sequence(FASTA, gene):
    try:
        seq = FASTA[gene['chrom']][gene['start']-1:gene['end']]
        if gene['strand'] == '-':
            seq = reverse_complement(seq)
        return seq  
    except Exception as e:
        print("An error occurred:", e) 
        
        return None

def clean_seq(seq):
	unwanted = "'[]/+.!@#$;:!*%)(&^~="
	return ''.join( c for c in seq.strip() if c not in unwanted )

def filter_intergenic_regions(genes_df, intergenic_fasta, output_folder):
    # Read intergenic sequences from FASTA file
    intergenic_sequences = readMultiFasta(intergenic_fasta)
    
    # Filter intergenic sequences corresponding to specific genes
    filtered_intergenic_sequences = {}
    for index, gene in genes_df.iterrows():
        locus_tag = gene['locus_tag']
        if locus_tag in intergenic_sequences:
            filtered_intergenic_sequences[locus_tag] = intergenic_sequences[locus_tag]
    
    # Write filtered intergenic sequences to a new FASTA file
    output_fasta = os.path.join(output_folder, 'filtered_intergenic_regions.fna')
    with open(output_fasta, 'w') as f:
        for locus_tag, sequence in filtered_intergenic_sequences.items():
            f.write(f'>{locus_tag}\n{sequence}\n')
    
    return intergenic_of_interest

def drop_rows(df, column_name):
    cleaned_df = df.dropna(subset=[column_name])
    return cleaned_df

def filter_fasta(fasta_file, locus_tags):
    filtered_intergenic = {}
    
    with open(fasta_file, 'r') as f:
        current_header = None
        current_sequence = ''
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header and current_header in locus_tags:
                    filtered_intergenic[current_header] = current_sequence
                current_header = line[1:]
                current_sequence=''
            else: 
                current_sequence += line
        if current_header and current_header in locus_tags:
           filtered_intergenic[current_header] = current_sequence
    return filtered_intergenic

def convert_to_fasta(intergenic_dict):
    fasta_string = ''
    for locus_tag, sequence in intergenic_dict.items():
        fasta_string += f'>{locus_tag}\n{sequence}\n'
    return fasta_string


#----------------------------------Main code-----------------------------------


# 1. Read genome sequence data
write_log("loading fasta")
FASTA = readMultiFasta(args.genome+'.fna')

if (os.path.exists(args.genome+'.g2d.fna') is False):
	write_log("writing fna to G2D format")
	f = open(args.outdir+'/'+genomename+'.g2d.fna', "w")
	for chrom in FASTA:                                                        
		f.write('>'+chrom+'\n'+FASTA[chrom]+'\n')  # use clean headers
	f.close
    
# 2. Read GFF and convert to df with prefered headers, convert datatypes to integers and add entries to the description column
write_log("loading gff")
GFF_df = pd.read_csv(args.genome+'.gff', header = None, comment= '#', sep = '\t', names = gff_header)
convert_dict = {'start': int, 'end': int}                                      
GFF_df = GFF_df.astype(convert_dict)
GFF_add_description_entries(GFF_df)
        
# 3. Write GFF as Table if not exists
filename=args.outdir+'/'+genomename+'.g2d.table'
if (os.path.exists(filename) is False):
  write_log("writing gff to table")
  header = ["chrom","db", "type", "start", "end", "strand", 'ID','Name','locus_tag','old_locus_tag' ]
  GFF_df.sort_values(by=['chrom']).to_csv(filename, index = False, sep ='\t', columns = header)
  GFF_df.sort_values(by=['start'])
    
# 4. keep only row that have 'gene' listed as type
genes  = GFF_df.loc[GFF_df['type'] == 'gene']
if (genes.shape[0] < 2): 
	write_log('Error: no genes are found in the GFF file')
	exit()
    
# 5. Create an intergenic gff file containing the beginning and end of each intergenic region
write_log("creating intergenic gff")
intergenic = ''
fnn = ''
prevGene = None  # Initialize prevGene as None
intergenic_data = []  # Initialize a list to store intergenic data

for index, gene in genes.sort_values(by=['chrom','start']).iterrows():
    seq = get_sequence(FASTA, gene)  # Get sequence for the current gene
    fnn += '>' + str(gene['locus_tag']) + '\n' + str(seq) + '\n'  # Append sequence to fnn

    if prevGene is not None and gene['chrom'] == prevGene['chrom']:
        intergenicSeq = FASTA[gene['chrom']][prevGene['end']:gene['start']]
        if min_intergenic_len < len(intergenicSeq) < max_intergenic_len:
            if prevGene['strand'] == '-':
                intergenic += '>' + prevGene['locus_tag'] + '\n' + reverse_complement(intergenicSeq) + '\n'
            if gene['strand'] == '+':
                intergenic += '>' + gene['locus_tag'] + '\n' + intergenicSeq + '\n'

            # Store intergenic data
            intergenic_data.append({
                'chrom': gene['chrom'],
                'locus_tag': gene['locus_tag'],
                'start': prevGene['end'],
                'end': gene['start'],
                'strand': gene['strand']
            })

    prevGene = gene

# Create DataFrame from intergenic data
intergenic_df = pd.DataFrame(intergenic_data)
output_file = os.path.join(args.outdir, "intergenic.gff")
intergenic_df.to_csv(output_file, sep = '\t', index = False)


# 6. Selecting intergenic regions of interest
write_log("selecting intergenic regions of interest")
intergenic_of_interest = ''
fnn_2 = ''
prevGene = pd.DataFrame()
first = True
for index, gene in genes.sort_values(by=['chrom','start']).iterrows():
	seq  = get_sequence(FASTA, gene)
	fnn_2 += '>'+str(gene['locus_tag'])+'\n'+str(seq)+'\n'
	if first: first = False
	else:
		if gene['chrom'] == prevGene['chrom']: 
			intergenicSeq = FASTA[gene['chrom']][prevGene['end']:gene['start']]
			if len(intergenicSeq) > min_intergenic_len and len(intergenicSeq) < max_intergenic_len:
				if prevGene['strand'] == '-': intergenic_of_interest += '>'+prevGene['locus_tag']+'\n'+reverse_complement(intergenicSeq)+'\n'
				if gene['strand']     == '+': intergenic_of_interest += '>'+gene['locus_tag']+'\n'+intergenicSeq+'\n'
	prevGene=gene
f = open(args.outdir+'/'+genomename+'.g2d.fnn', "w")
f.write(fnn_2)
f.close
with open(args.outdir+'/'+genomename+'.g2d.intergenic.ffn', "w") as f:
    f.write(intergenic_of_interest)
    f.flush()
    os.fsync(f.fileno())

genes_of_interest_df = pd.read_csv(args.genes_of_interest, comment= '#', sep = ',')
genes_of_interest_df = pd.merge(genes_of_interest_df, genes[['Name', 'locus_tag', 'start', 'end']], on='Name', how='left') # Creates a dataframe with the genes of interest and their locus_tags, joines them on the "Name" column
genes_of_interest_df = drop_rows(genes_of_interest_df, 'locus_tag') # Discards genes that don't have a locus tag
locus_tags_of_interest = ','.join(genes_of_interest_df['locus_tag'].tolist()) # Converts the column "locus_tag" to a list so it can be proccesed by the filter_fasta function

# Take from the ouput directory the file containing intergenic regions and filter de fasta sequence
intergenic_regions =args.outdir+'/'+genomename+'.g2d.intergenic.ffn'
if (os.path.exists(intergenic_regions) is True):
    write_log("filtering for intergenic regions infront of transcription factor regulated genes")
    intergenic_filtered = filter_fasta(intergenic_regions, locus_tags_of_interest)
    intergenic_filtered = convert_to_fasta(intergenic_filtered) # To convert the created library into a fasta format used be meme

with open(args.outdir+'/'+genomename+'.g2d.intergenic_filtered.ffn', "w") as f: # Store the ouput in the ouput directery given in the parser "-o"
    f.write(intergenic_filtered)
    f.flush()
    os.fsync(f.fileno()) 	
    




#------------------------------------------------------------------------------


'''
def filter_fasta(fasta_file, locus_tags):
    filtered_intergenic = {}
    
    with open(fasta_file, 'r') as f:
        current_header = None
        current_sequence = ''
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header and current_header in locus_tags:
                    filtered_reads[current_header] = current_sequenc
                current_header = line[1:]
                current_sequence=''
            else: 
                current_sequence += line
        if current_header and current_header in locus_tags:
           filtered_reads[current_header] = current_sequence
    return filtered_reads

for index, row in df_fna.iterrows():
    sequence = row['sequence']
    genic = sequence[start:end]  # Extract the letters 'rop'
    print(genic)

def
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id in sequence_ids:
            sequences[record.id] = str(record.seq)
    return sequences

df_fna.read = True

extracted_sequences = extract_sequences(df_fna, df_gff)

def fna_to_dataframe(fasta):
    # Initialize an empty list to store sequences
    sequences = []

    # Read the FASTA file and generate sequences
    fasta_sequence = SeqIO.parse(open(fasta_file), 'fasta')
    for fasta in fasta_sequence: 
        sequences.append([fasta.id, str(fasta.seq)])
    
    # Create a DataFrame from the list of sequences
    df_fasta = pd.DataFrame(sequences, columns=['identifier', 'sequence'])
    
    return df_fasta

def parse_gff_line(line):
    fields = line.strip().split('\t')
    seqid, source, feature, start, end, score, strand, phase, attributes = fields
    return {
        'chrom': seqid,
        'db': source,
        'type': feature,
        'start': int(start),
        'end': int(end),
        'name': score,
        'strand': strand,
        'score': phase,
        'description': attributes
        }

def gff_to_dataframe(gff_file):
    with open(gff_file, 'r') as f:
        lines = f.readlines()
        
    gff_data = [parse_gff_line(line) for line in lines if not line.startswith('#')]
    df = pd.DataFrame(gff_data, columns=gff_header)
    return df

# Function to extract locus_tag
def extract_locus_tag(attribute_str):
    match = re.match(r'.*locus_tag=([^;]+).*', attribute_str)
    if match:
        return match.group(1)
    else:
        return None
    
# Function to extract Name
def extract_name(attribute_str):
    match = re.match(r'.*Name=([^;]+).*', attribute_str)
    if match:
        return match.group(1)
    else:
        return None   

# Function to drop rows containing missing values in a certain column
def drop_rows(df, column_name):
    cleaned_df = df.dropna(subset=[column_name])
    return cleaned_df

'''
