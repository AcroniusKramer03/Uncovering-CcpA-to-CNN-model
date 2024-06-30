# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 15:30:24 2024

@author: tsjit

python3 "/scratch/hb-molgen/TsIsZiJo/tsjits/scripts/regulon_intergenic.py" -bacteria NAME_OF_BACTERIA -genome PATH/DIRECTLY/TO/EITHER/GFF/OR/FNA (without .gff or fna) -regulon PATH/TO/BACTERIA.regulon.tab" -outdir PATH/TO/OUTPUT/DIRECTORY

"""

#----------------------------------MODULES-------------------------------------
import pandas as pd
import re
import argparse
import os
import subprocess

#------------------------------PARSE PARAMETERS--------------------------------
parser = argparse.ArgumentParser(description='Annotation of Genome')
parser.add_argument('-bacteria', dest='bacteria', help='Name of bacteria, determines the names of output files')
parser.add_argument('-genome', dest='genome', help='fasta and general feature format of bacteria without .fna and .gff', nargs='?', default='')
parser.add_argument('-regulon', dest='regulon', help='The transcription factor regulon discovered by comparing uniprodIDs', nargs='?', default='')
parser.add_argument('-outdir', dest='outdir', help='output directory', nargs='?', default='.')
parser.add_argument('--version', action='version', version='Acronius Tsjits Kramer, version 2.0, 24/5 2024')
args = parser.parse_args()

promgramdir = os.path.dirname(os.path.realpath(__file__))
genomename = args.bacteria

#----------------------------------ARGUMENTS-----------------------------------
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
                if current_header:
                    # Extract the locus tag from the header
                    locus_tag = current_header.split('|')[0]
                    if locus_tag in locus_tags:
                        filtered_intergenic[current_header] = current_sequence
                current_header = line[1:]
                current_sequence = ''
            else:
                current_sequence += line
        # Check the last entry
        if current_header:
            locus_tag = current_header.split('|')[0]
            if locus_tag in locus_tags:
                filtered_intergenic[current_header] = current_sequence
    return filtered_intergenic

def convert_to_fasta(intergenic_dict):
    fasta_string = ''
    for locus_tag, sequence in intergenic_dict.items():
        fasta_string += f'>{locus_tag}\n{sequence}\n'
    return fasta_string


#----------------------------------MAIN CODE-----------------------------------
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

# 6. Selecting intergenic regions
write_log("selecting intergenic regions of interest")
intergenic_of_interest = ''
fnn_2 = ''
prevGene = pd.DataFrame()
first = True
for index, gene in genes.sort_values(by=['chrom','start']).iterrows():
    seq  = get_sequence(FASTA, gene)
    fnn_2 += '>'+str(gene['locus_tag'])+'\n'+str(seq)+'\n'
    if first: 
        first = False
    else:
        if gene['chrom'] == prevGene['chrom']: 
            intergenicSeq = FASTA[gene['chrom']][prevGene['end']:gene['start']]
    
            if len(intergenicSeq) > min_intergenic_len and len(intergenicSeq) < max_intergenic_len:
                if prevGene['strand'] == '-': 
                    intergenic_of_interest += '>' +prevGene['locus_tag']+'|'+str(prevGene['chrom'])+'|'+str(prevGene['end']) +'|'+str(gene['start'])+'|-\n'+reverse_complement(intergenicSeq)+'\n'
                if (gene['strand'] == '+'):
                    intergenic_of_interest += '>'+gene['locus_tag']+'|'+str(gene['chrom'])+'|'+str(prevGene['end'])+'|'+str(gene['start'])+'|+\n'+intergenicSeq+'\n'
    prevGene=gene

f = open(args.outdir+'/'+genomename+'.g2d.fnn', "w")
f.write(fnn_2)
f.close
with open(args.outdir+'/'+genomename+'.g2d.intergenic.ffn', "w") as f:
    f.write(intergenic_of_interest)
    f.flush()
    os.fsync(f.fileno())

# 7. Filtering intergenic regions of interest 
regulon_df = pd.read_csv(args.regulon, comment= '#', sep = '\t')
regulon_locustags = ','.join(regulon_df['locus-tag'].tolist()) # Converts the column "locus_tag" to a list so it can be proccesed by the filter_fasta function

# Take from the ouput directory the file containing intergenic regions and filter de fasta sequence
intergenic_regions =args.outdir+'/'+genomename+'.g2d.intergenic.ffn'
if (os.path.exists(intergenic_regions) is True):
    write_log("filtering for intergenic regions infront of transcription factor regulated genes")
    intergenic_filtered = filter_fasta(intergenic_regions, regulon_locustags)
    intergenic_filtered = convert_to_fasta(intergenic_filtered) # To convert the created library into a fasta format used be meme

local_output_file = args.outdir+'/'+genomename+'.regulon_intergenic.faa'
with open(local_output_file, 'w') as f:
    f.write(intergenic_filtered) 	
    




#------------------------------------------------------------------------------
'''    
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
output_file = os.path.join(args.outdir, genomename+".intergenic.gff")
intergenic_df.to_csv(output_file, sep = '\t', index = False)
'''
