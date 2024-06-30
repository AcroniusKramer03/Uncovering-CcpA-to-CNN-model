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

python3 ProPr_makeCNNmodel_PrepareSeqList.py -ratio 1 -lof examples/ListOfGenomes.txt -padleft 40 -padright 5 -outdir examples

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
parser = argparse.ArgumentParser(description='Train the model')
parser.add_argument('-outdir', dest='outdir', help='output folder', nargs='?', default='.')
parser.add_argument('-lof', dest='lof', help='Table with List of Genomes', default='ListOfGenomes.txt')
parser.add_argument('-tr', dest='training', help='training seqs outputfilename [training_seqs.txt]', default='training_seqs.txt')
parser.add_argument('-bg', dest='background', help='background seqs outputfilename [background_seqs.txt]', default='background_seqs.txt')
parser.add_argument('-fraglen',dest='FragLen', help='Length of the random DNA fragments [51]', nargs='?', default=51)
parser.add_argument('-padleft',dest='padLeft', help='Left padding bases [5]', nargs='?', default=5)  
parser.add_argument('-padright',dest='padRight', help='Right padding bases [5]', nargs='?', default=5)  
parser.add_argument('-forceGC',dest='forceGC', help='Get GC% from DNA [0] or force GC e.g., enter 40', nargs='?', default=0)  
parser.add_argument('-ratio',dest='ratio', help='Ratio background:training seqs [100]', nargs='?', default=100)  
parser.add_argument('--version', action='version', version='Anne de Jong, version 1.0, April 2024')
args = parser.parse_args()

args.FragLen  = int(args.FragLen)
args.padLeft  = int(args.padLeft)
args.padRight = int(args.padRight)
args.ratio = int(args.ratio)
# NOTE, makeManyRandoms is slow; default = 100 which will give for 2000 promoters 100 x 2000= 200.000 non-promoters

gff_header = ["chrom","db", "type", "start", "end", "name", "strand", "score", "description"]
bed_header = ["chrom","start", "end", 'name','score', "strand"]
bases = ['G', 'A', 'T', 'C']

''' --------------------------------------  DEFINE FUNCTIONS ---------------------------------------------------'''

def getCleanSeq(fna_file):
	# get DNA and replace N or n by G or g to prevent error in training; G is considered as most save replacment
	DNA = ''
	with open(fna_file) as lines:
		for line in	lines: 
			if line[0] != ">": 
				DNA += line.strip()
	DNA = re.sub('[BDEFHIJKLMNOPQRSUVWXYZ]','G', DNA.upper())
	return DNA	

def randomseq(dnaset, seqs):
	result = ''
	for i in range(args.FragLen):
		r = random.randint(0, len(dnaset)-1)
		result += dnaset[r]
	if result in seqs: result = randomseq(dnaset, seqs)
	return result

def makeManyRandoms(ACGTcontent, seqs):
	# make 100 bases with proper GC content
	dnaset = "A" * ACGTcontent[0] + "C" *ACGTcontent[1] + "G" * ACGTcontent[2] + "T" * ACGTcontent[3]
	RandomSeq = []
	for i in range(0, args.ratio*len(seqs)): RandomSeq.append(randomseq(dnaset, seqs))
	return RandomSeq

def reverse_complement(dna):
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	return ''.join([complement[base] for base in dna[::-1]])

def write_log(S):
	# write to console and logfile
	print(S)
	f = open(args.outdir+'/'+"ProPr.log", "a")
	f.write(S + '\n')
	f.close()

def Anne_getPercentage(genome):
	# return a list of percentages of A,C,G,T
	percentageList = []
	for base in ['A','C','G','T']:	percentageList.append(round(( genome.count(base) / len(genome) )*100))
	return(percentageList)

def stats_GCcontent(genomes):
	acgt_dict = {}
	for index, genome in genomes.iterrows():
		DNA = getCleanSeq(genome['FilePrefix']+'.fna')
		ACGTcontent = Anne_getPercentage(DNA) #  return a list of percentages of A,C,G,T
		ACGTcontent.append(ACGTcontent[1] + ACGTcontent[2])  # append GC%
		ACGTcontent.append(len(DNA))  # append length
		acgt_dict[genome['Name']] = ACGTcontent
	acgt_df = pd.DataFrame(acgt_dict).T  # Convert to Pandas DataFrame
	acgt_df.columns = ['A', 'C', 'G', 'T', 'GC', 'length']
	acgt_df.to_csv(args.outdir + '/' + 'statistics.csv', index=True, header=True) # Save to csv
  
def adjusted_FramentLength(seq):
	# resize fragment if needed to args.FragLen
	if len(seq) == args.FragLen:
		return seq
	else:
		if len(seq)>args.FragLen:
			print('Too large, fragment will be resized')
			excess_chars = len(seq) - args.FragLen
			trim_size = excess_chars // 2
			return seq[trim_size:-(trim_size if excess_chars % 2 == 0 else trim_size + 1)]
		else:
			N = args.FragLen - len(seq) 
			R = ''.join(random.choice(bases) for _ in range(N))
			return R+seq
			


''' ======================================   MAIN ================================================'''
write_log('Extracting Training Sequences...')

# Get the list of genomes
genomes = pd.read_csv(args.lof, sep='\t', comment="#")
stats_GCcontent(genomes)

training_list = []
background_list = []


def AddPadding(row):
    if row['strand'] == '+':
        row['start'] = int(row['start']) - args.padLeft
        row['end'] = int(row['end']) + args.padRight
    else:
        row['start'] = int(row['start']) - args.padRight
        row['end'] = int(row['end']) + args.padLeft
    return row


for index, genome in genomes.iterrows():
	training_fragments_table = []
	write_log(genome['ID'] + ' ' + genome['Name'])
	DNA = getCleanSeq(genome['FilePrefix']+'.fna')
	bed_df = pd.read_csv(genome['FilePrefix']+'.bed', header=None, sep='\t', comment="#", names=bed_header)
	bed_df = bed_df.dropna()
	write_log('Number of training seqs in ' + genome['Name'] + ' = ' + str(len(bed_df)))
	print(bed_df)
	bed_df = bed_df.apply(AddPadding, axis=1)

	for index, row in bed_df.iterrows():
		seq = DNA[int(row['start']):int(row['end'])]
		if row['strand'] == '-': seq = reverse_complement(seq)
		training_list.append(adjusted_FramentLength(seq))
		training_fragments_table.append(str(row['start']) + '\t' +row['strand']+'\t' +seq) 
	
	write_log('Promoter sequences written to: '+genome['Name'] +'.Training.promoter_seq.txt')
	with open(args.outdir+'/'+genome['Name'] +'.Training.promoter_seq.txt', 'w') as f:
		f.write('\n'.join(training_fragments_table))  
	
	# Create the background sequences based on genome ATGC content 
	if args.forceGC == 0:
		ACGTcontent = Anne_getPercentage(DNA) # percentages of A,C,G,T derived from DNA
	else:									  # GC% from args.forceGC
		ACGTcontent = []
		ACGTcontent.append(round((100-args.forceGC)/2))	# A
		ACGTcontent.append(round(args.forceGC/2))		# C	
		ACGTcontent.append(round(args.forceGC/2))		# G
		ACGTcontent.append(round((100-args.forceGC)/2))	# T

	background_list += makeManyRandoms(ACGTcontent, training_list)


print(training_list[0])	
print(background_list[0])	

# Report and write seqs
with open(args.outdir + '/' + args.training, mode='w') as f:
	f.write('\n'.join(training_list))
	f.write('\n')
with open(args.outdir + '/' + args.background, mode='w') as f:
	f.write('\n'.join(background_list))
	f.write('\n')

write_log('Total number of training seqs:   ' + str(len(training_list)))
write_log('Total number of background seqs: ' + str(len(background_list)))

