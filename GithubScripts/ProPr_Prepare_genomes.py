"""
Anne de Jong
2024 April

Simplified version of python3 /data/FACoPv2/03.FACoPv2_annotate_genomes.py

genome=/scratch/hb-molgen/TsIsZiJo/Anne/query_set1/L_lactis_MG1363_ASM942v1_100k
python3 /userapps/hb-molgen/ProPr/ProPr_Prepare_genomes.py -genome $genome -o /scratch/hb-molgen/TsIsZiJo/Anne/results

"""

import sys
import pandas as pd
import argparse 
import subprocess
import os
import glob
import re


parser = argparse.ArgumentParser(description='FACoPv2 Annotation of Genome')
parser.add_argument('-genome', dest='genome', help='Full Genome Basename without .fna and .gff', nargs='?', default='')
parser.add_argument('-o', dest='outdir', help='Output folder', nargs='?', default='.')
parser.add_argument('--version', action='version', version='Anne de Jong, version 2.0, Feb 2024')
args = parser.parse_args()

programdir = os.path.dirname(os.path.realpath(__file__)) # get the root of the package
genomename = os.path.basename(args.genome)

ProTransTerm = 'python3 /userapps/hb-molgen/ProTransTerm/ProTransTerm.py'

gff_header = ["chrom","db", "type", "start", "end", "name", "strand", "score", "description"]
gff_descriptions = ['ID','Name','locus_tag','old_locus_tag']
min_intergenic_len = 15
max_intergenic_len = 600
max_operon_intergenic = 150


##################################################################################################
##                             Functions                                                        ##
##################################################################################################

def run_cmd(cmd):
	print(cmd)
	result = ''
	result = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
	print("Stdout:", result.stdout)
	#print("Stderr:", result.stderr)

def write_log(text):
	f = open(args.outdir+'/sessionprogress', "a")
	f.write(text+'<br>')
	f.close
	print(text)

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

def GFF_add_decription_entries(GFF):
	for description in gff_descriptions: GFF[description]=""  # add entries to GFF
	for index, row in GFF.iterrows():
		items=row['description'].split(";")
		for item in items:
			entry=item.split("=")
			if entry and (entry[0] in gff_descriptions) : # add description entry if entry name exists
				GFF.loc[index,entry[0]] = entry[1]
	return GFF	

def reverse_complement(dna):
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
	#dna = re.sub('[BDEFHIJKLMNOPQRSUVWXYZ]','G', dna)  # only keep GATC all other chars will be a G
	return ''.join([complement[base] for base in dna[::-1]])

def get_sequence(FASTA, gene):
	seq = FASTA[gene['chrom']][gene['start']-1:gene['end']]
	if gene['strand'] == '-': seq = reverse_complement(seq)
	return seq

def translate(seq):
	protein = ""
	for i in range(0, len(seq), 3) : 
		codon = seq[i:i+3]
		if codon in CODONS: protein += CODONS[codon]
	return protein

def load_codons(codon_table_file):
	codons={}
	with open(codon_table_file) as f:
		for line in f:
			items = line.rstrip().split("\t")
			if items[1]: codons[items[0]] = items[1]
	return codons			

def clean_seq(seq):
	unwanted = "'[]/+.!@#$;:!*%)(&^~="
	return ''.join( c for c in seq.strip() if c not in unwanted )

	
##################################################################################################
##                             MAIN                                                             ##
##################################################################################################


CODONS = load_codons(programdir+'/data/codon_table.txt')

write_log('Annotating '+args.genome)
if os.path.getsize(args.genome+'.gff') > 10000000: 
	write_log('Wrong GFF file')
	exit()    # skip large gff files, which should not be in, but NCBI is not perfect


# 1. Read genome sequence data
write_log('\tload FASTA')
FASTA = readMultiFasta(args.genome+'.fna')

if (os.path.exists(args.genome+'.g2d.fna') is False):
	write_log("\twrite FNA to G2D format")
	f = open(args.outdir+'/'+genomename+'.g2d.fna', "w")
	for chrom in FASTA:  # write the old ptt file for each fasta entry . This is needed for the old, but good, TranstermHP
		f.write('>'+chrom+'\n'+FASTA[chrom]+'\n')  # use clean headers
	f.close
	
write_log("\tload GFF")
GFF = pd.read_csv(args.genome+'.gff', header=None,  comment='#',sep='\t', names=gff_header)
convert_dict = { 'start': int, 'end': int }
GFF = GFF.astype(convert_dict)  # be sure that start and end are integers
GFF_add_decription_entries(GFF)


# 2. Write GFF as Table if not exists
filename=args.outdir+'/'+genomename+'.g2d.table'
if (os.path.exists(filename) is False):
	write_log("\tGFF to table")
	header = ["chrom","db", "type", "start", "end", "strand", 'ID','Name','locus_tag','old_locus_tag' ]
	GFF.sort_values(by=['start']).to_csv(filename, index = False, sep ='\t', columns = header)


genes_df = GFF.loc[GFF['type'] == 'gene']
if (genes_df.shape[0] < 2): 
	write_log('Error: no genes are found in the GFF file')
	exit()


# 3. Write genes, intergenic and proteins in G2D FASTA format if not exists
write_log("\tFNN, FAA and intergenic to G2D format")
fnn = ''
faa = ''
intergenic = ''
prevGene = pd.DataFrame()
first = True
for index, gene in genes_df.sort_values(by=['chrom','start']).iterrows():
	seq  = get_sequence(FASTA, gene)
	prot = translate(seq)
	fnn += '>'+gene['locus_tag']+'\n'+seq+'\n'
	faa += '>'+gene['locus_tag']+'\n'+prot+'\n'
	if first: first = False
	else:
		if gene['chrom'] == prevGene['chrom']: 
			intergenicSeq = FASTA[gene['chrom']][prevGene['end']:gene['start']]
			if len(intergenicSeq) > min_intergenic_len and len(intergenicSeq) < max_intergenic_len:
				if prevGene['strand'] == '-': intergenic += '>'+prevGene['locus_tag']+'\n'+reverse_complement(intergenicSeq)+'\n'
				if gene['strand']     == '+': intergenic += '>'+gene['locus_tag']+'\n'+intergenicSeq+'\n'
	prevGene=gene
f = open(args.outdir+'/'+genomename+'.g2d.fnn', "w")
f.write(fnn)
f.close
f = open(args.outdir+'/'+genomename+'.g2d.faa', "w")
f.write(faa)
f.close
with open(args.outdir+'/'+genomename+'.g2d.intergenic.ffn', "w") as f:
    f.write(intergenic)
    f.flush()
    os.fsync(f.fileno()) # Ensure that all data is written to disk and wait until the file is completely written
	

# 4. Transcription Terminator prediction
write_log("\tTranscription Terminator prediction")
TransTermFilename = args.outdir+'/'+genomename+'.g2d.transterm'
cmd = ProTransTerm + ' -i '+ args.outdir+'/'+genomename+'.g2d.intergenic.ffn -o '+TransTermFilename
run_cmd(cmd)


# 5. Operon prediction depends on Transcription Terminator prediction 
if (os.path.exists(TransTermFilename) is False):
	write_log('\tError: Transcription Terminator file needed for Operon prediction')
else:	
	write_log('\tOperon prediction')
	operons = pd.DataFrame()
	TransTerm = pd.read_csv(TransTermFilename, comment='#',sep='\t')
	operon=1
	row = pd.Series(dtype='string')
	first = True
	for index, gene in genes_df.sort_values(by=['chrom','start']).iterrows():
		if first: first = False
		else:
			if   gene['chrom']  != prevGene['chrom']:  operon+=1                   # other chrom				
			elif gene['strand'] != prevGene['strand']: operon+=1                   # gene in other strand
			elif gene['start']-prevGene['end'] > max_operon_intergenic: operon+=1  # large gap
			else:                                                                  # terminator present
				dfSize = TransTerm[TransTerm['Position'].between(prevGene['end'], gene['start'])].size
				if (dfSize>0): operon+=1 
		row['locus_tag']   = gene['locus_tag']
		row['operonID']    = 'operon_'+"{:04d}".format(operon)
		row['description'] = gene['chrom']+';'+gene['description']
		operons = pd.concat([operons, row.to_frame().T], ignore_index=True)
		prevGene = gene	
	operons.sort_values(by=['locus_tag']).to_csv(args.outdir+'/'+genomename+'.g2d.operons', index = False, sep ='\t', header=True)

