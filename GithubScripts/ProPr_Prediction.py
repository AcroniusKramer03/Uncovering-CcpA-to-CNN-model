'''
Author: Anne de Jong 
date: March 2024
PyVersion: 3.7.6

Predict promoters in any DNA using best model
1) Best model selected on the basis of GC%
Update March 2024, Model made with new Keras library. Models are renames from .h5 to .keras

module load TensorFlow/2.11.0-foss-2022a

query=/scratch/hb-molgen/TsIsZiJo/Anne/query_set1/L_lactis_MG1363_ASM942v1_100k
python3 /userapps/hb-molgen/ProPr/ProPr_Prediction.py -sessiondir /scratch/hb-molgen/TsIsZiJo/Anne -query $query -modelName FGC3040_36

'''


import os
import sys
import re
import argparse
import pandas as pd
import numpy as np
from statistics import median_high
from tensorflow.keras.models import load_model
from scipy.cluster.hierarchy import ward, fcluster
from scipy.cluster.hierarchy import fclusterdata
from scipy.spatial.distance import pdist
from sklearn.metrics import confusion_matrix


# ---------------------------------------------------------------- parse parameters -------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description='Promoter Prediction')
parser.add_argument('-sessiondir', dest='sessiondir', help='Session dir, aka output dir', nargs='?', default='.')
parser.add_argument('-query', dest='query', help='Full path and name genome basename;  without .fna or .gff extension')
parser.add_argument('-modelName',dest='modelName', help='Full path and filename of the model', nargs='?', default='FGC3545_41')  # NEW
parser.add_argument('-promlen',dest='PromLen', help='Length of Promoter', nargs='?', default=71)
parser.add_argument('-prime',dest='prime', help='5 prime region after TSS', nargs='?', default=0)
parser.add_argument('-pval',dest='pvalue', help='p-value cutoff for initial prediction', nargs='?', default=0.99)
parser.add_argument('-out',dest='outPrefix', help='Prefix for Output files', nargs='?', default='ProPr')
parser.add_argument('--version', action='version', version='Anne de Jong, version 2.1, March 2024')
args = parser.parse_args()

#  ---- be sure Arguments are handled as int or real value ------
args.PromLen = int(args.PromLen)
args.prime   = int(args.prime)
args.pvalue  = float(args.pvalue)

# ----- For clustering -----------------------------
MIN_CLUSTER_SIZE = 1
WINDOW_SIZE = 10
PROBABILITY_CUTOFF = 0.99

# ------ For analysis ------------------------------
SHIFT_FOR_NEXT_WINDOW = 1

#  gff header for genome annotation 
gff_header = ["chrom","db", "type", "start", "end", "name", "strand", "score", "description"]


''' ==============================  DEFINING FUNCTIONS =============================================== '''

def Anne_one_hot_encode(seq):
    mapping = dict(zip("ACGT", range(4)))
    seq2 = [mapping[i] for i in seq]
    return np.eye(4)[seq2]

def remPrevFormatting(fna_file):
    f = open(fna_file)
    output = ''
    for line in f:
        if line[0] != ">":
            output = output + line.strip()
    for i in range(0, len(output)):
        if output[i] == 'N' or output[i] == 'n':
                output = output.replace(output[i], 'G')
    return output

def getCleanSeq(fna_file):
	# get DNA and replace N or n by G or g to prevent error in training; G is considered as most save replacment
	DNA = ''
	with open(fna_file) as lines:
		for line in	lines: 
			if line[0] != ">": 
				DNA += line.strip()
	DNA = re.sub('[BDEFHIJKLMNOPQRSUVWXYZ]','G', DNA.upper())
	return DNA	

def getFastaKey(fna_file):
	# result: the header between > and the first space
	fline=open(fna_file).readline().rstrip()
	keys = fline.split()
	return keys[0][1:]	

def makeQuerySet(DNA, window_shift):
	query_set = []
	for i in range(0, len(DNA)-args.PromLen-args.prime, window_shift):
		query_set.append(DNA[i:i+args.PromLen+args.prime])
	return query_set

def GCpercent(DNA):
	DNA = DNA.upper()
	return 100 * (DNA.count('G') + DNA.count('C')) / len(DNA)

def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])

def write_log(S):
	# write to console and logfile
	print(S)
	f = open(args.sessiondir+'/00.ppp.log', "a")
	f.write(S + '\n')
	f.close()	

def adjust_TSS():
	# Anne; adjusted_TSS added: Due to clustering the TSS shifts a few bases. here a correction is made on the basis of -10 position: TGxTATAAT
	# Expected TSS position = len(seq)-16
	# derived from test module BestTATAAT.py
	for index, row in promoter_db_centre.iterrows():
		seq= row['sequence']
		maxscore = 0
		bestpos = 0
		startpos = len(seq)-22  # TSS - 22; Search for -10 from this position to the end of the string
		for i in range(startpos, len(seq)-9):  # screen area for TGxTATAAT from len sequence -22 to end
			score=0 ;
			TGNTATAAT = seq[i:(i+9)]
			score -= (2 * TGNTATAAT.count('C'))
			score -= (2 * TGNTATAAT.count('G'))
			if (TGNTATAAT[0:2]=='TG'): score += 4
			if (TGNTATAAT[3:9]=='TATAAT'): score += 10
			for j in range(3, 8):
				duo=TGNTATAAT[j:(j+2)]
				if   (duo=='AA'): score += 3
				elif (duo=='TA'): score += 3
				elif (duo=='AT'): score += 3
				elif (duo=='TT'): score += 2
			if (score>maxscore):
				maxscore = score
				bestpos = i
		promoter_db_centre.loc[index, 'min10seq'] = seq[bestpos:(bestpos+9)]
		promoter_db_centre.loc[index, 'min10score'] = str(maxscore) ;

		# shift the TSS and sequence of the promoter
		TSSshift = 0
		if (bestpos > startpos):  TSSshift = int(16 - (len(seq) - bestpos))  # TGxTATAAT is TSS-16
		if (promoter_db_centre.loc[index, 'strand'] == '+'): 
			adjTSS = row.TSS + TSSshift ;
			PromSeq = DNA_sense[(adjTSS - args.PromLen):adjTSS]; 
		else:  
			adjTSS = row.TSS - TSSshift ;  
			PromSeq = reverse_complement(DNA_sense[adjTSS:(adjTSS + args.PromLen)]); 
		promoter_db_centre.loc[index, 'adjTSS'] = adjTSS ;
		promoter_db_centre.loc[index, 'sequence'] = PromSeq ;
		print(PromSeq+ ' ' +seq[bestpos:(bestpos+9)]+' pos='+str(bestpos)+' score='+str(maxscore)+' strand='+row.strand)
	

''' ==============================  LOAD MODEL ========================================================== '''

programdir = os.getcwd() 
programdir = os.path.dirname(os.path.realpath(__file__))
cnn_lstm = programdir+'/models/'+args.modelName+'/cnn_lstm.keras'
write_log('==>  LOAD MODEL ') 
write_log('  model file ='+cnn_lstm)
model = load_model(cnn_lstm)
write_log(' ')


''' ==============================  ENCODE DNA SEQUENCE ========================================================== '''
write_log('==> ENCODE DNA SEQUENCE    ')

DNA_sense = getCleanSeq(args.query+'.fna')
DNA_antisense = reverse_complement(DNA_sense)
write_log('  Number of bases in Clean DNA sequence: '+str(len(DNA_sense)))

test_sense_sequences     = makeQuerySet(DNA_sense,  1)
test_antisense_sequences = makeQuerySet(DNA_antisense,  1)

write_log("  Encode sense strand")
input_sense_features = []
for sequence in test_sense_sequences:
	input_sense_features.append(Anne_one_hot_encode(sequence))

write_log("  Encode anti-sense strand")
input_antisense_features = []
for sequence in test_antisense_sequences:
	input_antisense_features.append(Anne_one_hot_encode(sequence))

write_log('  Stacking features')
np.set_printoptions(threshold=40)
sense_features     = np.stack(input_sense_features)
antisense_features = np.stack(input_antisense_features)

write_log('  Sense and anti-sense-feature ready for prediction')
write_log(' ')


''' ==============================  MAKING PREDICTIONS ====================================================== '''
write_log('==> MAKING PREDICTIONS  ')

#Create pandas DataFrame as container for promoter data
promoter_db = pd.DataFrame(columns=['position','score','strand','sequence'])
promoter_db.set_index('position')

write_log('  ===> Prediction on sense strand')
predicted_sense_labels = model.predict(sense_features)  # make a numpy.ndarray; original label and predicted label

# save the first 1000 for evaluation
#np.savetxt(args.sessiondir+'/'+"predicted_sense_labels.csv", predicted_sense_labels[0:1000], delimiter="\t")

# Get promoters from sense strand with prediction pvalue > str(args.pvalue)
predicted_sense_promoter_list = []
probabilityValueSense = []
for i in range(0,len(predicted_sense_labels)):
	if (predicted_sense_labels[i][1]) > args.pvalue:
		probabilityValueSense.append(str(predicted_sense_labels[i][1]))
		predicted_sense_promoter_list.append(test_sense_sequences[i])  # Get the DNA sequence
		new_row = {'position': i, 'score': predicted_sense_labels[i][1], 'strand': '+', 'sequence': test_sense_sequences[i]}
		promoter_db = pd.concat([promoter_db, pd.DataFrame([new_row])], ignore_index=True)

write_log('    Number of putative promoters sense strand      : ' + str(len(predicted_sense_promoter_list)))

write_log('  ===> anti-sense strand')
predicted_antisense_labels = model.predict(antisense_features)  # make a numpy.ndarray; original label and predicted label

# Get promoters from anti-sense strand with prediction pvalue > str(args.pvalue)
predicted_antisense_promoter_list = []
probabilityValueAntisense = []
for i in range(0,len(predicted_antisense_labels)):
	if (predicted_antisense_labels[i][1]) > args.pvalue:
		probabilityValueAntisense.append(str(predicted_antisense_labels[i][1]))
		predicted_antisense_promoter_list.append(test_antisense_sequences[i])   # Get the DNA sequence
		new_row = {'position': len(DNA_sense) - i, 'score': predicted_antisense_labels[i][1], 'strand': '-', 'sequence': test_antisense_sequences[i]}
		promoter_db = pd.concat([promoter_db, pd.DataFrame([new_row])], ignore_index=True)

write_log('    Number of putative promoters anti-sense strand : ' + str(len(predicted_antisense_promoter_list)))
write_log('    Total number of putative predicted promoters: ' + str(len(promoter_db)))

if len(promoter_db) == 0: # write empty result file
	write_log('No promoters found')
	pd.DataFrame(columns=["No promoters found"]).to_csv(args.sessiondir+'/'+args.outPrefix+'.Promoters.txt')
	os._exit(1)		   

if len(promoter_db) >50000: # write empty result file
	write_log('Error; Too many promoters found. Exceeding limit of 50,000. Data size= '+str(len(promoter_db)))
	pd.DataFrame(columns=["Error; Too many promoters found"]).to_csv(args.sessiondir+'/'+args.outPrefix+'.Promoters.txt')
	os._exit(1)	

write_log(' ')


''' ===============================  PERFORM HIERARCHICAL CLUSTERING ========================================= '''
write_log(' ==> PERFORM HIERARCHICAL CLUSTERING  ')

# Clustering written by Daniel Kaptijn, modified by Anne
# For debugging: Optional reload prediction results DataFrame promoter_db
# For debugging: promoter_db_header=['position','score','strand','sequence']
# For debugging: promoter_db = pd.read_csv(args.sessiondir+'/'+args.outPrefix+'.Promoters.initial.txt',  sep ='\t' ,header=0,  names=promoter_db_header)

if len(predicted_sense_promoter_list)+len(predicted_antisense_promoter_list) == 0:
	with open(args.sessiondir+'/'+args.outPrefix+'.Promoters.txt', "w") as f:
		f.write("No promoters found")
	sys.exit()	
	
promoter_db['centre'] = np.NaN

if len(predicted_sense_promoter_list) > 1:
	write_log('  Cluster Sense Strand')
	Xs = promoter_db[['position','strand']].copy()	
	Xs = Xs[Xs['strand']=='+']
	Xs = Xs.drop(labels='strand', axis=1)
	Xs['2D'] = 0
	sense_pred = fclusterdata(Xs, t=WINDOW_SIZE, criterion='distance')
	sense_dict = {}
	for i in range(0, len(sense_pred)):
		key = int(sense_pred[i])       # cluster ID
		value = int(Xs['position'][i])  # TSS position
		if key not in sense_dict.keys():
			sense_dict[key] = [value]
		else:
			new_value = [i for i in sense_dict[key]]
			new_value.append(value)
			sense_dict[key] = new_value
	for i in sense_dict.keys():
		if len(sense_dict[i]) >= MIN_CLUSTER_SIZE:
			cluster_centre = median_high(sense_dict[i])
			promoter_db.loc[(promoter_db['position'] == cluster_centre) & (promoter_db['strand'] == '+'), 'centre'] = cluster_centre
		else:
			promoter_db.loc[(promoter_db['position'] == sense_dict[i]) & (promoter_db['strand'] == '+'), 'centre'] = sense_dict[i]

if len(predicted_antisense_promoter_list) > 1:
	write_log('  Cluster Anti-Sense Strand')
	Xs = promoter_db[['position','strand']].copy()
	Xs = Xs[Xs['strand']=='-']
	Xs = Xs.drop(labels='strand', axis=1)
	Xs = Xs.iloc[::-1].reset_index(drop=True)		# reverse order and reindex for anti-sense strand
	Xs['2D'] = 0
	antisense_pred = fclusterdata(Xs, t=WINDOW_SIZE, criterion='distance')
	antisense_dict = {}
	for i in range(0, len(antisense_pred)):
		key = int(antisense_pred[i])
		value = int(Xs['position'][i])
		if key not in antisense_dict.keys():
			antisense_dict[key] = [value]
		else:
			new_value = [i for i in antisense_dict[key]]
			new_value.append(value)
			antisense_dict[key] = new_value

	for i in antisense_dict.keys():
		if len(antisense_dict[i]) >= MIN_CLUSTER_SIZE:
			cluster_centre = median_high(antisense_dict[i])
			promoter_db.loc[(promoter_db['position'] == cluster_centre) & (promoter_db['strand'] == '-'), 'centre'] = cluster_centre
		else:
			promoter_db.loc[(promoter_db['position'] == antisense_dict[i]) & (promoter_db['strand'] == '-'), 'centre'] = antisense_dict[i]

write_log('    ===========================================================')
write_log('    Initial number of predictions:           ' +str(len(promoter_db)))
write_log('    Number of predictions after clustering:  ' +str(promoter_db.centre.count()))
write_log('    Saving clustered data: '+ args.outPrefix+'.Promoter.db.clustered.txt')
write_log('    ===========================================================')
promoter_db.to_csv(args.sessiondir+'/'+args.outPrefix+'.Promoter.db.clustered.txt', index = False, sep ='\t')


''' ==============================  CENTRE PROMOTERS ========================================= '''
# make a new df with centre promoters only
promoter_db['score'] = promoter_db['score'].apply(lambda x: round(x*10-9, 2))
promoter_db_centre = promoter_db
promoter_db_centre.dropna(inplace=True)  # remove all rows with NaN => these are the non-centre promoters
del promoter_db_centre['centre']
write_log(' ')


''' ==============================  ADD TSS POSITION ========================================= '''
write_log('==> ADD TSS POSITION  ')

# TSS added: The position is the start of the promoter. TSS is position + length of promoter
promoter_db_centre.loc[promoter_db_centre['strand'] == '+', 'TSS'] = promoter_db_centre['position'] + args.PromLen 
promoter_db_centre.loc[promoter_db_centre['strand'] == '-', 'TSS'] = promoter_db_centre['position'] - args.PromLen 
promoter_db_centre.sort_values(by="TSS", ascending=True, inplace=True)
promoter_db_centre.reset_index(drop=True, inplace=True)


''' ==============================  ADD Corrected TSS POSITION ========================================= '''
write_log('  Adjust TSS POSITION and -10 score ')
adjust_TSS()
write_log(' ')


''' ==============================  CLASSIFY PROMOTERS ========================================= '''
write_log('==> CLASSIFY PROMOTERS  ')

# Anne, Classification added: Use the original features annotation (from gff) to map and classify TSS position
genes_gff = pd.read_csv(args.query+".gff", header=None,  comment='#',sep='\t', names=gff_header)
convert_dict = { 'start': int, 'end': int }
genes_gff = genes_gff.astype(convert_dict)  # be sure that start and end are integers

# Here we use the start and end of genes from the GFF
promoter_db_centre['class'] = 'intergenic'
for index, feature in genes_gff.loc[genes_gff['type'] == 'gene'].iterrows():
	promoter_db_centre.loc[(promoter_db_centre['TSS'] > feature.start) & (promoter_db_centre['TSS'] < feature.end), 'class'] = 'feature'

write_log('  Number of rows with class intergenic: '+ str(len(promoter_db_centre[promoter_db_centre['class'] == 'intergenic'])))
write_log('  Number of rows with class feature: '+    str(len(promoter_db_centre[promoter_db_centre['class'] == 'feature'])))
write_log(' ')


''' ==============================  EXPORT PROMOTER DATA ========================================= '''
write_log('==> EXPORT PROMOTER DATA  ')
FastaKey = getFastaKey(args.query+'.fna')

write_log('  EXPORT promoters as GFF ')
def Export_GFFrows(index, row):
	line = FastaKey + '\tProPr\tpromoter' 
	if row['strand'] == '+':
		line += '\t' + str(row['position']) + '\t' + str(row['TSS'])
	else:
		line += '\t' + str(row['TSS']) + '\t' + str(row['position'])
	Name = 'prom_' + format(index, '04d')
	line += '\t'+Name+'\t'+row['strand']+'\t'+str(row['score'])
	min10class = 'strong' if int(row['min10score']) >= 5 else 'weak'
	line += '\tID='+Name+';locus_tag='+Name+';TSS='+str(row['TSS'])+';class='+row['class']+';min10seq='+row['min10seq']+';min10score='+row['min10score']+';min10class='+min10class+';promseq='+row['sequence']
	return line + '\n'

lines = '#'+"\t".join(gff_header) + '\n'
lines += promoter_db_centre.sort_values(by=['TSS']).apply(lambda x: Export_GFFrows(x.name, x), axis=1).str.cat(sep='')
with open(args.sessiondir+'/'+args.outPrefix+'.Promoters.gff', 'w') as f: f.write(lines)

write_log('  EXPORT promoters as TABLE')	
def Export_promoterRows(index, row):
	my_list = []
	my_list.append(FastaKey)
	my_list.append('prom_' + format(index, '04d'))
	if row['strand'] == '+':
		my_list.append(row['position'])
		my_list.append(row['TSS'])
	else:
		my_list.append(row['TSS'])
		my_list.append(row['position'])
	my_list.append(row['adjTSS'])	
	my_list.append(row['strand'])
	my_list.append(row['score'])
	my_list.append(row['min10seq'])
	my_list.append(row['min10score'])
	my_list.append(row['class'])
	my_list.append(row['sequence'])
	return '\t'.join(str(item) for item in my_list) + '\n'

table_header = ["chrom", "ID", "start", "end", "adjTSS", "strand", "score", "min10seq", "min10score",  "class", "sequence"]
lines = "\t".join(table_header) + '\n'
lines += promoter_db_centre.sort_values(by=['TSS']).apply(lambda x: Export_promoterRows(x.name, x), axis=1).str.cat(sep='')
with open(args.sessiondir+'/'+args.outPrefix+'.Promoters.txt', 'w') as f: f.write(lines)

