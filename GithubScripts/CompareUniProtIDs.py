'''
Author: Anne de Jong 
date: Sept 2023

python3 /data/sigma54/CompareUniProtIDs.py -sessiondir $SESSIONDIR -query query.g2d.diamond.tab -out query.regulon.tab -regulon $PROGRAMDIR/database/Sigma54RegulonProteins.UniProtID 


'''

import argparse
import pandas as pd


# ---------------------------------------------------------------- parse parameters -------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description='Add Regulon to UniProt')
parser.add_argument('-sessiondir', dest='sessiondir', help='Sessionfolder')
parser.add_argument('-query', dest='query', help='diamond result filename of the query: query.g2d.diamond.tab')
parser.add_argument('-out', dest='out', help='output filename : query.regulon.tab')
parser.add_argument('-regulon', dest='regulon', help='Full path and name of the list with UniProt RegulonIDs')
parser.add_argument('--version', action='version', version='Anne de Jong, version 1.0, Sept 2023')
args = parser.parse_args()

print("Compare UniProt IDs to predict putative regulon members");
gff_header = ["chrom","db", "type", "start", "end", "name", "strand", "score", "description"]
diamond_header = ["locus-tag","uniprotID","indentity","length","mismatches","gaps","qStart","qEnd","tStart","tEnd","evalue","bitscore"]


query=pd.read_csv(args.sessiondir+'/'+args.query, header=None, comment='#',sep='\t', names=diamond_header)
regulon=pd.read_csv(args.regulon, header=None, comment='#',sep='\t', names=["uniprotID"])

merged = pd.merge(query, regulon, how='inner', on=['uniprotID'])

filename=args.sessiondir+'/'+args.out
print("Regulon genes/proteins written to "+filename)
merged.sort_values(by=['locus-tag']).to_csv(args.sessiondir+'/'+args.out, index = False, sep ='\t', header=True)

#  load gff  
#gff = pd.read_csv(args.gff, header=None,  comment='#',sep='\t', names=gff_header)
#convert_dict = { 'start': int, 'end': int }
#gff = gff.astype(convert_dict)  # be sure that start and end are integers
#gff_genes  = gff.loc[gff['type'] == 'gene']
#gff_genes.sort_values(by=['start']).to_csv(args.gff, index = False, sep ='\t', columns=gff_header, header=None)



