"""Normalize counts from featureCounts

Given the GEO dataset ID, this script loads the corresponding matrix 
that contains the raw expressing counts from featureCounts. Based on 
the gene lengths defined with the help of a biomart metadata table, 
expression counts are normalized by TPM.

Parameters
----------
geo_id : str
	the GEO dataset, either GSE77314 or GSE126848

date:	2019-08-27
author:	Steffen Albrecht
	
"""

from sys import *
import pandas as pd
import numpy as np

geo_id = argv[1]

# load the counts and the biomart gene annotations
data = pd.read_csv('../%s/featureCounts/counts_TP.tsv'%(geo_id), sep='\t')
genes = pd.read_csv('gene_result.txt', sep='\t')

# calculate and collect gene lengths
lens = []
wc = 0
for index, row in genes.iterrows():
	l = np.nan
	try:
		gene_start = int(row['start_position_on_the_genomic_accession'])
		gene_end = int(row['end_position_on_the_genomic_accession'])
		l = int(gene_end - gene_start)
	except:
		wc += 1
	lens.append(l)
print('%d warnings...'%wc)
id_len = dict(zip(genes['GeneID'], lens))

# collect all genes for which the length could be calculated
# (metadata was available)
genes_with_len = []
lens = []
for gene in data['gene_id']:
	l = id_len.get(gene, np.nan)
	if np.isnan(l):
		genes_with_len.append(False)
	else:
		genes_with_len.append(True)
		lens.append(l)
print('%d are missing'%(sum(map(lambda x: not x, genes_with_len))))
data = data.loc[genes_with_len]

# normalize with TPM
for col in data.columns.values:
	if col in set(['gene_id']):
		continue
	by_len = []
	for index, row in data.iterrows():
		count = data[col][index]
		gene = row['gene_id']
		l_kb = float(id_len[gene]) / 1000.0
		by_len.append(float(count) / l_kb)
	mc = float(sum(by_len)) / 1000000.0
	normalized_col = [float(value) / mc for value in by_len]
	print(normalized_col[:10])
	data[col] = normalized_col

# output the normalized values as a matrix and the same matrix transposed
data.to_csv('../%s/featureCounts/tpm_TP.tsv'%(geo_id), 
			sep='\t', index=False)
m = []
for col in data:
	col_name = col
	if col_name == 'gene_id':
		col_name = 'sra_id'
	m.append([col_name] + list(data[col]))
with open('../%s/featureCounts/tpm.tsv'%(geo_id), 'w') as f:
	for row in m:
		f.write('\t'.join(map(str, row)) + '\n')


















