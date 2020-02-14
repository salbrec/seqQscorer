"""Plot PCA for a given GEO dataset

Based on the tpm expression values of the given GEO dataset this script
post-processes the expression values and applies the PCA. The first two
principal components are used to plot samples in two dimensions. As
described in the paper, the quality probabilities are added to 10% of
the samples with the highest probabilities.

Parameters
----------
geo_id : str
	the GEO dataset, either GSE77314 or GSE126848

date:	2019-08-30
author:	Steffen Albrecht
	
"""

from sys import *
import os
import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import preprocessing
from sklearn.decomposition import PCA

# get GEO ID and create file paths
geo_id = argv[1]
tpm_data_file = './%s/featureCounts/tpm.tsv'%(geo_id)
sra_data_file = './%s/SRA_list_comp.tsv'%(geo_id)
quality_probs = './%s/quality_prbabilities.tsv'%(geo_id)

# read in the data needed: TPM values, SRA metadata, quality probabilities
tpm_data = pd.read_csv(tpm_data_file, sep='\t')
sra_data = pd.read_csv(sra_data_file, sep='\t')
quality_data = pd.read_csv(quality_probs, sep='\t', names=['sra', 'probability', 'descr'])

# extract column names that describe a gene ID
sample_col = tpm_data.columns[0]
sample_ids = list(tpm_data[sample_col])
gene_cols = list(filter(lambda x: not x == sample_col, tpm_data.columns.values))

# separate expression values, log transform, and apply scaling (standardize)
X = tpm_data[gene_cols].values
for i in range(X.shape[1]):
	X[:,i] = list(map(lambda x: math.log(x+1, 2), X[:,i]))
X = preprocessing.scale(X, axis=0)

# apply default PCA on the expression values
pca = PCA(n_components=2)
dim_red_results = pca.fit_transform(X)
tpm_data['comp1'] = dim_red_results[:,0]
tpm_data['comp2'] = dim_red_results[:,1]

# add sample information from the SRA metadata
hue_map = dict(zip(sra_data[sample_col], sra_data['disease']))
tpm_data['Disease'] = [hue_map[sra] for sra in tpm_data[sample_col]]
if geo_id == 'GSE126848':
	style_map = dict(zip(sra_data[sample_col], sra_data['gender']))
	tpm_data['Gender'] = [style_map[sra] for sra in tpm_data[sample_col]]

# add quality probabilities to the dataset
q_map = dict(zip(quality_data['sra'], quality_data['probability']))
tpm_data['Quality'] = [q_map[sra] for sra in tpm_data[sample_col]]

# define values that format the plot
fs_rank = 14.0
rank_diff = 4.0
star_x_shift = 5
star_y_shift = 0
fs_legend_title = 14
fs_ticks = 10
fs_axis_labels = 15

# plot quality probabilities of 10% with highest values
threshold = list(sorted(tpm_data['Quality']))[-int(0.1 * tpm_data.shape[0])]
if True:
	for index, row in tpm_data.iterrows():
		if row['Quality'] >= threshold:
			x = row['comp1']
			y = row['comp2']
			plt.text(x=x+5, y=y, s = '%.2f'%(row['Quality']), size = 12)

# create and format the PCA plot
ax = None
if geo_id == 'GSE126848':
	ax = sns.scatterplot(x='comp1',y='comp2',hue='Disease', style='Gender', data=tpm_data, **{'s':75})
else:
	ax = sns.scatterplot(x='comp1',y='comp2',hue='Disease', data=tpm_data, **{'s':75})

plt.legend(loc='upper left', title_fontsize=16, fontsize=14,
			 framealpha=0.5, frameon=True)
plt.xticks(fontsize=fs_ticks)
plt.yticks(fontsize=fs_ticks)
plt.xlabel('Principal Component 1', fontsize=fs_axis_labels)
plt.ylabel('Principal Component 2', fontsize=fs_axis_labels)
title = 'PCA on GEO dataset'
ax.set_title(title, fontsize=16, fontstyle='normal')
fig = ax.get_figure()
fig.set_size_inches(6, 6)
fig.savefig('./PCA.svg')










