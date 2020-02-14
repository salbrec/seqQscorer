"""One-feature prediction analysis

Within this script, the predictive performance of single features is evaluated
by the area under ROC curve. Later the one-feature performance is compared to
the multi-feature performance achieved by machine learning algorithms 

Methods
-------

print_nice_table()
	function that prints a formated table on the console
get_auROC(values_raw, labels)
	function that calculates the area under ROC curve

date:	2019-08-02
author:	Steffen Albrecht

"""

from sys import *
import random
import os
import numpy as np
import pandas as pd
import pickle
from sklearn.metrics import roc_auc_score
import utils.gs_utils as utils
from terminaltables import AsciiTable

def print_nice_table(table):
	print_table = AsciiTable(table)
	print(print_table.table)

def get_auROC(values_raw, labels):
	auROC = -1.0
	try:
		normed = [float(v)/max(values_raw) for v in values_raw]
		AUCs = []
		AUCs.append(roc_auc_score(labels, normed))
		normed.reverse()
		AUCs.append(roc_auc_score(labels, normed))
		auROC = max(AUCs)
	except Exception:
		print('Error / Warning in claculating the area under ROC curve')
	return auROC

feature_abbr_ls, feature_abbr_sl = utils.get_abbr_maps()

# used to manually select the different settings, settings
# with the same run type can be applied together. Tables were
# manually merged later 
settings = []
#settings.append(('None', 'None', 'None'))
#settings.append(('human', 'RNA-seq', 'single-ended'))
#settings.append(('human', 'ChIP-seq', 'single-ended'))
#settings.append(('mouse', 'ChIP-seq', 'single-ended'))
settings.append(('human', 'DNase-seq', 'paired-ended'))
settings.append(('mouse', 'DNase-seq', 'paired-ended'))

auROCs = {}
cases = []
all_features = []
for setting in settings:
	# define metadata and load the dataset
	organism = setting[0]
	assay = setting[1]
	run_type = setting[2]
	data = pickle.load(open('./datasets/%s_%s_%s'%(assay, organism, run_type),'rb'))
	
	# define the feature set short names depending on the run type
	metrics = ['FQC']
	if run_type == 'None':
		metrics.append('BowM')
	if run_type == 'single-ended':
		metrics.append('BowS')
	if run_type == 'paired-ended':
		metrics.append('BowP')
	metrics.append('RA')
	metrics.append('TSS')
	all_columns = []
	# collect all column names representing the feature names
	for short_name in metrics:
		long_name = feature_abbr_sl[short_name]
		col_names = list(filter(lambda x: x[:len(long_name)] == long_name, 
						  data.columns.values))
		all_columns += col_names
	# fromat a string that describes the surrent setting
	case = '%s %s (%d)'%(organism.title(), assay, data.shape[0])
	cases.append(case)
	auROCs[case] = {}
	# separate feature values and class values 
	X = data[all_columns]
	y = np.array([1 if s == 0 else 0 for s in data['status']])
	
	# for each feature, calculate the predictive performance described
	# by the area under ROC curve
	for feature in all_columns:
		values = np.array(X[feature])
		if not feature in all_features:
			all_features.append(feature)
		auROC = get_auROC(values, y)
		auROCs[case][feature] = auROC

# collect all auROC values within a table
header = ['Feature'] + cases
table = [header]
for feature in all_features:
	temp = [feature]
	for case in cases:
		auROC = auROCs[case].get(feature, 0.0)
		temp.append('%.3f'%(auROC))
	table.append(temp)

# print all auROC values and write them into a text file
print_nice_table(table)
with open('./auROCs_single_features.tsv', 'w') as f:
	for row in table:
		f.write('\t'.join(row).replace('.', ',') + '\n')
























