"""Generic model performance

Within this script, the predictive performance of the generic model is
evaluated when it is applied to subsets specific the different
assay-species specific subsets.

Methods
-------

print_nice_table()
	function that prints a formated table on the console

date:	2019-06-10
author:	Steffen Albrecht
	
"""


from sys import *
import random
import requests, json
import os
import importlib.util
import numpy as np
import pandas as pd
import warnings
import datetime
import pickle
from sklearn.metrics import roc_auc_score
from sklearn.utils import shuffle
import utils.gs_utils as utils
from terminaltables import AsciiTable

with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	from sklearn.model_selection import StratifiedKFold
	from sklearn.model_selection import GridSearchCV
	from sklearn.model_selection import cross_validate
	from sklearn.impute import SimpleImputer
	from sklearn.model_selection import cross_val_predict

def print_nice_table(table):
	print_table = AsciiTable(table)
	print(print_table.table)

# define variables for randomization and the feature set name maps
global the_RS, SEED
SEED = 1
the_RS = np.random.RandomState(SEED)
feature_abbr_ls, feature_abbr_sl = utils.get_abbr_maps()

# load the dataset
data = pd.read_csv('./datasets/clean_dataset_2072.csv', sep='\t')
data = shuffle(data, random_state=the_RS)
data = shuffle(data, random_state=the_RS)
feature_sets = ['FQC', 'BowM', 'RA', 'TSS', 'FQC-BowM-RA-TSS']

# create metadata lists used to associate single samples
data_accession = list(data['accession'])
data_org = list(data['organism'])
data_assay = list(data['assay'])
data_status = list(data['status'])
data_rt = list(data['runType'])

# get information about the dataset
revoked = len(list(filter(lambda x: x == 0, list(data['status']))))
released = len(list(filter(lambda x: x == 1, list(data['status']))))
n_samples = data.shape[0]
rev_perc = float(revoked) / float(n_samples) * 100.0

print('Dataset has %d samples in total'%(n_samples))
print('\t %d are revoked (%.1f%%)'%(revoked, rev_perc))
print('\t %d are released'%(released))
print('')

# define all different combinations based on the metadata
# some combinations are skipped because of the limited data avilable
orgas = ['human', 'mouse']
assays = ['ChIP-seq', 'DNase-seq', 'RNA-seq']
run_types = ['single-ended', 'paired-ended']
combis = ['Various_Various_Various']
combis += [orga+'_Various_Various' for orga in orgas]
combis += ['Various_'+assay+'_Various' for assay in assays]
combis += ['Various_Various_'+rt for rt in run_types]
for orga in orgas:
	for assay in assays:
		if orga == 'mouse' and assay == 'RNA-seq':
			continue
		combis.append('%s_%s_Various'%(orga, assay))
for orga in orgas:
	for assay in assays:
		for rt in run_types:
			if orga == 'mouse' and assay == 'RNA-seq':
				continue
			if orga == 'human' and assay == 'RNA-seq' and rt == 'paired-ended':
				continue
			if orga == 'mouse' and assay == 'ChIP-seq' and rt == 'paired-ended':
				continue
			if orga == 'mouse' and assay == 'DNase-seq' and rt == 'single-ended':
				continue
			if orga == 'human' and assay == 'DNase-seq' and rt == 'single-ended':
				continue
			combis.append('%s_%s_%s'%(orga, assay, rt))

# prepare the table for the final results
table = [['Organism', 'Assay', 'Run Type'] + feature_sets]
features_len = {}
features_auROC = {}

for features in feature_sets:
	features_len[features] = {}
	features_auROC[features] = {}
	# get metrics splitted in a list from command line 
	# (short names are always given)
	short_feature_names = features.split('-')
	# collect relevant col names used within grid search
	gs_feature_cols = []
	for short_name in short_feature_names:
		long_name = feature_abbr_sl[short_name]
		col_names = list(filter(lambda x: x[:len(long_name)] == long_name, 
						  data.columns.values))
		gs_feature_cols += col_names

	# imputation of missing values
	imp = SimpleImputer(missing_values=np.nan,strategy='median')
	for col in gs_feature_cols:
		data[col] = imp.fit_transform(data[[col]]).ravel()
	
	# load the generic model depending specific to the the feature set
	model_file_path = './generic_models/None_None_None_%s.model'%(features) 
	classifier = pickle.load(open(model_file_path, 'rb'))
	
	
	# intiate the objects needed by the grid search
	tenFold = StratifiedKFold(n_splits=10, random_state=the_RS, shuffle=True)
	
	# separate feature values and class values
	X = np.array(data.iloc[:, [True if col in gs_feature_cols else False 
							for col in data.columns.values]])
	y = np.array(data['status'])
	
	# initialize and run the ten-fold cross-validation
	trues = dict((combi, []) for combi in combis)
	prob_real = dict((combi, {'prob':[], 'real':[]}) for combi in combis)
	probas = cross_val_predict(classifier, X=X, y=y, 
							method='predict_proba', cv=tenFold)
	
	# for each sample, the real class label and the class probability 
	# predicted by the generic model is collected according to its metadata 
	for i in range(len(probas)):
		accession = data_accession[i]
		assay = data_assay[i]
		orga = data_org[i]
		rt = data_rt[i]
		real = data_status[i]
		prob = probas[i][1]
		pred = 1 if prob > 0.5 else 0
		trues['Various_Various_Various'].append(1 if real == pred else 0)
		trues['%s_Various_Various'%(orga)].append(1 if real == pred else 0)
		trues['Various_%s_Various'%(assay)].append(1 if real == pred else 0)
		trues['Various_Various_%s'%(rt)].append(1 if real == pred else 0)
		trues['%s_%s_Various'%(orga, assay)].append(1 if real == pred else 0)
		
		prob_real['Various_Various_Various']['prob'].append(prob)
		prob_real['%s_Various_Various'%(orga)]['prob'].append(prob)
		prob_real['Various_%s_Various'%(assay)]['prob'].append(prob)
		prob_real['Various_Various_%s'%(rt)]['prob'].append(prob)
		prob_real['%s_%s_Various'%(orga, assay)]['prob'].append(prob)
		
		prob_real['Various_Various_Various']['real'].append(real)
		prob_real['%s_Various_Various'%(orga)]['real'].append(real)
		prob_real['Various_%s_Various'%(assay)]['real'].append(real)
		prob_real['Various_Various_%s'%(rt)]['real'].append(real)
		prob_real['%s_%s_Various'%(orga, assay)]['real'].append(real)
		
		case = '%s_%s_%s'%(orga, assay, rt)
		if not case in trues:
			continue
		trues[case].append(1 if real == pred else 0)
		
		prob_real[case]['prob'].append(prob)
		prob_real[case]['real'].append(real)
	
	# after all probabilities and real values are collected,
	# the area under ROC curve is calculated for each combination
	for combi in combis:
		features_len[features][combi] = '%d'%(len(trues[combi]))
		features_auROC[features][combi] = roc_auc_score(prob_real[combi]['real'],
												  prob_real[combi]['prob'])
	print('%s was computed'%(features))

# create the final table, print it to the console and write into a text file
auROC_table = [['Organism', 'Assay', 'Run Type'] + feature_sets]
for combi in combis:
	row = list(combi.split('_'))
	for feature in feature_sets:
		row.append('%.3f'%features_auROC[feature][combi])
	auROC_table.append(row)

print_nice_table(auROC_table)

with open('auROCs_generic_model.tsv', 'w') as f:
	for row in auROC_table:
		f.write('\t'.join(row) + '\n')





















