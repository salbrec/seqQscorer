"""Cross-species generalization

Within this script, the 5-fold cross-validation training and testing 
sets are created and serialized 

Parameters
----------
assay : str
	ChIP-seq or DNase-seq

date:	2019-04-08
author:	Steffen Albrecht

"""

from sys import *
import random
import os
import numpy as np
import pandas as pd
import copy
import importlib.util
import warnings

with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	from sklearn.model_selection import StratifiedKFold
	from sklearn.impute import SimpleImputer

global utils
spec = importlib.util.spec_from_file_location("utils.py", "../utils/utils.py")
utils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(utils)

global ml_utils
spec = importlib.util.spec_from_file_location("ml_utils.py", "../utils/ml_utils.py")
ml_utils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ml_utils)

import grid_util as grid_collection

# define variables for random and 
global the_RS, SEED
SEED = 1
the_RS = np.random.RandomState(SEED)
feature_abbr_ls, feature_abbr_sl = utils.get_abbr_maps()
classifiers, algorithms, param_grid, param_grid_test = grid_collection.get_grid()
organisms = ['human', 'mouse']
run_types = {'ChIP-seq': 'single-ended', 'DNase-seq': 'paired-ended'}

assay = argv[1]

# add another mode: using only samples having all features: FQC + Bow + RA
# for paired-ended cases using only FQC, half of the samples are kicked out 
run_type = run_types[assay]

data_sets = {}
gs_feature_cols = []

for organism in organisms:
	# load and specify dataset
	data = pd.read_csv('../datasets/raw_fastq_samples.csv', sep=',')
	if organism != 'None':
		data = data.loc[data['organism'] == organism]
	if assay != 'None':
		data = data.loc[data['assay'] == assay]
	if run_type != 'None':
		data = data.loc[data['runType'] == run_type]
	
	# delete those samples that have no value at all for a given quality metric
	# the setting "all_features_avail" might be interesting here
	defining_the_set_features = ['FQC']
	if run_type == 'None':
		defining_the_set_features.append('BowM')
	if run_type == 'single-ended':
		defining_the_set_features.append('BowS')
	if run_type == 'paired-ended':
		defining_the_set_features.append('BowP')
	defining_the_set_features.append('RA')
	defining_the_set_features.append('TSS')
	
	# collect column names relevant for the grid search
	gs_feature_cols = []
	for short_name in defining_the_set_features:
		long_name = feature_abbr_sl[short_name]
		col_names = list(filter(lambda x: x[:len(long_name)] == long_name, data.columns.values))
		data = data.dropna(axis='index', how='all', subset=col_names)
		gs_feature_cols += col_names
	
	# imputation of missing values
	imp = SimpleImputer(missing_values=np.nan,strategy='median')
	for col in gs_feature_cols:
		data[col] = imp.fit_transform(data[[col]]).ravel()
	
	# depending on the features selected, data can be unbalanced again
	data = ml_utils.filter_50_50(data, seed=SEED)
	data_sets[organism] = data.copy()

# for each organism (human and mouse), create startified and balanced 
# training and testing sets
n_min = min([data_sets[orga].shape[0] for orga in data_sets])
for organism in data_sets:
	if data_sets[organism].shape[0] > n_min:
		data_sets[organism] = ml_utils.filter_50_50(data_sets[organism],seed=SEED,n=n_min)
for organism in data_sets:
	data = data_sets[organism]
	print(organism, data.shape, sum(data['status'] == 1), sum(data['status'] == 0))
for organism in data_sets:
	data = data_sets[organism]
	X = np.array(data.iloc[:, [True if col in gs_feature_cols else False for col in data.columns.values]])
	y = np.array(data['status'])
	# create and serialize stratified training and testing sets 
	# based on a 5-fold cross-validation strategy
	stratCV = StratifiedKFold(n_splits=5, random_state=the_RS, shuffle=True)
	for i, (train_iloc, test_iloc) in enumerate(stratCV.split(X, y)):
		train = data.iloc[train_iloc]
		test  = data.iloc[test_iloc]
		train.to_csv('./datasets/%s_%s_%d_train.csv'%(assay, organism, i), sep='\t')
		test.to_csv('./datasets/%s_%s_%d_test.csv'%(assay, organism, i), sep='\t')





































