"""Cross-species generalization

A grid search applied to the datasets, specific for this analysis
to find the best performing classification model

Parameters
----------
assay : str
	ChIP-seq or DNase-seq
organism : str
	human, mouse
fold : int
	fold-index of the 5-fold cross-validation
features : str 
	abbreviations for feature sets separated by "-"

date:	2019-04-08
author:	Steffen Albrecht

"""

from sys import *
import random
import requests, json
import os
import importlib.util
import numpy as np
import pandas as pd
import copy
import warnings
import datetime
import pickle
import sklearn.pipeline as pipeline

with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	from sklearn.model_selection import StratifiedKFold
	from sklearn.model_selection import GridSearchCV
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

# parse optional command line arguments
testing_mode = '-test' in argv
n_cores = 1 if testing_mode else 64
for i in range(len(argv)):
	if argv[i] == '-cores':
		n_cores = int(argv[i+1])

# define some global variables and all the settings for the grid search
global the_RS, SEED
SEED = 1
the_RS = np.random.RandomState(SEED)
feature_abbr_ls, feature_abbr_sl = utils.get_abbr_maps()
classifiers, algorithms, param_grid, param_grid_test = grid_collection.get_grid()

if testing_mode:
	param_grid = param_grid_test

assay = argv[1]
organism = argv[2]
fold = int(argv[3])
features = argv[4]

data = pd.read_csv('./datasets/%s_%s_%d_train.csv'%(assay, organism, fold), sep='\t')
run_types = {'ChIP-seq': 'single-ended', 'DNase-seq': 'paired-ended'}

# collect relevant column names
short_feature_names = features.split('-')
gs_feature_cols = []
for short_name in short_feature_names:
	long_name = feature_abbr_sl[short_name]
	col_names = list(filter(lambda x: x[:len(long_name)] == long_name, data.columns.values))
	gs_feature_cols += col_names

# separate feature values and class values
gs_feature_cols = sorted(gs_feature_cols)
X = np.array(data.iloc[:, [True if col in gs_feature_cols else False for col in data.columns.values]])
y = np.array(data['status'])
scorer = ml_utils.get_scorer_dict()

# run the grid search for each algorithm to identify the best
# algorithm with its best setting
best = []
for algo in algorithms:
	print(algo)
	print(param_grid[algo])
	
	stratTFcv = StratifiedKFold(n_splits=10, random_state=the_RS, shuffle=True)
	measures = dict((score, None) for score in scorer.keys())
	target_dir = './models/%s_%s_%d_%s_%s/'%(assay, organism, fold, features, algo)
	if not os.path.exists(target_dir):
		os.mkdir(target_dir)
	if os.path.exists('%smeasures.dict'%(target_dir)):
		if os.path.exists('%sbest_estimator.model'%(target_dir)):
			continue
	
	try:
		grid = GridSearchCV(algorithms[algo], param_grid[algo], 
						scoring=scorer, n_jobs=n_cores, cv=stratTFcv, 
						refit='auROC', return_train_score=True)
		grid.fit(X,y)
		
		pickle.dump(grid.best_estimator_, open('%sbest_estimator.model'%(target_dir), 'wb'))
		grid_results = pd.DataFrame(grid.cv_results_)

		# collect the scoring values
		best_setting = grid_results.iloc[grid.best_index_,:]
		for score in sorted(scorer.keys()):
			measures[score] = best_setting['mean_test_'+score]
		
		pickle.dump(measures, open('%smeasures.dict'%(target_dir), 'wb'))
		
		with open('%smeasures.txt'%(target_dir), 'w') as f:
			for score in sorted(scorer.keys()):
				f.write('%s\t%.5f\n'%(score, measures[score]))
		
		best.append((algo, grid.best_estimator_, measures))
	except:
		print('%s did not work properly'%algo)

# serialize best classification model
best = sorted(best, key=lambda x: x[2]['auROC'], reverse=True)
target_dir = './best_models/%s_%s_%d_%s/'%(assay, organism, fold, features)
if not os.path.exists(target_dir):
	os.mkdir(target_dir)
pickle.dump(best[0][1], open('%sbest_estimator.model'%(target_dir), 'wb'))
pickle.dump(best[0][2], open('%smeasures.dict'%(target_dir), 'wb'))

# write all validation measures into a text file  
with open('%smeasures.txt'%(target_dir), 'w') as f:
	for score in sorted(scorer.keys()):
		f.write('%s\t%.5f\n'%(score, best[0][2][score]))


