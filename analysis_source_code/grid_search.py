"""Grid Search

This script runs the grid search for one of the ten classification algorithms
and for all algorithm-specific parameter setting combinations

Parameters
----------
organism : str
	human, mouse, or None
assay : str
	ChIP-seq, DNase-seq, RNA-seq, or None
runType : str
	single-ended, paired-ended, or None
features : str
	abbreviations for feature sets separated by "-"
	FQC, [BowM, BowS, BowP], RA, TSS
clf_setting : str
	classifier and feature selection setting separated by "-"
		classifier code: e.g. RFC for RandomForest
		feature selection method: chi2 or None
		percentage of features to keep: 100, 75, 50, or 25

-test (optional)
	runs a small grid for testing
-cores n_cores : int (optional)
	to specify the number of cores used to run the grid search

date:	2019-02-10
author:	Steffen Albrecht
example:
		python grid_search.py human ChIP-seq single-ended FQC RFC-None-100
		
"""

from sys import *
import random
import os
import importlib.util
import numpy as np
import pandas as pd
import copy
import warnings
import datetime
import pickle
import sklearn.pipeline as pipeline
from sklearn.feature_selection import SelectKBest, chi2

with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	from sklearn.model_selection import StratifiedKFold
	from sklearn.model_selection import GridSearchCV
	from sklearn import preprocessing

import utils.gs_utils as utils
import utils.ml_utils as ml_utils
import utils.grid as grid_collection

# define global variables for random and get the grid search settings
global the_RS, SEED
SEED = 1
the_RS = np.random.RandomState(SEED)
feature_abbr_ls, feature_abbr_sl = utils.get_abbr_maps()
classifiers, algorithms, param_grid, param_grid_test = grid_collection.get_grid()

# initialize log string and stop the time
start_time = datetime.datetime.now()
global log_str
log_str = '%s\n\n'%(' '.join(argv[1:]))
feature_selection_methods = {'chi2': chi2}

# command line arguments
organism = argv[1]
assay = argv[2]
runType = argv[3]
features = argv[4]
clf_setting = argv[5]
clf = clf_setting.split('-')[0]
fs_method = clf_setting.split('-')[1]
fs_K = int(clf_setting.split('-')[2])

# parse optional command line parameters
testing_mode = '-test' in argv
n_cores = 1 if testing_mode else 64
for i in range(len(argv)):
	if argv[i] == '-cores':
		n_cores = int(argv[i+1])

# define and create the folder that will contain the results
clf_dir = './gs_results'
if testing_mode:
	clf_dir += '/test'
case = '%s_%s_%s_%s'%(organism, assay, runType, features)
clf_dir = '%s/%s/%s/'%(clf_dir, case, clf_setting)
if not os.path.exists(clf_dir):
	os.makedirs(clf_dir)

# load the dataset specific for the given setting
data = pickle.load(open('./datasets/%s_%s_%s'%(assay,organism,runType),'rb'))

# split short names of feature sets into a list 
short_feature_names = features.split('-')
# collect relevant column names used within grid search
gs_feature_cols = []
for short_name in short_feature_names:
	long_name = feature_abbr_sl[short_name]
	col_names = list(filter(lambda x: x[:len(long_name)] == long_name, 
						 data.columns.values))
	gs_feature_cols += col_names

# get information about the dataset
revoked = len(list(filter(lambda x: x == 1, list(data['status']))))
released = len(list(filter(lambda x: x == 0, list(data['status']))))
n_samples = data.shape[0]

## initialize the objects needed by the grid search
tenFold = StratifiedKFold(n_splits=10, random_state=the_RS, shuffle=True)
param_grids = param_grid
scorer = ml_utils.get_scorer_dict()

# simplify the whole grid search (testing)
if testing_mode:
	tenFold = StratifiedKFold(n_splits=3, random_state=the_RS, shuffle=True)
	param_grids = param_grid_test

# measures will contain dataset information and validation measures (scores)
measures = dict((score, None) for score in scorer.keys())
measures['n'] = n_samples
measures['revoked'] = revoked
measures['released'] = released
measures['rev_perc'] = None
rev_perc = float(revoked) / float(n_samples) * 100.0
measures['rev_perc'] = rev_perc

# sort the column names to have it also sorted in the serialized models 
gs_feature_cols = sorted(gs_feature_cols)
# features X and class-labels y
X = np.array(data[gs_feature_cols])
y = np.array(data['status'])

# initialize the pipeline that runs the feature selection 
# and model training nested within the cross-validation
all_params = {}
for param in param_grids[clf]:
	all_params['clf__'+param] = param_grids[clf][param]
pipe = None
if fs_method == 'None':
	pipe = pipeline.Pipeline([('clf', algorithms[clf])])
else:
	k = int(float(fs_K) / 100.0 * len(gs_feature_cols))
	pipe = pipeline.Pipeline([('fs', SelectKBest()), ('clf', algorithms[clf])])
	all_params['fs__score_func'] = [feature_selection_methods[fs_method]]
	all_params['fs__k'] = [k]

# initialize and run the grid search
grid = GridSearchCV(pipe, all_params, 
					scoring=scorer, n_jobs=n_cores, cv=tenFold, 
					refit='auROC', return_train_score=True)
grid.fit(X, y)

# serialize all the objects: best model, the entire grid, 
# and all grid search results
pickle.dump(grid.best_estimator_, open('%sbest_estimator.model'%(clf_dir),
									   'wb'))
pickle.dump(grid, open('%sentire_grid.grid'%(clf_dir), 'wb'))
grid_results = pd.DataFrame(grid.cv_results_)
grid_results.to_csv('%sgrid_results.csv'%(clf_dir), sep='\t')

# collect the mean scoring values
best_setting = grid_results.iloc[grid.best_index_,:]
for score in sorted(scorer.keys()):
	measures[score] = best_setting['mean_test_'+score]

# finish the logging, print it to the console and write it into a text file
log_str += 'Best Parameters:\n'
log_str += '%s\n\n'%('\n'.join(['%s:\t%s'%(param, value) for 
								param, value in grid.best_params_.items()]))
end_time = datetime.datetime.now()
time = str(end_time - start_time)
time = time[:time.find('.')]
for m in measures:
	log_str += '%s\t%s\n'%(m, str(measures[m]))
log_str += '\nTime needed: %s\n'%time
print(log_str)
with open('%slog.txt'%(clf_dir), 'w') as f:
	f.write(log_str)

# serialize measures
pickle.dump(measures, open('%smeasures.dict'%(clf_dir), 'wb'))




