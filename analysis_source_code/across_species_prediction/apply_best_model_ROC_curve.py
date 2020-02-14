"""Cross-species generalization

Based on the best performing model and the different testing sets, this
script validates the model on the set of samples from the "unseen species".
Predictive performance is evaluated by ROC curves that are also plotted.

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
from itertools import combinations
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import interp

with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	from sklearn.metrics import roc_curve, auc

global utils
spec = importlib.util.spec_from_file_location("utils.py", "../utils/utils.py")
utils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(utils)

global ml_utils
spec = importlib.util.spec_from_file_location("ml_utils.py", "../utils/ml_utils.py")
ml_utils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ml_utils)

import grid_util as grid_collection

# define global variables for random and get the grid search settings
global the_RS, SEED
SEED = 1
the_RS = np.random.RandomState(SEED)
feature_abbr_ls, feature_abbr_sl = utils.get_abbr_maps()
classifiers, algorithms, param_grid, param_grid_test = grid_collection.get_grid()

# define settings for the datasets and feature sets
assays = ['ChIP-seq', 'DNase-seq']
organisms = ['human', 'mouse']
run_types = {'ChIP-seq': 'single-ended', 'DNase-seq': 'paired-ended'}
data_sets = {}
gs_feature_cols = []
metrics = ['FQC', 'Bow', 'RA', 'TSS']
feature_combis = sum([list(map(list, combinations(metrics,i))) for i in range(1,len(metrics)+1)], [])

# for both assays, apply best models, evaluate, and collect ROC-curve data
curves = {}
for assay in ['ChIP-seq', 'DNase-seq']:
	print(assay, '\n')
	curves[assay] = {}
	for fc in feature_combis:
		features = '-'.join(fc)
		if assay == 'ChIP-seq':
			features = features.replace('Bow','BowS')
		if assay == 'DNase-seq':
			features = features.replace('Bow','BowP')
		curves[assay][features] = []
		for orga_train in organisms:
			for orga_test in organisms:
				# used to collect ROC-curve data
				temp_curves = []
				tprs = []
				mean_fpr = np.linspace(0, 1, 100)
				
				for fold in range(5):
					# load the best classification model
					model_file = './best_models/%s_%s_%d_%s/best_estimator.model'%(assay, orga_train, fold, features)
					model = pickle.load(open(model_file, 'rb'))
				
					test_data = pd.read_csv('./datasets/%s_%s_%d_test.csv'%(assay, orga_test, fold), sep='\t')
					
					# delete those samples that have no value at all for a given quality metric
					gs_feature_cols = []
					for short_name in features.split('-'):
						long_name = feature_abbr_sl[short_name]
						col_names = list(filter(lambda x: x[:len(long_name)] == long_name, test_data.columns.values))
						gs_feature_cols += col_names
					
					# separate feature values and class values
					X = np.array(test_data.iloc[:, [True if col in gs_feature_cols else False for col in test_data.columns.values]])
					y = np.array(test_data['status'])
					
					# for the current fold, create the ROC-curve
					probas = model.predict_proba(X)
					fpr, tpr, thresholds = roc_curve(y, probas[:,1],pos_label=1)
					auROC = auc(fpr, tpr)
					tprs.append(interp(mean_fpr, fpr, tpr))
					tprs[-1][0] = 0.0
					
					print(orga_train, '->', orga_test, fold, '%.3f'%auROC)
					temp_curves.append((orga_train, orga_test, fpr, tpr, auROC))
				median_curve = sorted(temp_curves, key=lambda x: x[-1])[2]
				print(orga_train, '->', orga_test, 'median:\t', '%.3f'%median_curve[-1], '\tmean:  ', np.mean([c[-1] for c in temp_curves]))
				
				# collect data for the curve
				mean_tpr = np.mean(tprs, axis=0)
				mean_tpr[-1] = 1.0
				curve_data = (orga_train, orga_test, mean_fpr, mean_tpr, auc(mean_fpr, mean_tpr))
				curves[assay][features].append(curve_data)
				print('')

# define plot properties and plot the ROC curves
fig_size_w = 8.2
fig_size_h = 11.6
curves_width = 2
colors = {}
colors['mouse'] = 'green'
colors['human'] = 'darkblue'
short_org = {'human':'HS', 'mouse':'MM'}
fs_legend = 12
fs_labels = 12
fs_ticks = 10
ticks = [x/10 for x in range(0,12,2)]
def plot_ROC_curve(assay, features, prefix_index, curves):
	fig = plt.figure(figsize=(5,5))
	for curve in sorted(curves, key=lambda x: x[-1], reverse=True):
		label = '%s -> %s'%(short_org[curve[0]], short_org[curve[1]])
		label += ' %.2f'%(curve[-1])
		ls = '-' if curve[0] == curve[1] else '--'
		plt.plot(curve[2], curve[3], 
				color=colors[curve[1]], 
				linewidth=curves_width, linestyle=ls,
				label=label)
	plt.plot([0, 1], [0, 1], 'k--', lw=1)
	plt.xlim([0.0, 1.01])
	plt.ylim([0.0, 1.01])
	plt.xlabel('FPR', fontsize=fs_labels)
	plt.ylabel('TPR', fontsize=fs_labels)
	plt.xticks(fontsize=fs_ticks, ticks=ticks)
	plt.yticks(fontsize=fs_ticks, ticks=ticks)
	plt.title('%s\n%s'%(assay, features))
	plt.legend(loc="lower right", fontsize=fs_legend)
	plt.tight_layout()
	fig.savefig('./ROC_curves/%s_%d_%s.svg'%(assay, prefix_index, features))
for assay in assays:
	for prefix_index, fc in enumerate(feature_combis):
		features = '-'.join(fc)
		if assay == 'ChIP-seq':
			features = features.replace('Bow','BowS')
		if assay == 'DNase-seq':
			features = features.replace('Bow','BowP')
		plot_ROC_curve(assay, features, prefix_index, curves[assay][features])
for assay in assays:
	for fc in feature_combis:
		features = '-'.join(fc)
		print(assay, features)



































