"""Train a new model on your own labeled data.

The user must provied preprocessed NGS samples and quality labels
in best case manually curated. This script can then be used to train
a new model on this data. The model can be used as well with seqQscorer 
(--model). A use case is provided in the github repository together with
comprehensive descriptions and example runs in the main README.

date:	2020-12-04
author:	Steffen Albrecht

"""

from sys import *
import os
import pickle
import pandas as pd
import numpy as np
import json
import random
import argparse

import warnings
warnings.filterwarnings("ignore")

# import project utils
import utils.Exceptions as myExceptions
import utils.utils as utils
import utils.parser as parser
import utils.custom_metrics as cm

from sklearn.model_selection import cross_validate, StratifiedKFold
from sklearn.metrics import roc_auc_score, precision_recall_curve, auc, precision_score, recall_score, f1_score, accuracy_score

# parse command line arguments
script_dir = './'
if argv[0].find('/') >= 0:
	script_dir = argv[0][: - argv[0][::-1].find('/')]
utils_dir = '%sutils/'%(script_dir)

argsParser = argparse.ArgumentParser(description='Train a new model based on your labeled input data.')
argsParser.add_argument('--training', '-t', type=str, required=True, help='Input directory containing the feature set files.')
argsParser.add_argument('--labels', '-l', type=str, required=True, help='Table (tab-separated) that contains the labels. It has to contain one column called "sampleID" for the IDs that have to fit to the files within "training". A second column should be therer called "quality" describing low and high quality by 1 and 0, respectively. It is also possible to use "--column" to chose a column within "labels" different to "quality".')
argsParser.add_argument('--model', '-m', type=str, required=True, help='File path for saving the model.')
argsParser.add_argument('--column', '-c', type=str, default='quality', required=False, help='Specifies the column within the file with the labels.')
argsParser.add_argument('--cvFolds', '-cv', type=int, default=10, help='Specifies the number of folds for the cross-validation.')
argsParser.add_argument('--seed', '-rs', type=int, default=1, help='Used to set the random state of the classifier and the cross-validation.')
argsParser.add_argument('--species', '-s', type=str, default='generic', 
						choices=['generic','human', 'mouse'],  help='Species specifying the model used.')
argsParser.add_argument('--assay', '-a', type=str, default='generic', 
						choices=['generic','ChIP-seq','DNase-seq','RNA-seq'], help='Assay specifying the model used.')
argsParser.add_argument('--runtype', '-r', type=str, default='generic',
						choices=['generic','single-end','paired-end'], help='Run-Type specifying the model used.')
argsParser.add_argument('--bestCalib', action='store_true', help='Classifier setting is used that achieved the lowest brier score, hence the best calibration of the probabilities.')
argsParser.add_argument('--inputOut', '-io', type=str, default=None,
						help='To specify an out file that will contain the parsed input together with the quaity labels. The output file is exactly the data used to train the model Output will be tab-separated.')
argsParser.add_argument('--noRAW', action='store_true', help='Ignore all RAW features.')
argsParser.add_argument('--noMAP', action='store_true', help='Ignore all MAP features.')
argsParser.add_argument('--noLOC', action='store_true', help='Ignore all LOC features.')
argsParser.add_argument('--noTSS', action='store_true', help='Ignore all TSS features.')
argsParser.add_argument('--useRF', '-rf', type=str, default=None, help='Simply use Random Forest. Specify its parameters via a ":"-separated parameter: "criterion:maxDepth:maxFeatures:nTrees"')
args = argsParser.parse_args()

feature_sets = ['RAW','MAP','LOC','TSS']

# restrict feature sets used according to given optional parameters
if args.noRAW:
	feature_sets.remove('RAW')
if args.noMAP:
	feature_sets.remove('MAP')
if args.noLOC:
	feature_sets.remove('LOC')
if args.noTSS:
	feature_sets.remove('TSS')

model_sel_metric = 'brier' if args.bestCalib else 'auROC'

# initiate the classification model and other data needed
species, assay, run_type = args.species, args.assay, args.runtype
best_clf, feature_selection, selection, parameters, auROC, brier = utils.get_best_classifier(utils_dir, species, assay, 
															run_type, feature_sets, '_noFS', model_sel_metric)

if best_clf == None:
	message = '''\nPlease check the given setting:
	assay:\t\t%s\n\tspecies:\t%s\n\trun-type:\t%s\n'''%(assay, species, run_type)
	message += '\tfeature sets:\t%s\n'%('-'.join(feature_sets))
	message += 'There is no recommended algorithm-setting suggestion for this setting,\n'
	message += 'the recommendation for the generic data is used to proceed.\n'
	print(message)
	species, assay, run_type = 'generic', 'generic', 'generic'

application_case = '%s_%s_%s_%s'%(species, assay, run_type, '-'.join(feature_sets))
application_case += '%s%s'%(model_sel_metric, '_noFS')
model_file_path = '%smodels/%s.model'%(script_dir, application_case)

# load median values organized by subset, needed to impute missing values
medians = pickle.load(open('%sutils/medians.dict'%(script_dir), 'rb'))
medians = medians[species][assay][run_type]

# parse given input files
input_data, feature_columns = parser.generate_input_data(args.training, feature_sets, run_type, medians)

# if the particula model is used for the very first time, it is trained and serialized
best_clf, feature_selection, selection, parameters, auROC, brier = utils.get_best_classifier(utils_dir, 
																	 species, assay, run_type, feature_sets, 
																	 '_noFS', model_sel_metric)

# setup the classifier
clf_setup = None
if args.useRF == None:
	print('Using the %s classifier with appropriate parameters'%(utils.clf_full_names(best_clf)))
	print('that performed well for the given species-assay-runtype specification.\n')
	clf = utils.get_clf_algos()[best_clf]
	if not best_clf in ['GNB','KNN']:
		parameters['random_state'] = args.seed
	clf_setup = clf.set_params(**parameters)
else:
	print('Using the Random Forest Classifier specified by the user.\n')
	clf = utils.get_clf_algos()['RFC']
	rf_args = args.useRF.split(':')
	
	max_depth = None
	if rf_args[1] != 'None' and rf_args[1] != '-1':
		max_depth = int(rf_args[1])
	
	max_features = 'auto'
	if rf_args[2] not in ['None', '-1', 'auto', 'sqrt', 'log2']:
		max_features = float(rf_args[2])
	else:
		if rf_args[2] in ['auto', 'sqrt', 'log2']:
			max_features = rf_args[2]
	
	params = {'criterion':rf_args[0], 'max_depth':max_depth}
	params.update({'max_features':max_features, 'n_estimators':int(rf_args[3])})
	params['random_state'] = args.seed
	
	clf_setup = clf.set_params(**params)

# collect labels for the samples in input data
labels = pd.read_csv(args.labels, sep='\t')
labels.dropna(subset=[args.column], inplace=True)
label_map = dict( zip( labels['sampleID'], labels[args.column] ) )
y = []
for sample in input_data['sampleID']:
	if sample in label_map:
		y.append(label_map[sample])
y = np.array(y)

# filter input data for samples without a label
input_data = input_data.loc[ [True if sample in label_map else False for sample in input_data['sampleID']]]

input_data = input_data[feature_columns]

X = np.array(input_data)

print('Input data has %d samples and %d features.'%(input_data.shape))
print('Fraction of positive labels: %.2f'%(np.mean(y)))

# training the model
print('Training and cross-validating the model now...')
stratifiedFold = StratifiedKFold(n_splits=args.cvFolds, random_state=args.seed, shuffle=True)

cv_auROC = []
cv_auPRC = []

# collect important evaluation measures during the grid search to 
# provide more information about which decision threshold to use
metrics = {}
metrics['Precision'] = dict( (dt, []) for dt in range(1,10) )
metrics['Recall']    = dict( (dt, []) for dt in range(1,10) )
metrics['F1']        = dict( (dt, []) for dt in range(1,10) )
metrics['Accuracy']  = dict( (dt, []) for dt in range(1,10) )
for train, test in stratifiedFold.split(X, y):
	model = clf_setup.fit(X[train],y[train])
	
	probas = model.predict_proba(X[test])
	
	auROC = roc_auc_score(y[test], probas[:,1])
	prec, rec, dts = precision_recall_curve(y[test], probas[:,1], pos_label=1)
	auPRC = auc(rec, prec)
	
	cv_auROC.append(auROC)
	cv_auPRC.append(auPRC)
	
	# for the different DTs
	for dt in range(1,10):
		threshold = float(dt) / 10.0
		y_pred = [1 if prob > threshold else 0 for prob in probas[:,1]]
		
		metrics['Precision'][dt] = precision_score(y[test], y_pred) 
		metrics['Recall'][dt] = recall_score(y[test], y_pred)
		metrics['F1'][dt] = f1_score(y[test], y_pred)
		metrics['Accuracy'][dt] = accuracy_score(y[test], y_pred)

print('The cross-validated performance metrics are:')
print('\tauROC: %.3f'%np.mean(cv_auROC))
print('\tauPRC: %.3f'%np.mean(cv_auPRC))
print('')
print('See also more metrics for different decision thresholds:')

header = ['Decision Threshold:']
header += list( map(lambda x: '%.1f'%(float(x)/10.0), range(1,10) ))
table = [ header ]
for metric in metrics:
	row = [metric]
	row += [ '%.3f'%(np.mean(metrics[metric][dt])) for dt in range(1,10) ]
	table.append(row)
utils.print_nice_table(table)

# train again all all samples, then serialize
model = clf_setup.fit(X,y)
pickle.dump(model, open(args.model, 'wb'))

# write input data to a file
if args.inputOut != None:
	input_data['quality'] = y
	input_data.to_csv(args.inputOut, sep='\t', index=False)





















