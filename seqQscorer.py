"""Machine Learning Quality Assessment of NGS Data

Main script for seqQscorer. "python seqQscorer.py --help" will display
a formatted help text on the console. A comprehensive description is provided
in the GitHub README that includes examples as well.

date:	2020-11-26
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

# parse command line arguments
script_dir = './'
if argv[0].find('/') >= 0:
	script_dir = argv[0][: - argv[0][::-1].find('/')]
utils_dir = '%sutils/'%(script_dir)

argsParser = argparse.ArgumentParser(description='seqQscorer - A machine learning application for quality assessment of NGS data')
argsParser.add_argument('--indir', '-i', type=str, required=True, help='Input directory containing the feature set files')
argsParser.add_argument('--species', '-s', type=str, default='generic', 
						choices=['generic','human', 'mouse'],  help='Species specifying the model used.')
argsParser.add_argument('--assay', '-a', type=str, default='generic', 
						choices=['generic','ChIP-seq','DNase-seq','RNA-seq'], help='Assay specifying the model used.')
argsParser.add_argument('--runtype', '-r', type=str, default='generic',
						choices=['generic','single-ended','paired-ended'], help='Run-Type specifying the model used.')
argsParser.add_argument('--model', '-m', type=str, default=None, help='Path to a serialized model, trained on own data. If used, the parameters --species, --assay, and --runtype have no impact on the classification model.')
argsParser.add_argument('--noRAW', action='store_true', help='Ignore all RAW features.')
argsParser.add_argument('--noMAP', action='store_true', help='Ignore all MAP features.')
argsParser.add_argument('--noLOC', action='store_true', help='Ignore all LOC features.')
argsParser.add_argument('--noTSS', action='store_true', help='Ignore all TSS features.')
argsParser.add_argument('--noFS', action='store_true', help='Switch off feature selection. (has only an impact if the best performance was achieved with chi2 or RFE)')
argsParser.add_argument('--bestCalib', action='store_true', help='Classifier setting is used that achieved the lowest brier score, hence the best calibration of the probabilities.')
argsParser.add_argument('--peaktype', '-pt', type=str, default=None, choices=['narrow','broad'], help='Optionally specify the peak-type for ChIP-seq data.')
argsParser.add_argument('--probOut', '-po', type=str, default=None,
						help='To specify an output file for the probabilities. Output will be tab-separated.')
argsParser.add_argument('--compOut', '-co', type=str, default=None,
						help='To specify an out file for the comprehensive output. Output will be kind of tab-separated.')
argsParser.add_argument('--inputOut', '-io', type=str, default=None,
						help='To specify an out file that will contain the parsed input. Output will be tab-separated.')
args = argsParser.parse_args()

if not os.path.isdir(args.indir):
	raise myExceptions.WrongFeatureInputException(
						'"%s" is not a directory'%(args.indir))

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
fs_suffix = '_noFS' if args.noFS else ''

# initiate the classification model and other data needed
species, assay, run_type = args.species, args.assay, args.runtype

if args.peaktype != None:
	if args.assay != 'ChIP-seq' or args.runtype != 'single-ended':
		raise myExceptions.WrongSettingException(
			'Peak-type specification can only be used for single-ended ChIP-seq.')
	assay = args.peaktype + assay

best_clf, feature_selection, selection, parameters, auROC, brier = utils.get_best_classifier(utils_dir, species, assay, 
															run_type, feature_sets, fs_suffix, model_sel_metric)

if best_clf == None and args.model == None:
	message = '''\nPlease check the given setting:
	assay:\t\t%s\n\tspecies:\t%s\n\trun-type:\t%s\n'''%(assay, species, run_type)
	message += '\tfeature sets:\t%s\n'%('-'.join(feature_sets))
	message += 'A specialized model for this setting is not available,\n'
	message += 'the generic model is used to proceed.\n'
	print(message)
	species, assay, run_type = 'generic', 'generic', 'generic'

application_case = '%s_%s_%s_%s'%(species, assay, run_type, '-'.join(feature_sets))
application_case += '_%s%s'%(model_sel_metric, fs_suffix)
model_file_path = '%smodels/%s.model'%(script_dir, application_case)

if args.model != None:
	print('An external model is provided.')
	model_file_path = args.model
	try: 
		pickle.load(open(model_file_path, 'rb'))
	except:
		raise myExceptions.IncorrectModelException(
			'The provided model from file "%s" could not be loaded.'%(model_file_path))

# load median values organized by subset, needed to impute missing values
medians = pickle.load(open('%sutils/medians.dict'%(script_dir), 'rb'))
medians = medians[species][assay][run_type]

# parse given input files
input_data, feature_columns = parser.generate_input_data(args.indir, feature_sets, run_type, medians)

# if the particula model is used for the very first time, it is trained and serialized
best_clf, feature_selection, selection, parameters, auROC, brier = utils.get_best_classifier(utils_dir, 
																	 species, assay, run_type, feature_sets, 
																	 fs_suffix, model_sel_metric)

if not os.path.exists(model_file_path) and args.model == None:
	print('\nThe required model was not used so far.')
	print('It needs to be trained and serialized...')
	
	clf = utils.get_clf_algos()[best_clf]
	if not best_clf in ['GNB','KNN']:
		parameters['random_state'] = 1
	
	clf_setup = clf.set_params(**parameters)
	
	data_file_path = '%sutils/datasets/%s_%s_%s.tsv'%(script_dir, assay, species, run_type)
	train_data = pd.read_csv(data_file_path, sep='\t')
	
	y = np.array(train_data['status'])
	train_data = train_data[feature_columns]
	if selection != None:
		train_data = train_data.loc[:,selection]
	X = np.array(train_data)
	
	model = clf_setup.fit(X,y)
	pickle.dump(model, open(model_file_path, 'wb'))
	
	print('... training and serialization is done!')
	print('The model is instantly available from now!')

# load the model
model = pickle.load(open(model_file_path, 'rb'))

# prepare input data format
input_values = input_data[feature_columns]
if selection != None:
	input_values = input_values.loc[:,selection]

# apply model on given samples to get the probabilities
probabilities = model.predict_proba(np.array(input_values))
fileIDs = list(input_data['sampleID'])
fileID_score = list(zip(fileIDs, [prob[1] for prob in probabilities]))

if args.model == None:
	print('\nThe best predictive performance was achived by %s'%(utils.clf_full_names(best_clf)))
	print('%s feature selection is applied (using %s%s of the features)'%(
		tuple(feature_selection.split('-') + ['%'])))
	print('Within the cross-validated gird search this model achived:')
	print('\tauROC: %s'%(auROC))
	print('\tBrier: %s'%(brier))
	print('\nThe model used, achieved these measures for different decision thresholds applied\non the probabilities within the grid-search (ten-fold cross-validation):')
	table = utils.read_in_measure_table(utils_dir, species, assay, run_type, 
										feature_sets, fs_suffix, model_sel_metric)
	utils.print_nice_table(table)
	print('')

# print the scores to the console
probas_str = ''
for fileID, score in sorted(fileID_score, key=lambda x: x[1]):
	probas_str += '%s\t%f\n'%(fileID, score)
	print(fileID, '%.3f '%(score), '(probability for being of low quality)', sep='\t')

# write probabilities to file if a file-path is given
if args.probOut != None:
	try:
		open(args.probOut, 'w').write(probas_str)
	except:
		raise myExceptions.WrongOutputFileException(
			'Unable to write the probabilities to file!')

# write comprehensive output to file if a file-path is given
if args.compOut != None:
	comp_out  = 'Model trained by: %s\n'%(utils.clf_full_names(best_clf))
	comp_out += '%s feature selection applied\n'%(feature_selection.split('-')[0])
	comp_out += '%s %s of the features are used\n'%(feature_selection.split('-')[1], '%')
	comp_out += 'auROC: %s\n'%(auROC)
	comp_out += 'Brier: %s\n\n'%(brier)
	comp_out += 'Metric table:\n'
	for row in table:
		comp_out += '\t'.join(row) + '\n'
	comp_out += '\n' + probas_str
	
	try:
		open(args.compOut, 'w').write(comp_out)
	except:
		raise myExceptions.WrongOutputFileException(
			'Unable to write the comprehensive output to file!')

# write the parsed input into a file, if a pth is given
if args.inputOut != None:
	try:
		input_data.to_csv(args.inputOut, sep='\t', index=False)
	except:
		raise myExceptions.WrongOutputFileException(
			'Unable to write the parsed input to file')




























































































