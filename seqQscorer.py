"""Machine Learning Quality Assessment of NGS Data

Main script for seqQscorer. "python seqQscorer.py --help" will display
a formatted help text on the console. A comprehensive description is provided
in the GitHub README that includes examples as well.

date:	2019-05-12
author:	Steffen Albrecht

"""

from sys import *
import os
import pickle
import pandas as pd
import numpy as np
import random
import argparse
from sklearn.impute import SimpleImputer

# import project utils
import utils.parser as parser
import utils.Exceptions as myExceptions
import utils.utils as utils

feature_path = {}
feature_path['fastqc'] = None
feature_path['bowtie'] = None
feature_path['readsAnno'] = None
feature_path['tss'] = None

species = 'None'
assay = 'None'
run_type = 'None'
features = []

out_file_path = None

# parse command line arguments
script_directory = './'
if argv[0].find('/') >= 0:
	script_directory = argv[0][: - argv[0][::-1].find('/')]

if '--help' in argv:
	print(utils.get_help_text())
	raise SystemExit(0)

if '--rt' in argv:
	run_type = argv[ argv.index('--rt') + 1 ]

if '--spec' in argv:
	species = argv[ argv.index('--spec') + 1 ]

if '--assay' in argv:
	assay = argv[ argv.index('--assay') + 1 ]

if '--raw' in argv:
	feature_path['fastqc'] = argv[ argv.index('--raw') + 1 ]
	features.append('RAW')

if '--map' in argv:
	feature_path['bowtie'] = argv[ argv.index('--map') + 1 ]
	if run_type != 'single-ended' and run_type != 'paired-ended':
		features.append('MAPm')
	if run_type == 'single-ended':
		features.append('MAPs')
	if run_type == 'paired-ended':
		features.append('MAPp')

if '--loc' in argv:
	feature_path['readsAnno'] = argv[ argv.index('--loc') + 1 ]
	features.append('LOC')

if '--tss' in argv:
	feature_path['tss'] = argv[ argv.index('--tss') + 1 ]
	features.append('TSS')

if '--out' in argv:
	out_file_path = argv[ argv.index('--out') + 1 ]

# collect all file paths for the given input
feature_files = {}
for feature_type in feature_path.keys():
	path = feature_path[feature_type]
	if path == None:
		continue
	feature_files[feature_type] = []
	if os.path.isdir(path):
		for subdir, dirs, files in os.walk(path):
			for feature_file in files:
				if os.path.isfile(path + feature_file):
					feature_files[feature_type].append(path + feature_file)
				else:
					raise myExceptions.WrongFeatureInputException(
						'"%s" in directory "%s" is not a file.'%(
							feature_file, path))
	else:
		if os.path.isfile(path):
			feature_files[feature_type].append(path)
		else:
			raise myExceptions.WrongFeatureInputException(
				'"%s" is neither a file nor a directory.'%(path))

# initialize the specific parser functions
parsers = {}
parsers['fastqc'] = parser.get_FastQC_features
parsers['bowtie'] = parser.get_Bowtie_features
parsers['readsAnno'] = parser.get_readsAnno_features
parsers['tss'] = parser.get_TSS_features

# parse the given input files to extract all feature values
all_IDs = set()
all_feature_names = set()
fileID_features = {}
for feature_type, file_paths in feature_files.items():
	for feature_file in file_paths:
		folder, file_name_ending, fileID = utils.get_path_info(feature_file)
		if not fileID in fileID_features:
			fileID_features[fileID] = {}
		temp_features = parsers[feature_type](feature_file)
		all_IDs.add(fileID)
		all_feature_names |= set(temp_features.keys())
		fileID_features[fileID].update(temp_features)

# load the appropriate model
case = '%s_%s_%s_%s'%(str(species), str(assay), str(run_type), '-'.join(features))
best_model_file_path = '%smodels/%s.model'%(script_directory, case)
if not os.path.exists(best_model_file_path):
	message = '''\nPlease check the given setting:
	assay:\t\t%s\n\tspecies:\t%s\n\trun-type:\t%s\n'''%(assay, species, run_type)
	message += '\tfeature sets:\t%s\n'%('-'.join(features))
	message += 'A specialized model for this setting is not available,\n'
	message += 'the generic model is used to proceed.\n'
	print(message)
	features = ['MAPm' if f in ['MAPs', 'MAPp'] else f for f in features]
	species, assay, run_type = None, None, None
	case = '%s_%s_%s_%s'%(str(species), str(assay), str(run_type), '-'.join(features))
	best_model_file_path = '%smodels/%s.model'%(script_directory, case)

# prepare a pandas DataFrame containing the input data
all_feature_names = sorted(list(all_feature_names), reverse=True)
input_data = dict((feature_name, []) for feature_name in all_feature_names) 
input_data['accession'] = []
for accession in fileID_features:
	input_data['accession'].append(accession)
	for feature_name in all_feature_names:
		input_data[feature_name].append(fileID_features[accession].get(feature_name, np.nan))
input_data = pd.DataFrame(input_data)

# load median values organized by subset, needed to impute missing values
medians = pickle.load(open('%sutils/medians.dict'%(script_directory), 'rb'))
medians = medians[str(species)][str(assay)][str(run_type)]

# define all columns that are considered as features
feature_abbr_sl = {'RAW': 'FastQC', 'MAPs': 'BowtieSE',
				   'MAPp': 'BowtiePE', 'MAPm': 'BowtieMI',
				   'LOC': 'readsAnno', 'TSS': 'TSS'}
feature_cols = []
for short_name in features:
	long_name = feature_abbr_sl[short_name]
	col_names = list(filter(lambda x: x[:len(long_name)] == long_name, medians.keys()))
	feature_cols += col_names
feature_cols = sorted(feature_cols) # the sorting is also done in the grid search 

# search for missing features
missing_features = list(filter(lambda x: not x in input_data.columns, feature_cols))
for feature in missing_features:
	median = medians[feature]
	print('\nWarning: No values at all for the following feature:')
	print('\t%s'%(feature))
	print('\nFeature was imputed for all samples with its mean value: "%s"\n'%(str(median)))
	input_data[feature] = input_data.shape[0] * median
input_data = input_data[feature_cols + ['accession']]

# imputation of missing values by the median value
for feature in feature_cols:
	input_data[feature] = [medians[feature] if np.isnan(value) else value for value in input_data[feature]]

# apply model on given samples to get the probabilities
model = pickle.load(open(best_model_file_path, 'rb'))
probabilities = model.predict_proba(input_data[feature_cols])
fileIDs = list(input_data['accession'])
fileID_score = list(zip(fileIDs, [prob[1] for prob in probabilities]))

# print the scores to the console
for fileID, score in sorted(fileID_score, key=lambda x: x[1]):
	print(fileID, '%.3f '%(score), '(probability for being of low quality)', sep='\t')

# write scores to file
if out_file_path != None:
	with open(out_file_path, 'w') as f:
		for fileID, score in sorted(fileID_score, key=lambda x: x[1]):
			f.write('%s\t%f\t(probability for being of low quality)\n'%(fileID, score))













