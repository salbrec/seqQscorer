"""Parser utils

This script provides functions used by seqQscorer to parse input files.

Methods
-------

get_FastQC_features(feature_file_path)
	parses the RAW features from the FastQC tool
parse_BowtieSE(lines)
	parses the MAP features from Bowtie2 from single-end sequencing samples
parse_BowtiePE(lines)
	parses the MAP features from Bowtie2 from paired-end sequencing samples
get_Bowtie_features(feature_file_path)
	directly used by seqQscorer, parses Bowtie2 input and defines the run type
get_readsAnno_features(feature_file_path)
	parses the LOC features from ChIPseeker
get_TSS_features(feature_file_path)
	parses the TSS features from ChIPpeakAnno
generate_input_data(indir, feature_sets, run_type, medians, noVerbose)
	given the input directory this function reads in the feature sets for 
	all samples provided by the user

date:	2020-11-02
author:	Steffen Albrecht

"""

import os
import numpy as np
import pandas as pd

global FastQC_value_map
FastQC_value_map = {'FAIL': 0, 'WARN': 1, 'PASS': 2}


def get_RAW_features(feature_file_path):
	features = {}
	with open(feature_file_path, 'r') as feature_file:
		for line in feature_file:
			line = line.strip().split('\t')
			feature_name = line[1].replace(' ', '_')
			value = FastQC_value_map.get(line[0], np.nan)
			features['FastQC_'+feature_name] = value
	return features

def parse_BowtieSE(lines):
	lines = lines.split('\n')
	features = {}
	features['BowtieSE_no_mapping'] = float(lines[2].split('(')[1].split('%')[0])
	features['BowtieSE_uniquely'] = float(lines[3].split('(')[1].split('%')[0])
	features['BowtieSE_multiple'] = float(lines[4].split('(')[1].split('%')[0])
	features['BowtieSE_overall'] = float(lines[5].split('%')[0])
	# for mixed Bowtie
	features['BowtieMI_no_mapping'] = features['BowtieSE_no_mapping']
	features['BowtieMI_uniquely'] = features['BowtieSE_uniquely']
	features['BowtieMI_multiple'] = features['BowtieSE_multiple']
	features['BowtieMI_overall'] = features['BowtieSE_overall']
	return features

def parse_BowtiePE(lines):
	lines = lines.split('\n')
	features = {}
	features['BowtiePE_con_no_mapping'] = float(lines[2].split('(')[1].split('%')[0])
	features['BowtiePE_con_uniquely'] = float(lines[3].split('(')[1].split('%')[0])
	features['BowtiePE_con_multiple'] = float(lines[4].split('(')[1].split('%')[0])
	features['BowtiePE_dis_uniquely'] = float(lines[7].split('(')[1].split('%')[0])
	features['BowtiePE_cod_no_mapping'] = float(lines[11].split('(')[1].split('%')[0])
	features['BowtiePE_cod_uniquely'] = float(lines[12].split('(')[1].split('%')[0])
	features['BowtiePE_cod_multiple'] = float(lines[13].split('(')[1].split('%')[0])
	features['BowtiePE_overall'] = float(lines[14].split('%')[0])
	# for mixed Bowtie
	features['BowtieMI_no_mapping'] = features['BowtiePE_con_no_mapping']
	features['BowtieMI_uniquely'] = features['BowtiePE_con_uniquely']
	features['BowtieMI_multiple'] = features['BowtiePE_con_multiple']
	features['BowtieMI_overall'] = features['BowtiePE_overall']
	# for SE Bowtie
	features['BowtieSE_no_mapping'] = features['BowtiePE_con_no_mapping']
	features['BowtieSE_uniquely'] = features['BowtiePE_con_uniquely']
	features['BowtieSE_multiple'] = features['BowtiePE_con_multiple']
	features['BowtieSE_overall'] = features['BowtiePE_overall']
	return features

def get_MAP_features(feature_file_path):
	lines = open(feature_file_path, 'r').read()
	if 'concordantly' in lines and 'discordantly' in lines:
		return parse_BowtiePE(lines)
	else:
		return parse_BowtieSE(lines)

def get_LOC_features(feature_file_path):
	features = {}
	with open(feature_file_path, 'r') as f:
		f.readline()
		for line in f:
			line = line.strip().split('\t')
			feature_name = line[1]
			feature_name = feature_name.replace('"', '')
			feature_name = feature_name.replace("'", '')
			feature_name = feature_name.replace(' (<=300)', '')
			feature_name = feature_name.replace(' ', '_')
			feature_name = 'readsAnno_'+feature_name
			features[feature_name] = float(line[2])
	return features

def get_TSS_features(feature_file_path):
	tss = pd.read_csv(feature_file_path, sep='\t')
	tss_dist = list(map(str,tss['tss_dist']))
	feature_names = ['TSS_'+name if name[0] == '-' else 'TSS_+'+name for name in tss_dist]
	feature_values = list(tss['perc'])
	return dict(zip(feature_names, feature_values))


def generate_input_data(indir, feature_sets, run_type, medians, noVerbose=True, restrict=None):
	print('Parsing input data...')
	# initialize the specific parser functions
	parsers = {}
	parsers['RAW'] = get_RAW_features
	parsers['MAP'] = get_MAP_features
	parsers['LOC'] = get_LOC_features
	parsers['TSS'] = get_TSS_features
	
	# parse input data and create an input data frame
	if indir[-1] != '/':
		indir += '/'
	
	parsed_input = {}
	for subdir, dirs, files in os.walk(indir):
		for feature_file in files:
			file_path = indir + feature_file
			sample_ID = feature_file[:-4]
			if restrict != None:
				if restrict != sample_ID:
					continue
			feature_set = feature_file[-3:]
			if os.path.exists(file_path):
				if feature_file[-4:] in [ '.'+fs for fs in feature_sets ]:
					if not sample_ID in parsed_input:
						parsed_input[sample_ID] = {}
					
					features = parsers[feature_set](file_path)
					parsed_input[sample_ID].update(features)
	
	feature_prefix = {'RAW': 'FastQC', 'MAP': 'BowtieMI',
					'LOC': 'readsAnno', 'TSS': 'TSS'}
	if run_type != 'generic':
		feature_prefix['MAP'] = 'BowtieSE' if run_type == 'single-ended' else 'BowtiePE'

	feature_cols = []
	for abbr in feature_sets:
		prefix = feature_prefix[abbr]
		col_names = list(filter(lambda x: x[:len(prefix)] == prefix, medians.keys()))
		feature_cols += col_names
	feature_cols = sorted(feature_cols)
	
	missing = False
	input_data = dict( (col, []) for col in ['sampleID'] + feature_cols )
	for sample in parsed_input:
		input_data['sampleID'].append(sample)
		for col in feature_cols:
			if col in parsed_input[sample]:
				input_data[col].append(parsed_input[sample][col])
			else:
				missing = True
				input_data[col].append(np.nan)
				if col in feature_cols and not noVerbose:
					print('\nWarning! The feature "%s" is missing for %s'%(col, sample))
	if missing and not noVerbose:
		print('\nMissing values will be imputed by median.')
		print('However, you might check your input data.\n')
	
	input_data = pd.DataFrame(input_data)
	input_data = input_data[['sampleID'] + feature_cols]
	
	# search for missing features
	missing_features = list(filter(lambda x: not x in input_data.columns, feature_cols))
	for feature in missing_features:
		median = medians[feature]
		if not noVerbose:
			print('\nWarning: No values at all for the following feature:')
			print('\t%s'%(feature))
			print('\nFeature was imputed for all samples with its median value: "%s"\n'%(str(median)))
		input_data[feature] = input_data.shape[0] * median

	# imputation of missing values by the median value
	for feature in feature_cols:
		input_data[feature] = [medians[feature] if np.isnan(value) else value for value in input_data[feature]]
	
	print('... input data loaded.\n')
	return input_data, feature_cols



























