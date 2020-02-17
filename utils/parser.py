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

date:	2019-05-12
author:	Steffen Albrecht

"""

import os
import numpy as np
import pandas as pd

global FastQC_value_map
FastQC_value_map = {'FAIL': 0, 'WARN': 1, 'PASS': 2}


def get_FastQC_features(feature_file_path):
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

def get_Bowtie_features(feature_file_path):
	lines = open(feature_file_path, 'r').read()
	if 'concordantly' in lines and 'discordantly' in lines:
		return parse_BowtiePE(lines)
	else:
		return parse_BowtieSE(lines)

def get_readsAnno_features(feature_file_path):
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





























