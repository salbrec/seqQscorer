"""With experiment differences

This script performs the analysis for the within-experiment differences. 
First, experiments are identified having both released and revoked files. 
Then a classification model is trained without a set of files from the 
same experiment. This model is used to calculate a class probability for 
the files that were excluded before. A boxplot, showing the distribution 
of these probabilities for released and revoked files, is directly created. 
All class probabilities are serialized for further plots.  

date:	2019-06-15
author:	Steffen Albrecht

"""

from sys import *
import random
import os
import numpy as np
import pandas as pd
import copy
import warnings
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
import utils.gs_utils as utils
with warnings.catch_warnings():
	from sklearn.impute import SimpleImputer

# set variables for random and load the feature set names
global the_RS, SEED
SEED = 1
the_RS = np.random.RandomState(SEED)
feature_abbr_ls, feature_abbr_sl = utils.get_abbr_maps()

# read in metadata for fastq files
files = pd.read_csv('./datasets/fastq_files.tsv', sep='\t')
acc_exp = dict(zip(files['accession'], files['experiment']))
acc_status = dict(zip(files['accession'], files['status']))

# load data and prepare directories that contain the results
all_data = pd.read_csv('./datasets/raw_fastq_samples.csv', sep=',')
the_features = 'FQC'
classifier_dir = '../models/'
tables_dir = './within_exp_tables/%s/'%(the_features)
if not os.path.exists(tables_dir):
	os.makedirs(tables_dir)
dot_plot_tables_dir = './within_exp_tables/tsv_tables/%s/'%(the_features)
if not os.path.exists(dot_plot_tables_dir):
	os.makedirs(dot_plot_tables_dir)

settings = []
settings.append(('human', 'ChIP-seq', 'single-ended'))
settings.append(('mouse', 'ChIP-seq', 'single-ended'))
settings.append(('human', 'DNase-seq', 'paired-ended'))
settings.append(('mouse', 'DNase-seq', 'paired-ended'))

# later used to create the boxplot
boxplot_data = {'Setting': [], 'prob': [], 'Status': []}
label_order = []

for setting in settings:
	# create a subset specific to the current setting
	organism = setting[0]
	assay = setting[1]
	run_type = setting[2]
	data = all_data.copy()
	if organism != 'None':
		data = data.loc[data['organism'] == organism]
	if assay != 'None':
		data = data.loc[data['assay'] == assay]
	if run_type != 'None':
		data = data.loc[data['runType'] == run_type]
	# adapt feature names to run type
	features = the_features
	if run_type == 'single-ended':
		features = features.replace('Bow','BowS')
	if run_type == 'paired-ended':
		features = features.replace('Bow','BowP')
	
	# collect relevant column names
	short_feature_names = features.split('-')
	classification_feature_cols = []
	for short_name in short_feature_names:
		long_name = feature_abbr_sl[short_name]
		col_names = list(filter(lambda x: x[:len(long_name)] == long_name, 
						  data.columns.values))
		classification_feature_cols += col_names

	# delete those samples that have no value at all for a given quality metric
	# impute missing values
	defining_the_set_features = short_feature_names
	defining_the_set_features = ['FQC']
	if run_type == 'None':
		defining_the_set_features.append('BowM')
	if run_type == 'single-ended':
		defining_the_set_features.append('BowS')
	if run_type == 'paired-ended':
		defining_the_set_features.append('BowP')
	defining_the_set_features.append('RA')
	defining_the_set_features.append('TSS')
	for short_name in defining_the_set_features:
		long_name = feature_abbr_sl[short_name]
		col_names = list(filter(lambda x: x[:len(long_name)] == long_name, 
						  data.columns.values))
		data = data.dropna(axis='index', how='all', subset=col_names)
	imp = SimpleImputer(missing_values=np.nan,strategy='median')
	for col in classification_feature_cols:
		data[col] = imp.fit_transform(data[[col]]).ravel()
	
	print(' '.join(setting))
	print('Shape of the data:', data.shape)
	
	# parse given experiments and select the relevant ones
	exp_info = {}
	for accession in data['accession']:
		exp = acc_exp[accession]
		file_status = acc_status[accession]
		if not exp in exp_info:
			exp_info[exp] = {'accession': exp, 'files': [], 'status': []}
		exp_info[exp]['files'].append(accession)
		exp_info[exp]['status'].append(file_status)
	exp_diff_status = list(filter(lambda exp: len(set(exp_info[exp]['status'])) > 1, 
							   exp_info.keys()))
	exp_diff_status_more_than_3 = list(filter(lambda exp: len(exp_info[exp]['files']) >= 4, 
										   exp_diff_status))
	
	# collect and format data that is later used to include comprehensive 
	# information within the boxplot x labels
	overview_of = exp_diff_status
	setting_plot = '\n'.join(setting[:-1])
	setting_plot += '\n# Exp: %d'%(len(overview_of))
	rel = 0
	rev = 0
	for exp in overview_of:
		for s in exp_info[exp]['status']:
			if s == 'released':
				rel += 1
			if s == 'revoked':
				rev += 1
	setting_plot += '\n# Files: %s'%(rel+rev)
	setting_plot += '\nRel: %d | Rev: %d'%(rel,rev)
	label_order.append(setting_plot)
	
	# separate feature values and class values
	X = np.array(data.iloc[:, [True if col in classification_feature_cols 
							else False for col in data.columns.values]])
	y = np.array(data['status'])
	
	# for each experiment of interest, a model is trained excluding its files
	# the model is then used to calculate a quality probability for these files
	table_str = ''
	dot_plot_table = ''
	for exp in sorted(overview_of, key=lambda x: abs(4 - len(exp_info[x]['files']))):
		train = [False if acc_exp[acc] == exp else True for acc in data['accession']]
		test  = [True if acc_exp[acc] == exp else False for acc in data['accession']]
		
		model_features = model_features.replace('FQC', 'RAW')
		model_file_path = '%s%s_%s_%s_%s.model'%(classifier_dir, organism, 
										   assay, run_type, model_features) 
		classifier = pickle.load(open(model_file_path, 'rb'))
		model = classifier.fit(X[train], y[train])
		probas = model.predict_proba(X[test])
		test_set = data.loc[test]
		sub_table = []
		for i in range(test_set.shape[0]):
			acc = test_set.iloc[i]['accession']
			status_int = test_set.iloc[i]['status']
			status = acc_status[acc]
			prob = probas[i][1]
			sub_table.append([acc, status, '%.3f'%prob, status_int])
			boxplot_data['Setting'].append(setting_plot)
			boxplot_data['prob'].append(prob)
			boxplot_data['Status'].append(status)
		exp_accession = exp.split('/')[-2]
		table_str += '%s\n'%(exp_accession)
		for row in sorted(sub_table, key=lambda x: x[-1]):
			table_str += '\t%s\t%s\t\t%s\n'%(tuple(row[:3]))
			dot_plot_table += '\t'.join([exp_accession] + row[:3]) + '\n'
	
	# save results
	table_str = table_str.replace('.', ',')
	with open('%s%s_%s_%s.txt'%(tables_dir, organism, assay, run_type) ,'w') as f:
		f.write(table_str)
	with open('%s%s_%s_%s.txt'%(dot_plot_tables_dir, organism, assay, run_type) ,'w') as f:
		f.write(dot_plot_table)
	print('')

# create and plot the boxplot
boxplot_data = pd.DataFrame(boxplot_data)
red   = '#b61a30ff'
green = '#28713cef'
colors = len(settings) * [green, red]
my_plot = sns.boxplot(x="Setting", y='prob', hue="Status",
            palette=colors,
            saturation=0.5,
            width=0.55,
            data=boxplot_data,
            order=label_order,
            hue_order = ['released', 'revoked'])
sns.despine(offset=10, trim=True)
my_plot.set_title('Probabilities of beeing Revoked\nMetrics: %s\n'%(the_features), 
				  fontsize=22)
my_plot.set_ylabel('Revoked Probability', fontsize=18)
my_plot.set_xlabel('', fontsize=16)
plt.xticks(rotation=0, fontsize=18)
plt.yticks(fontsize=16)
my_plot.legend(title = 'True Status', loc='upper right', fontsize=16, 
			   title_fontsize=18, frameon = True)
fig = my_plot.get_figure()
fig.set_size_inches(16, 14)
plt.tight_layout()
fig.savefig('./within_exp_tables/%s.png'%(the_features))




















