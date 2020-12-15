"""Creating quality reports based on the quality-labeled ENCODE datasets

"python guidelineReports.py --help" will display a formatted help text on 
the console. A comprehensive description is provided in the GitHub README 
that includes examples as well. In short, this script parses the given
input samples and creates an overview of several plots serving as orientation
where to position a sample with respect to different quality features.

date:	2020-12-05
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
import copy

import warnings
warnings.filterwarnings("ignore")

# import project utils
import utils.Exceptions as myExceptions
import utils.utils as utils
import utils.parser as parser

import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt

# parse command line arguments
script_dir = './'
if argv[0].find('/') >= 0:
	script_dir = argv[0][: - argv[0][::-1].find('/')]
utils_dir = '%sutils/'%(script_dir)

argsParser = argparse.ArgumentParser(description='inspectSamples - Creates a report as statistical guidelines for manual inspection.')
argsParser.add_argument('--indir', '-i', type=str, required=True, help='Input directory containing the feature set files. The feature set files are perfectly fomated by the script "deriveFeatures.py": the file names (until the ".") define the sample ID while the file endings define the corresponding feature set RAW, MAP, LOC, and TSS. By default seqQscorer applies the machine learning model to all samples from the given directory within milliseconds. However, it can be restricted to one sample using --sampleID.')
argsParser.add_argument('--species', '-s', type=str, default='generic', 
						choices=['generic','human', 'mouse'],  help='Species specifying the model used.')
argsParser.add_argument('--assay', '-a', type=str, default='generic', 
						choices=['generic','ChIP-seq','DNase-seq','RNA-seq'], help='Assay specifying the model used.')
argsParser.add_argument('--runtype', '-r', type=str, default='generic',
						choices=['generic','single-end','paired-end'], help='Run-Type specifying the model used.')
argsParser.add_argument('--outdir', '-o', type=str, default='./guideline_reports/', help='Output directory for the guideline reports. Default: "./guideline_reports/"')
argsParser.add_argument('--format', '-f', type=str, default='pdf', choices=['pdf','png', 'svg'], help='File format of the output guideline reports.')
argsParser.add_argument('--sampleID', '-id', type=str, default=None,
						help='Restrict application of seqQscorer to only one sample defined by the ID.')
args = argsParser.parse_args()

if not os.path.isdir(args.indir):
	raise myExceptions.WrongFeatureInputException(
						'"%s" is not a directory'%(args.indir))

outdir = args.outdir if args.outdir[-1] == '/' else args.outdir + '/'
if not os.path.exists(outdir):
	os.mkdir(outdir)

# load median values organized by subset, needed to impute missing values
medians = pickle.load(open('%sutils/medians.dict'%(script_dir), 'rb'))
medians = medians['generic']['generic']['generic']

# parse given input files
feature_sets = ['RAW','MAP','LOC','TSS']
input_data, feature_columns = parser.generate_input_data(args.indir, feature_sets, 'generic', medians, True, args.sampleID)

# load the generic dataset and specify it according to the user setting
data = pd.read_csv('./utils/datasets/generic_generic_generic.tsv', sep='\t')
if args.species != 'generic':
	print('Data was specified for species =', args.species)
	data = data.loc[data['organism'] == args.species]
if args.assay != 'generic':
	print('Data was specified for assay =', args.assay)
	data = data.loc[data['assay'] == args.assay]
if args.runtype != 'generic':
	print('Data was specified for runtype =', args.runtype)
	data = data.loc[data['runType'] == args.runtype]

if data.shape[0] == 0:
	raise myExceptions.WrongSettingException(
			'The data specification results in an empty ENCODE reference set.')

# settup some helper data structures and the plot settings
feature_prefix = {'RAW': 'FastQC', 'MAP': 'BowtieMI',
					'LOC': 'readsAnno', 'TSS': 'TSS'}
FastQC_value_map = {0:'FAIL', 1:'WARN', 2:'PASS'}
base_colors = ['#c51b8a','#9ecae1']
diamond_color = sns.color_palette(['#f26a24'], n_colors=1, desat=1.0)
fs_xlabel = 16
fs_ylabel = 16
diamond_size = 200
fig = plt.figure(figsize=(20,26))
title_base = 'Reference values were derived from the ENCODE dataset specified for\n'
title_base += 'species: %s, assay: %s, and run-type: %s (# samples: %d)'%(args.species, args.assay, args.runtype, data.shape[0])

for index, row in input_data.iterrows():
	sampleID = row['sampleID']
	title = title_base + '\nThe orange diamond shows the corresponding value for %s'%(sampleID)
	subplot_position = 1
	for fs in feature_sets:
		for feature in data.columns:
			if 'Basic_Statistic' in feature:
				continue
			if feature_prefix[fs] in feature:
				plt.subplot(17,2,subplot_position)
				
				# use feature name to create an appropriate x label for the subplot
				x_label = feature
				for prefix in feature_prefix.values():
					x_label = x_label.replace(prefix+'_', '')
				x_label = '%s: %s'%(fs, x_label.replace('_', ' '))
				
				if fs == 'RAW':
					ys = ['PASS', 'WARN', 'FAIL']
					
					# prepare information for stacked bar charts
					values_low = list(map(int, data.loc[data['status'] == 1][feature]))
					values_high = list(map(int, data.loc[data['status'] == 0][feature]))
					low_bars = [values_low.count(i) for i in range(3)]
					high_bars = [values_high.count(i) for i in range(3)]
					sums = [low_bars[i]+high_bars[i] for i in range(3)]
					low_bars = low_bars[::-1]
					high_bars = high_bars[::-1]
					
					# plot a stacked barchart that shows the number of low- and high-quality samples
					# for the given FastQC feature flag FAIL, WARN, or PASS
					plt.barh(ys, low_bars, left=[0,0,0], color='#b03083ff', edgecolor='white', label='Low-quality ENCODE reference',zorder=1)
					plt.barh(ys, high_bars, left=low_bars, color='#a6c7d9ff', edgecolor='white', label='High-quality ENCODE reference',zorder=1)
					
					plt.xlabel(x_label + ' (# samples)', fontsize=fs_xlabel)
					plt.yticks(ys, fontsize=12)
					
					# draw a diamond that shows the feature value of the given sample
					diamond_x = sums[row[feature]]
					plt.scatter(x=diamond_x, y=[FastQC_value_map[row[feature]]], label='Values for '+sampleID, c=diamond_color, s=diamond_size, marker='d', zorder=2)
					
					# The report title and the figure legend are attached to the two subplots on the top. 
					# These are recognized by the feature name and correspondingly, the information is added
					if 'Per_base_sequence_quality' in feature:
						plt.text(0, 4.2, s=title, fontsize = 22)
					if 'Per_tile_sequence_quality' in feature:
						plt.legend(loc=4, bbox_to_anchor=(1,1), title='', title_fontsize=18, fontsize=18)
					
				else:
					# substrackt data needed for the boxplots
					bp_data = {'Quality': ['Low' if s == 1 else 'High' for s in data['status']]}
					bp_data[x_label] = data[feature]
					bp_data['fake_y'] = data.shape[0] * ['low\nhigh']
					bp_data = pd.DataFrame(bp_data)
					
					# plot boxplot that describe the distribution of the feature values in the 
					# reference data
					ax = sns.boxplot(x=x_label, y='fake_y', hue='Quality',
								hue_order=['Low','High'],
								data=bp_data, 
								palette = base_colors,
								zorder=1,
								**{'showfliers':False})
					
					# draw a diamond that shows the feature value of the given sample
					plt.scatter(x=[row[feature]], y=[0], label='Values for '+sampleID, c=diamond_color, s=diamond_size, marker='d', zorder=2)
					
					ax.set_xlabel(x_label + ' (achieved value)', fontsize=fs_xlabel)
					ax.set_ylabel('', fontsize=1, rotation=0)
					ax.get_legend().set_visible(False)
					
					# adjust the y labels slightly
					locs, y_labels = plt.yticks(fontsize=fs_ylabel)
					for i, l in enumerate(y_labels):
						l.set_horizontalalignment('right')
						l.set_verticalalignment('center')
						l.set_position((-0.003, i))
						
				subplot_position += 1
	# finalize figure
	plt.tight_layout()
	plt.subplots_adjust(top=0.935, bottom=0.025, left=0.035, right=0.99, hspace=0.9, wspace=0.09)
	
	# save report according to the user settings (given directory and file format)
	plt.savefig('%s%s.%s'%(outdir, sampleID, args.format))
	plt.clf()
	print('\nReport created for ', sampleID)



















































































