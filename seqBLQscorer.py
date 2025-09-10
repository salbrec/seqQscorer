"""Machine Learning Quality Assessment of NGS Data using the ENCODE blocklist

Main script for seqQscorer. "python seqBLQscorer.py --help" will display
a formatted help text on the console. A comprehensive description is provided
in the GitHub README that includes examples as well. In short, seqQscorer
uses dedicated machine learning algorithms to train classification models
to perform automatic NGS quality control.

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
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier

import warnings
warnings.filterwarnings("ignore")

# import project utils
import utils.Exceptions as myExceptions

SEED = 482020

# parse command line arguments
script_dir = './'
if argv[0].find('/') >= 0:
	script_dir = argv[0][: - argv[0][::-1].find('/')]
utils_dir = '%sutils/'%(script_dir)

argsParser = argparse.ArgumentParser(description='seqBLQscorer - A machine learning application for quality assessment of NGS data')
argsParser.add_argument('--indir', '-i', type=str, required=True, 
						help='Input directory containing the feature set files. The feature set files are perfectly fomated by the script "deriveBLfeatures.py": the file names (until the ".") define the sample ID while the file endings define the corresponding feature set RAW, MAP, LOC, and TSS. By default seqQscorer applies the machine learning model to all samples from the given directory within milliseconds. However, it can be restricted to one sample using --sampleID.')
argsParser.add_argument('--assembly', type=str, required=True, choices=['hg38', 'mm10'], 
						help='Species assembly defines the blocklist to be used.')
argsParser.add_argument('--assay', type=str, required=True, 
						choices=['ChIP-seq','DNase-seq','RNA-seq'], 
						help='Assay specifying the model used.')
argsParser.add_argument('--runtype', '-r', type=str, required=True,
						choices=['se','pe'], help='Run-Type specifying the model used.')
argsParser.add_argument('--mapping', '-m', type=str, default='bl', 
						choices=['BL', 'WG'], 
						help='Specify if a blocklist mapping (BL) or whole-genome (WG) mapping has been used to generate the features.')
argsParser.add_argument('--probOut', '-po', type=str, default=None,
						help='To specify an output file for the probabilities. Output will be tab-separated.')
argsParser.add_argument('--inputOut', '-io', type=str, default=None,
						help='To specify an out file that will contain the parsed input. Output will be tab-separated.')
argsParser.add_argument('--notVerbose', '-nv', action='store_true', help='Turn off verboseness, without being quiet.')
argsParser.add_argument('--sampleID', '-id', type=str, default=None,
						help='Restrict application of seqBLQscorer to only one sample defined by the ID.')

args = argsParser.parse_args()

if not os.path.isdir(args.indir):
	raise myExceptions.WrongFeatureInputException(
						'"%s" is not a directory'%(args.indir))


# define the data underlying the model
assay, assembly, runtype = args.assay, args.assembly, args.runtype
meta = pd.read_csv('%sdata/seqBLQ_meta.csv'%(utils_dir))
meta = meta.loc[ (meta['assembly'] == assembly) & (meta['assay'] == assay) & (meta['runtype'] == runtype) ]

# read in the blocklist features for the seqQscorer dataset and merge with the meta data
bl_features = pd.read_csv('%sdata/%sMap_%s.csv'%(utils_dir, args.mapping.lower(), args.assembly))
data = pd.merge(meta, bl_features, on='accession', how='inner')

# prepare data for scikit-learn model training
y = np.array(data['status'])
X = data.drop(labels=['old_accession','accession','assay','runtype','status','assembly','bl_fraqs'], axis=1)

# prepare "indir" and read in samples provided as input for the model
indir = args.indir
if indir[-1] != '/':
	indir += '/'

feature_names = list(X.columns)
samples = dict( (fn, []) for fn in feature_names )
sample_IDs = []

for subdir, dirs, files in os.walk(indir):
	for potential_feature_file in files:
		if subdir != indir:
			continue
		sample_IDs.append(potential_feature_file.split('.')[0])
		bl_values = pd.read_csv(indir + potential_feature_file)
		value_map = dict(zip( bl_values['blID'], bl_values['count'] ))
		for fn in feature_names:
			samples[fn].append( value_map.get(fn, 0.0) )
samples = pd.DataFrame(samples)

# fit PCA and determine the number of principle components (n_PCs) needed
pca = PCA(n_components=min(X.shape))
pca.fit(X)
for n_PCs in range( 2, min(X.shape) ):
	if sum(pca.explained_variance_ratio_[:n_PCs]) > 0.9999:
		break
if not args.notVerbose:
	print('\nTraining data has been prepared with %d samples and %d blocklist features.'%(X.shape))
	print('Input has been parsed with %d samples and %d blocklist features.'%(samples.shape))
	print('PCA has been applied. %d PCs will be used.'%(n_PCs))

# apply PCA dimensionality reduction and keep only "n_PCs" components
X = pca.transform(X)
X = X[:,:n_PCs]
samples = pca.transform(samples)
samples = samples[:,:n_PCs]

if not args.notVerbose:
	print('Training classifier on %d samples...'%(X.shape[0]), end=' ')
clf = RandomForestClassifier(random_state=SEED)
clf.fit(X, y)

if not args.notVerbose:
	print("Model trained!\nApplying trained model on the samples povided in '%s'"%(indir))
probas = clf.predict_proba( samples )
lowQ_probas = pd.DataFrame( {'Filename_sample_ID':sample_IDs, 'lowQ_Probability':list(probas[:,1])} )
lowQ_probas.sort_values(by='lowQ_Probability', ascending=False, inplace=True)

print('\n\nThese are the results, ranked by low-quality probability (the higher, the lower is the expected quality)\n')
print(lowQ_probas)

if args.probOut != None:
	sep = ','
	if args.probOut.endswith('.tsv') or args.probOut.endswith('.txt'):
		sep = '\t'
	lowQ_probas.to_csv(args.probOut, sep=sep, index=False)

































































































