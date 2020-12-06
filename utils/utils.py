"""seqQscorer utils

Util functions used by seqQscorer.

Methods
-------

get_path_info()
	given a file path, this function extracts the files directory, the 
	base file name, and its ending
get_help_text()
	returns a formatted string describing the seqQscorer help text
clf_full_names(abbr)
	given one of the abbreviations, this function returns the full name of
	the algorithm. Used to clarify the terminal output
get_best_classifier(utils_dir, species, assay, run_type, feature_sets, fs_suffix, metric)
	given the user specifications from seqQscorer, this function parses a text table
	in order to return the classifier and feature selection specifications that are most 
	recommendable for the application
read_in_measure_table(utils_dir, species, assay, run_type, feature_sets, fs_suffix, metric)
	seqQscorer prints a table with machine learning evaluation measures for different decision 
	thresholds. The source file with this information is parsed by this function
def get_clf_algos()
	this function creates and returns a dictionary of default classifier configuratons

date:	2020-11-10
author:	Steffen Albrecht

"""

from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.neural_network import MLPClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import AdaBoostClassifier
from sklearn.tree import ExtraTreeClassifier

import json
from terminaltables import AsciiTable

def print_nice_table(table):
	print_table = AsciiTable(table)
	print(print_table.table)

def get_path_info(file_path):
	if file_path.find('/') < 0:
		folder = './'
		file_name_ending = file_path
	else:
		folder = file_path[:-file_path[::-1].find('/')]
		file_name_ending = file_path[-file_path[::-1].find('/'):]
	file_name = file_name_ending[:file_name_ending.find('.')]
	return folder, file_name_ending, file_name

def clf_full_names(abbr):
	name_map = {'RFC':'Random Forest', 'MLP':'Multi-Layer Perceptron', 'GBC':'Gradient Boosting', 
			 'LRN':'Logistic Regression', 'SVC':'Support Vector Machine', 'GNB':'Naive Bayes', 
			 'ADT':'Adaboost with Decision Tree', 'ETC':'Randomized Decision Tree', 
			 'KNN':'K-nearest Neighbor'}
	return name_map[abbr]

def get_best_classifier(utils_dir, species, assay, run_type, feature_sets, fs_suffix, metric):
	case = '%s_%s_%s_%s'%(species, assay, run_type, '-'.join(feature_sets))
	
	classifier = None
	feature_selection = None
	selection_str = ''
	selection = None
	params = None
	auROC = None
	brier = None
	table_file_path = '%stables/best_algo_params_%s%s.tsv'%(utils_dir, metric, fs_suffix)
	with open(table_file_path, 'r') as f:
		for line in f:
			line = line.strip().split('\t')
			if line[0] == case:
				classifier = line[1]
				feature_selection = line[2]
				selection_str = line[3]
				params = json.loads(line[4])
				auROC = line[5]
				brier = line[6]
				break
	
	if feature_selection != 'No-100':
		selection = [True if b == '1' else False for b in selection_str.split(',')]
	
	return classifier, feature_selection, selection, params, auROC, brier

def read_in_measure_table(utils_dir, species, assay, run_type, feature_sets, fs_suffix, metric):
	setting = '%s_%s_%s'%(species, assay, run_type)
	
	table_str = None
	table_file_path = '%stables/dec_threshold_tabs_%s%s.tsv'%(utils_dir, metric, fs_suffix)
	with open(table_file_path, 'r') as f:
		for line in f:
			line = line.strip().split('\t')
			if line[0] == setting and line[1] == '-'.join(feature_sets):
				table_str = line[2]
				break
	table = []
	for row in table_str.split('|'):
		table.append( list(row.split(',')) )
	return table

def get_clf_algos():
	algorithms = {}
	algorithms['RFC'] = RandomForestClassifier(random_state=1)
	algorithms['GBC'] = GradientBoostingClassifier(random_state=1)
	algorithms['LRN'] = LogisticRegression(random_state=1)
	algorithms['SVC'] = SVC(kernel='rbf', probability=True, random_state=1)
	algorithms['GNB'] = GaussianNB()
	algorithms['MLP'] = MLPClassifier(random_state=1)
	DTC = DecisionTreeClassifier(random_state = 1, max_features = "auto", class_weight = "balanced",max_depth = None)
	algorithms['ADT'] = AdaBoostClassifier(base_estimator = DTC)
	algorithms['ETC'] = ExtraTreeClassifier(random_state=1)
	return algorithms


























