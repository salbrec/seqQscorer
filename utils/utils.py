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

def getFileName(file_path):
	file_name = file_path[ -file_path[::-1].find('/') : ]
	if file_name.endswith('.gz'):
		file_name = file_name[:-3]
	if file_name.endswith('.fastq'):
		file_name = file_name[:-6]
	elif file_name.endswith('.fq'):
		file_name = file_name[:-3]
	return file_name

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

def read_Bowtie_stats(fp):
    """
    Reads a Bowtie2 report and returns the feature values in a dict.

    Args:
        fp (str): File path to FastQC summary.txt.

    Returns:
        dict: a map of MAP feature names and feature values.
        str: unparsed file content or "not exist" message.
    """
    if not os.path.exists(fp):
        return None, 'File does not exist.'
    lines = open(fp,'r').read().split('\n')
    if len(lines) != 7:
        return None, lines
    stats = {'total': int(lines[0].split()[0]) }
    stats['unpaired'] = int(lines[1].split()[0])
    stats['0times'] = int(lines[2].split()[0])
    stats['1time'] = int(lines[3].split()[0])
    stats['multi'] = int(lines[4].split()[0])
    stats['overall'] = stats['1time'] + stats['multi']
    percentages = {}
    for key in stats:
        if key != 'total':
            percentages['perc_'+key] = stats[key] / stats['total'] * 100.0
    stats.update(percentages)
    return stats, lines

def get_file_length(fp):
    """
    Use the linux function to get the number of lines in a file.

    Args:
        fp (str): File path

    Returns:
        int: number of lines
    """
    wc = None
    try:
        call = 'wc -l %s'%(fp)
        res = subprocess.check_output(call, shell=True, text=True)
        wc = int(res.split()[0])
    except:
        print('Could not get the wc -l here:',fp)
    return wc

def read_blocklist(bl_file):
    """
    Reads in a blocklist file and uses a certain format that allows for 
    using the regions efficiently in read counting procedure

    Args:
        bl_file (str): File path to the blocklist regions file.

    Returns:
        dict: dictionary containing lists of regions for each chromosome used as key
    """
    df_blocklist = pd.read_csv(bl_file, sep='\t', names=['chr','start','end','ID'])
    blocklist = {}
    for index, row in df_blocklist.iterrows():
        if not row['chr'] in blocklist:
            blocklist[row['chr']] = []
        bl_type = row['ID'].split('_')[0][2:]
        blocklist[row['chr']].append( (index+1, row['chr'], row['start'], row['end'],
            bl_type_map[bl_type], row['ID']) )
    return blocklist

def count_reads_in_regions(summits, regions, chrom_size_map, bl_mapping):
    """
    For each region, counting the summits within the region. 
    An overlap is only given if the summit is in the region, hence, 
    if more than half of the read overlaps with the blocklist region
    Args:
        summits (dict): dictionary with chromosome as key. The values
                        are lists of summits for the chromosome. The summits 
                        describe the center of the reads in this application.
        regions (dict): dictionary with chromosome as key. The values 
                        are lists of blocklist regions. 
        chrom_size_map (dict): dictionary with chromosome as key. The values
                                describe the chromosome length.
    Returns:
        pd.DataFrame: 
    """
    bincov = {'binID':[], 'chr':[], 'start':[], 'end':[], 'count':[],
        'blID':[], 'blType':[]}

    if bl_mapping:
        for chrom in chrom_size_map:
            if not chrom in regions:
                continue
            for binID, reg_chrom, reg_start, reg_end, reg_type, reg_ID in regions[chrom]:
                if reg_ID in summits:
                    bincov['binID'].append( binID ); bincov['chr'].append( reg_chrom );
                    bincov['start'].append( reg_start ); bincov['end'].append( reg_end );
                    bincov['count'].append( len(summits[reg_ID]) ); 
                    bincov['blType'].append( reg_type ); bincov['blID'].append( reg_ID );
        return pd.DataFrame(bincov)
    else:
        for chrom in chrom_size_map:
            if not chrom in summits or not chrom in regions:
                continue

            for binID, reg_chrom, reg_start, reg_end, reg_type, reg_ID in regions[chrom]:
                if chrom != reg_chrom:
                    print('!!! Something is WRONG here: ', reg_chrom, chrom)
                    return None
            
                count = len(list(filter( lambda x: x > reg_start and x < reg_end, summits[chrom] )))
            
                if count != 0:
                    bincov['binID'].append( binID ); bincov['chr'].append( reg_chrom );
                    bincov['start'].append( reg_start ); bincov['end'].append( reg_end );
                    bincov['count'].append( count ); bincov['blType'].append( reg_type );
                    bincov['blID'].append( reg_ID )
        return pd.DataFrame(bincov)























