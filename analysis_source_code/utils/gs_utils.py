"""Grid search utils

Util functions used by the grid search script

Methods
-------

get_abbr_maps()
	provides dictionaries to map between short and long names 
	that describe the different feature sets

date:	2019-02-10
author:	Steffen Albrecht

"""

def get_abbr_maps():
	feature_abbr_ls = {}
	feature_abbr_ls['FastQC'] = 'FQC'
	feature_abbr_ls['BowtieSE'] = 'BowS'
	feature_abbr_ls['BowtiePE'] = 'BowP'
	feature_abbr_ls['BowtieMI'] = 'BowM'
	feature_abbr_ls['readsAnno'] = 'RA'
	feature_abbr_ls['TSS'] = 'TSS'
	
	feature_abbr_sl = {}
	feature_abbr_sl['FQC'] = 'FastQC'
	feature_abbr_sl['BowS'] = 'BowtieSE'
	feature_abbr_sl['BowP'] = 'BowtiePE'
	feature_abbr_sl['BowM'] = 'BowtieMI'
	feature_abbr_sl['Bow'] = 'Bowtie'
	feature_abbr_sl['RA'] = 'readsAnno'
	feature_abbr_sl['TSS'] = 'TSS'
	
	return feature_abbr_ls, feature_abbr_sl



	