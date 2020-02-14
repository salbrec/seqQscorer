"""Machine learning utils

This script provides functions that calculate certain measures used to 
validate classification models. 

Methods
-------

get_auPRC_released(real_y, probas)
	Area under precision recall curve for class label "released"
get_auPRC_revoked(real_y, probas)
	Area under precision recall curve for class label "revoked"
my_F1(precision, recall)
	F1 measure
topF1_Released(real_y, probas)
	Searches for the decision threshold that leads to the highest F1 
	F1 measure is returned and calculated for class label "released"
topF1_TH_Released(real_y, probas)
	Searches for the decision threshold that leads to the highest F1 
	decision threshold is returned, F1 measure is calculated for 
	class label "released"
topF1_Revoked(real_y, probas)
	Searches for the decision threshold that leads to the highest F1 
	F1 measure is returned and calculated for class label "revoked"
topF1_TH_Revoked(real_y, probas)
	Searches for the decision threshold that leads to the highest F1 
	decision threshold is returned, F1 measure is calculated for 
	class label "revoked"
F1_released(real_y, probas)
	returns F1 measure for class label "released"
F1_revoked(real_y, probas)
	returns F1 measure for class label "revoked"
get_scorer_dict()
	initializes and returns a dictionary with all validation 
	measures (scores) used within the grid search to validate 
	the classification models

date:	2019-02-10
author:	Steffen Albrecht

"""

import numpy as np
import pandas as pd
import copy
from sklearn.model_selection import StratifiedKFold
import warnings
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	from sklearn.metrics import roc_curve, precision_recall_curve, auc, f1_score, accuracy_score, roc_auc_score, make_scorer
	from sklearn import preprocessing

def get_auPRC_released(real_y, probas):
	if sum(real_y) == len(real_y) or len(real_y) == 0:
		return -1.0
	if len(probas) == 0 or len(real_y) == 0:
		return -1.0
	precision, recall, thresholds = precision_recall_curve(real_y, 1.0-probas, pos_label=0)
	return auc(recall, precision)

def get_auPRC_revoked(real_y, probas):
	if sum(real_y) == len(real_y) or len(real_y) == 0:
		return -1.0
	if len(probas) == 0 or len(real_y) == 0:
		return -1.0
	precision, recall, thresholds = precision_recall_curve(real_y, probas, pos_label=1)
	return auc(recall, precision)

def my_F1(precision, recall):
	return 2.0 * (precision * recall) / (precision + recall)

def topF1_Released(real_y, probas):
	if sum(real_y) == len(real_y) or len(real_y) == 0:
		return -1.0
	if len(probas) == 0 or len(real_y) == 0:
		return -1.0
	precision, recall, thresholds = precision_recall_curve(real_y, 1.0-probas, pos_label=0)
	threshold_f1 = []
	for i in range(len(thresholds)):
		threshold_f1.append((thresholds[i], my_F1(precision[i], recall[i])))
	threshold_f1 = sorted(threshold_f1, key = lambda x: x[1], reverse=True)
	return threshold_f1[0][1]

def topF1_TH_Released(real_y, probas):
	if sum(real_y) == len(real_y) or len(real_y) == 0:
		return -1.0
	if len(probas) == 0 or len(real_y) == 0:
		return -1.0
	precision, recall, thresholds = precision_recall_curve(real_y, 1.0-probas, pos_label=0)
	threshold_f1 = []
	for i in range(len(thresholds)):
		threshold_f1.append((thresholds[i], my_F1(precision[i], recall[i])))
	threshold_f1 = sorted(threshold_f1, key = lambda x: x[1], reverse=True)
	return threshold_f1[0][0]

def topF1_Revoked(real_y, probas):
	if sum(real_y) == len(real_y) or len(real_y) == 0:
		return -1.0
	if len(probas) == 0 or len(real_y) == 0:
		return -1.0
	precision, recall, thresholds = precision_recall_curve(real_y, probas, pos_label=1)
	threshold_f1 = []
	for i in range(len(thresholds)):
		threshold_f1.append((thresholds[i], my_F1(precision[i], recall[i])))
	threshold_f1 = sorted(threshold_f1, key = lambda x: x[1], reverse=True)
	return(threshold_f1[0][1])

def topF1_TH_Revoked(real_y, probas):
	if sum(real_y) == len(real_y) or len(real_y) == 0:
		return -1.0
	if len(probas) == 0 or len(real_y) == 0:
		return -1.0
	precision, recall, thresholds = precision_recall_curve(real_y, probas, pos_label=1)
	threshold_f1 = []
	for i in range(len(thresholds)):
		threshold_f1.append((thresholds[i], my_F1(precision[i], recall[i])))
	threshold_f1 = sorted(threshold_f1, key = lambda x: x[1], reverse=True)
	return(threshold_f1[0][0])

def F1_released(real_y, probas):
	return f1_score(real_y, [1 if p > 0.5 else 0 for p in probas], pos_label=0)

def F1_revoked(real_y, probas):
	return f1_score(real_y, [1 if p > 0.5 else 0 for p in probas], pos_label=1)

def topACC(real_y, probas):
	accs = []
	for threshold in [float(i) / 100.0 for i in range(0, 101)]:
		predicted = [1 if p > threshold else 0 for p in probas]
		correct = [True if predicted[i] == real_y[i] else False for i in range(len(predicted))]
		accs.append(float(sum(correct)) / float(len(correct)))
	return max(correct)

def get_scorer_dict():
	scorer = {'auROC': 'roc_auc', 'ACC': 'accuracy', 
		   'Prec': 'precision', 'Recall': 'recall'}
	scorer['F1_rev'] = make_scorer(F1_revoked, greater_is_better=True, needs_proba=True)
	scorer['F1_def'] = 'f1'
	scorer['F1_rel'] = make_scorer(F1_released, greater_is_better=True, needs_proba=True)
	
	scorer['auPRC_Rev'] = make_scorer(get_auPRC_revoked, greater_is_better=True, needs_proba=True)
	scorer['auPRC_Rel'] = make_scorer(get_auPRC_released, greater_is_better=True, needs_proba=True)
	
	scorer['topF1_Rev'] = make_scorer(topF1_Revoked, greater_is_better=True, needs_proba=True)
	scorer['topF1_Rel'] = make_scorer(topF1_Released, greater_is_better=True, needs_proba=True)
	
	scorer['topF1_TH_Rev'] = make_scorer(topF1_TH_Revoked, greater_is_better=True, needs_proba=True)
	scorer['topF1_TH_Rel'] = make_scorer(topF1_TH_Released, greater_is_better=True, needs_proba=True)
	
	scorer['topACC'] = make_scorer(topACC, greater_is_better=True, needs_proba=True)
	
	return scorer

