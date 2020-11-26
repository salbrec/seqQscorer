import numpy as np
import pandas as pd
import copy
from sklearn.model_selection import StratifiedKFold
import warnings
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	from sklearn.metrics import roc_curve, precision_recall_curve, auc, f1_score, accuracy_score, roc_auc_score, make_scorer
	from sklearn import preprocessing

import random
from scipy import interp

def auPRC_lowQual(y_true, probas):
	if sum(y_true) == len(y_true) or len(y_true) == 0:
		return -1.0
	if len(probas) == 0 or len(y_true) == 0:
		return -1.0
	precision, recall, thresholds = precision_recall_curve(y_true, probas, pos_label=1)
	return auc(recall, precision)


def F1_dt(y_true, probas, dt):
	y_pred = [1 if p > dt else 0 for p in probas]
	return f1_score(y_true, y_pred, pos_label=1)

def get_scorers():
	scorers = {'auROC': 'roc_auc'}
	scorers['auPRC'] = make_scorer(auPRC_lowQual, greater_is_better=True, needs_proba=True)
	
	return scorers

