"""Cross-species generalization

Slightly differs to the main grid search. For this grid, XG-Boost 
was excluded because of problems with the serialization

Methods
-------

get_grid()
	defines the algorithm list and all parameters used within grid search

date:	2019-04-08
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

import numpy as np

def get_grid():
	SEED = 1
	the_RS = np.random.RandomState(SEED)
	
	classifiers = []
	algorithms = {}
	param_grid = {}
	param_grid_test = {}

	# Random Forest
	classifiers.append('RFC')
	algorithms['RFC'] = RandomForestClassifier(random_state=the_RS)
	param_grid['RFC'] = {'n_estimators': [100, 200, 500, 1000],
						'max_features': [None, 'auto', 'sqrt', 'log2'],
						'max_depth' : [None, 5,10,20],
						'criterion' : ['gini', 'entropy'],
						'warm_start': [False,True], 
						'oob_score' : [False, True]}
	param_grid_test['RFC'] = {'n_estimators': [50], 
							'max_depth': [4,5],
							'criterion' : ['gini'],
							'warm_start': [False], 
							'oob_score' : [False]}

	# Gradient Boosting Classifier
	classifiers.append('GBC')
	algorithms['GBC'] = GradientBoostingClassifier(random_state=the_RS)
	param_grid['GBC'] = {"loss":["deviance"],
						"learning_rate": [0.01, 0.05, 0.1, 0.2],
						"min_samples_split": np.linspace(0.1, 0.5, 5),
						"min_samples_leaf": np.linspace(0.1, 0.5, 5),
						"max_depth":[3,5,8],
						"max_features":["log2","sqrt"],
						"criterion": ["friedman_mse",  "mae"],
						"subsample":[0.5, 0.618, 0.85, 0.9, 1.0],
						"n_estimators":[50, 100, 200]}
	param_grid_test['GBC'] = {'n_estimators': [100], 
							"criterion": ["friedman_mse",  "mae"],
							"max_depth":[3,5,8],
							"max_features":["log2","sqrt"]}

	# Logistic Regression
	classifiers.append('LRN')
	algorithms['LRN'] = LogisticRegression(random_state=the_RS)
	param_grid['LRN'] = grid={"C":np.logspace(-3,3,7), 
							"penalty":["l1","l2"]}
	param_grid_test['LRN'] = grid={"C":[3], 
								"penalty":["l1"]}

	# Support Vector Machine
	classifiers.append('SVC')
	algorithms['SVC'] = SVC(kernel='rbf', probability=True, random_state=the_RS)
	param_grid['SVC'] = {'C':[0.001, 0.01, 0.1, 1, 10], 
						'gamma': [0.001, 0.01, 0.1, 1]}
	param_grid_test['SVC'] = {'C':[1], 
							'gamma': [0.001]}

	# KNeighborsClassifier
	classifiers.append('KNN')
	algorithms['KNN'] = KNeighborsClassifier()
	param_grid['KNN'] = {'n_neighbors': [5, 6, 7, 8, 9],
						'metric': ['minkowski', 'euclidean', 'manhattan'], 
						'weights': ['uniform', 'distance']}
	param_grid_test['KNN'] = {'n_neighbors': [3], 
							'weights': ['distance']}

	# GaussianNB
	classifiers.append('GNB')
	algorithms['GNB'] = GaussianNB()
	param_grid['GNB'] = {}
	param_grid_test['GNB'] = {}

	# MLPClassifier
	classifiers.append('MLP')
	algorithms['MLP'] = MLPClassifier(random_state=the_RS)
	param_grid['MLP'] = {'solver': ['lbfgs', 'sgd'], 
						'max_iter': [100,500,1000,1500,3000], 
						'alpha': 10.0 ** -np.arange(1, 7), 
						'hidden_layer_sizes':[(5,) , (10,) , (50,), (100,) , (5,5), (10,10), (50,10), (100,50), (200,100), (500,400)]}
	param_grid_test['MLP'] = {'solver': ['lbfgs'],
							'max_iter': [50]}

	# AdaBoost with DTL
	classifiers.append('ADT')
	DTC = DecisionTreeClassifier(random_state = the_RS, max_features = "auto", class_weight = "balanced",max_depth = None)
	algorithms['ADT'] = AdaBoostClassifier(base_estimator = DTC)
	param_grid['ADT'] = {'base_estimator__criterion':["gini", "entropy"],
						'base_estimator__splitter':["best", "random"],
						'n_estimators':[1, 2, 10]}
	param_grid_test['ADT'] = {'base_estimator__criterion':["entropy"],
							'n_estimators':[1]}

	# Extra Tree Classifier
	classifiers.append('ETC')
	algorithms['ETC'] = ExtraTreeClassifier(random_state=the_RS)
	param_grid['ETC'] = {'criterion':['gini', 'entropy'],
						'splitter':['best', 'random'],
						'max_depth': [None, 1, 2, 3, 4, 5],
						'max_features': [None, 1, 2, 3, 4]}
	param_grid_test['ETC'] = {'max_depth': [2],
					 'max_features': [2]}
	
	return classifiers, algorithms, param_grid, param_grid_test

