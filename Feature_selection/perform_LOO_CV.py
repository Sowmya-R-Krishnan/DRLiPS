#Program to calculate model metrics for the best feature combination with 10-fold cross-validation

import sys
import csv
import pickle
import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler
from sklearn.pipeline import make_pipeline
from sklearn import svm, metrics, tree, ensemble, naive_bayes, neighbors, discriminant_analysis
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB, ComplementNB
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import AdaBoostClassifier, RandomForestClassifier, GradientBoostingClassifier
from sklearn.metrics import make_scorer
from sklearn.model_selection import LeaveOneOut, cross_val_score, StratifiedKFold
from sklearn.metrics import roc_auc_score, f1_score, auc

import xgboost
from xgboost import XGBClassifier

data = sys.argv[1]
best_models = sys.argv[2]
nomit = sys.argv[3]
outpath = sys.argv[4]

nomit = int(nomit)
df = pd.read_csv(data, sep='\t', header=0)
df.fillna(0, inplace=True)
#print("Total no. of features = "+str(len(df.columns)-nomit-1))

feats_final = []
feats_final = list(df.columns)[nomit:-1]
print("Total no. of features = "+str(len(feats_final)))

best_df = pd.read_csv(best_models, sep='\t', header=None)
print("No. of best models = "+str(len(best_df.index)))
print("No. of features in best models = "+str(len(best_df.columns)-4))

roc_column = len(best_df.columns)-2
roc_val = list(best_df[roc_column])
top_model_roc = max(roc_val)
top_models_df = best_df[best_df[roc_column]==top_model_roc]
top_models_df.reset_index(inplace=True, drop=True)
print(top_models_df)

#ROC scoring metric for CV
def roc_func(data1, data2):
	roc_value = roc_auc_score(data1, data2, max_fpr=0.2)  #Setting max_fpr to compute partial AUC score
	return roc_value

#F1 scoring metric for CV
def f1_func(data1, data2):
	f1_value = f1_score(data1, data2)
	return f1_value
	
#AUC scoring metric for CV
def auc_func(data1, data2):
	auc_value = auc(data1, data2)
	return auc_value

for i, row in top_models_df.iterrows():
	top_model_specs = row.tolist()

	top_model_feat = top_model_specs[:-4]
	print("Top model "+str(i)+" features: "+','.join(top_model_feat))

	#Creating a dataset out of best model features for cross-validation
	X_with_label = df['Family_name'].to_numpy()
	X_dataset = df[top_model_feat].to_numpy()
	Y_dataset = df['class'].to_numpy()

	#model = svm.SVC()
	#model = GaussianNB()
	#model = RandomForestClassifier()
	#model = DecisionTreeClassifier()
	#model = KNeighborsClassifier()
	#model = AdaBoostClassifier()
	#model = LinearDiscriminantAnalysis()
	#model = GradientBoostingClassifier()
	#model = XGBClassifier()
	model = make_pipeline(MinMaxScaler(), svm.SVC(kernel="sigmoid"))

	model.fit(X_dataset, Y_dataset)
	Y_pred_train = model.predict(X_dataset)
	print("Training partial ROC: "+str(roc_auc_score(Y_dataset, Y_pred_train, max_fpr=0.2)))
	print("Training F1 score: "+str(f1_score(Y_dataset, Y_pred_train)))

	roc_scorer = make_scorer(roc_func)
	f1_scorer = make_scorer(f1_func)

	cv_roc = cross_val_score(model, X_dataset, Y_dataset, cv=10, scoring=roc_scorer)
	cv_f1 = cross_val_score(model, X_dataset, Y_dataset, cv=10, scoring=f1_scorer)
	print("10 fold CV results")
	print("-------------------")
	print("10-fold CV partial ROC: "+str(np.mean(cv_roc)))
	print("10-fold CV F1: "+str(np.mean(cv_f1)))

	strat_cv_roc = []
	strat_cv_f1 = []
	strat_cv_auc = []

	strat_cv_folds = StratifiedKFold(n_splits=10, shuffle=True)
	strat_cv_folds.get_n_splits(X_dataset, Y_dataset)
	cv_epoch = 1

	for j, (train_index, test_index) in enumerate(strat_cv_folds.split(X_dataset, Y_dataset)):
		X_train, X_test = X_dataset[train_index], X_dataset[test_index]
		Y_train, Y_test = Y_dataset[train_index], Y_dataset[test_index]

		model = make_pipeline(MinMaxScaler(), svm.SVC(kernel="sigmoid"))
		model.fit(X_train, Y_train)
		Y_pred = model.predict(X_test)
		
		roc = roc_auc_score(Y_test, Y_pred, max_fpr=0.2)
		f1 = f1_score(Y_test, Y_pred)

		strat_cv_roc.append(roc)
		strat_cv_f1.append(f1)
		cv_epoch = cv_epoch + 1

	print(" ")
	print("Stratified 10 fold CV results")
	print("-------------------")
	print("Average partial ROC score: "+str(np.mean(strat_cv_roc)))
	print("Average F1 score: "+str(np.mean(strat_cv_f1)))

	roc_val = []
	f1_val = []

	Y_pred_loo = []
	Y_test_loo = []

	loo_cv = LeaveOneOut()
	cv_epoch = 1
	for train_index, test_index in loo_cv.split(X_dataset, Y_dataset):
		X_train, X_test = X_dataset[train_index], X_dataset[test_index]
		Y_train, Y_test = Y_dataset[train_index], Y_dataset[test_index]

		#model = svm.SVC()
		#model = GaussianNB()
		#model = RandomForestClassifier()
		#model = DecisionTreeClassifier()
		#model = KNeighborsClassifier()
		#model = AdaBoostClassifier()
		#model = LinearDiscriminantAnalysis()
		#model = GradientBoostingClassifier()
		#model = XGBClassifier()
		model = make_pipeline(MinMaxScaler(), svm.SVC(kernel="sigmoid"))
		model.fit(X_train, Y_train)
		Y_pred = model.predict(X_test)

		Y_pred_loo.append(Y_pred[0])
		Y_test_loo.append(Y_test[0])
		cv_epoch = cv_epoch + 1

	loo_roc = roc_auc_score(Y_test_loo, Y_pred_loo, max_fpr=0.2)
	loo_f1 = f1_score(Y_test_loo, Y_pred_loo)

	print(" ")
	print("LOO CV results")
	print("-------------------")
	print("Average partial ROC score: "+str(loo_roc))
	print("Average F1 score: "+str(loo_f1))
	print(" ")

	break



















