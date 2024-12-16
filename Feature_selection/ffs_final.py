#Progressively add non-redundant features to the regression model

import sys
import csv
import pickle
import pandas as pd
import numpy as np
from itertools import combinations
from sklearn import svm, metrics, tree, ensemble, preprocessing, pipeline, naive_bayes
from sklearn.naive_bayes import GaussianNB, MultinomialNB, ComplementNB
from sklearn.preprocessing import MinMaxScaler
from sklearn.pipeline import make_pipeline
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import AdaBoostClassifier
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import accuracy_score, balanced_accuracy_score, roc_auc_score, f1_score
from tqdm import tqdm
from multiprocessing import Pool
import itertools
from collections import Counter

data = sys.argv[1]
two_combos = sys.argv[2]
nomit = sys.argv[3]
n_feat = sys.argv[4]
pass_feat = sys.argv[5]
outpath = sys.argv[6]

n_feat = int(n_feat)
#n_feat = 5   #Arbitrary initial value
nomit = int(nomit)

df = pd.read_csv(data, sep='\t', header=0)
#print("Total no. of features = "+str(len(df.columns)-nomit-1))

feats_final = []
if(pass_feat=="all"):
	feats_final = list(df.columns)[nomit:-1]
	#print(feats_final)
else:
	with open(pass_feat) as feats:
		for feat in feats.readlines():
			feat = feat.strip()
			feats_final.append(feat)

print("Total no. of features = "+str(len(feats_final)))

with open(two_combos, 'rb') as f:
	pair_corrs = pickle.load(f)

#Using set() instead of list() to speed up search process - Proven to be blazing fast (https://stackoverflow.com/questions/5993621/fastest-way-to-search-a-list-in-python)
pair_corrs = set(pair_corrs)

def corr_function(data, combination):
	corr = data[combination[0]].corr(data[combination[1]])
	return corr

def regression_model(data, combination, y):
	X = data[list(combination)]
	#model = svm.SVC(kernel="sigmoid")
	#model = DecisionTreeClassifier()
	model = make_pipeline(MinMaxScaler(), ComplementNB())
	model.fit(X, y)
	predicted = model.predict(X)
	#pearson_corr, p_value = pearsonr(y, predicted)
	acc = accuracy_score(y, predicted)
	bal_acc = balanced_accuracy_score(y, predicted)
	roc = roc_auc_score(y, predicted)
	f1 = f1_score(y, predicted)

	Y_pred_loo = []
	Y_test_loo = []
	loo_cv = LeaveOneOut()
	for train_index, test_index in loo_cv.split(X, y):
		#print(combination)
		X_train, X_test = X.iloc[train_index], X.iloc[test_index]
		Y_train, Y_test = y.iloc[train_index], y.iloc[test_index]

		#model = DecisionTreeClassifier()
		model = make_pipeline(MinMaxScaler(), ComplementNB())
		#model = AdaBoostClassifier()
		model.fit(X_train, Y_train)
		Y_pred = model.predict(X_test)

		Y_pred_loo.append(Y_pred[0])
		Y_test_loo.append(Y_test.iloc[0])

	loo_roc = roc_auc_score(Y_test_loo, Y_pred_loo)
	loo_f1 = f1_score(Y_test_loo, Y_pred_loo)

	return(combination, acc, bal_acc, roc, f1, loo_roc, loo_f1)

def find_pass_combos(combination, pair_corrs, data, y):
	feat_pairs = set(list(combinations(combination, 2)))
	feat_pass = [True for pair in feat_pairs if(pair in pair_corrs)]
	if(len(feat_pass)==len(feat_pairs)):
		return regression_model(data, combination, y)

y = df["class"]

#-------------------------------------------------------------------------------------------------------------------------------------
models = []
for n in range(n_feat, n_feat+1):
	models = []
	feat_combinations = set(list(combinations(feats_final, n)))
	print("No. of possible "+str(n)+" feature combinations: "+str(len(feat_combinations)))

	if(n>2):
		feat_combinations_pass = []
		for combination in feat_combinations:
			feat_pairs = set(list(combinations(combination, 2)))
			feat_pass = [True for pair in feat_pairs if(pair in pair_corrs)]
			if(len(feat_pass)==len(feat_pairs)):
				feat_combinations_pass.append(combination)
	else:
		feat_combinations_pass = [combination for combination in feat_combinations if(combination in pair_corrs)]

	print("Feasible "+str(n)+" feature combinations enumerated: "+str(len(feat_combinations_pass)))

	try:
		pool = Pool(6)
		results = pool.starmap(regression_model, zip(itertools.repeat(df), feat_combinations_pass, itertools.repeat(y)))
	finally:
		pool.close()
		pool.join()

	#print(results[0:10])  #[(('total_heavy_atoms', 'v2'), 0.7635135135135135, 0.7713509316770186, 0.7713509316770186, 0.7200000000000001), (('npr1', "theta'"), 0.6216216216216216, 0.5, 0.5, 0.0), (('npr1', 'phase_angle'), 0.6216216216216216, 0.5, 0.5, 0.0)]

	models = [(list(combo), results[i][1], results[i][2], results[i][3], results[i][4], results[i][5], results[i][6]) for i, combo in enumerate(feat_combinations_pass)]
	
	with open(outpath+"Best_"+str(n)+"_feature_combos.log", 'w') as f:
		for feat_combo, acc, bal_acc, roc, f1, loo_roc, loo_f1 in models:
			#print('\t'.join(feat_combo)+"\t"+str(acc)+"\t"+str(bal_acc)+"\t"+str(roc)+"\t"+str(f1), file=f)
			print('\t'.join(feat_combo)+"\t"+str(roc)+"\t"+str(f1)+"\t"+str(loo_roc)+"\t"+str(loo_f1), file=f)

	print(str(n)+" feature combinations done.")

			
		


























