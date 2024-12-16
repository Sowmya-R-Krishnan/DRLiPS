#Program to test the best models on the external test dataset collected

import sys
import csv
import pickle
import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler
from sklearn.pipeline import make_pipeline
from sklearn import svm, metrics
from sklearn.metrics import roc_auc_score, f1_score, balanced_accuracy_score, accuracy_score

test_data = sys.argv[1]
dataset = sys.argv[2]
best_models = sys.argv[3]
outpath = sys.argv[4]
nomit = sys.argv[5]

nomit = int(nomit)
df = pd.read_csv(dataset, sep='\t', header=0)
df.fillna(0, inplace=True)

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

external_df = pd.read_csv(test_data, sep="\t", header=0)

for i, row in top_models_df.iterrows():
	top_model_specs = row.tolist()
	top_model_feat = top_model_specs[:-4]
	print("Top model "+str(i)+":", top_model_specs)

	X_dataset = df[top_model_feat].to_numpy()
	Y_dataset = df['class'].to_numpy()
	model = make_pipeline(MinMaxScaler(), svm.SVC(kernel="sigmoid", probability=True))
	model.fit(X_dataset, Y_dataset)

	test_X = external_df[top_model_feat].to_numpy()
	test_Y = external_df['class'].to_numpy()
	test_labels = list(external_df["PDB_ID"])
	Y_pred = model.predict(test_X)

	roc_test = roc_auc_score(test_Y, Y_pred, max_fpr=0.2)
	f1_test = f1_score(test_Y, Y_pred)
	acc_test = accuracy_score(test_Y, Y_pred)
	bal_acc_test = balanced_accuracy_score(test_Y, Y_pred)

	print("External test set ROC: "+str(roc_test))
	print("External test set F1: "+str(f1_test))
	print("External test set accuracy: "+str(acc_test))
	print("External test set balanced accuracy: "+str(bal_acc_test))
	print("External test set balanced accuracy: "+str(bal_acc_test))





























