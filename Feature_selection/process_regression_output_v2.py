#Program to compare every successive set of feature combination results and extract the best models

import sys
import csv
import pandas as pd

data1 = sys.argv[1]
data2 = sys.argv[2]
least_featnum = sys.argv[3]
outfile = sys.argv[4]

print(data1, data2, least_featnum)

least_featnum = int(least_featnum)

df1 = pd.read_csv(data1, sep='\t', header=None)
df2 = pd.read_csv(data2, sep='\t', header=None)

print(df2.columns)
print(df2.head())

#min_model_metric = max(list(df1.iloc[:,-1:]))
min_model_metric_1 = max(list(df1[least_featnum+2]))
min_model_metric_2 = max(list(df1[least_featnum+3]))
print(min_model_metric_1, min_model_metric_2)
best_combo = df1.loc[df1[least_featnum+2]==min_model_metric_1]
best_combo = best_combo.loc[best_combo[least_featnum+3]==min_model_metric_2]

#print(min_model_metric)
print("Intial best model:", len(best_combo.index))
print(best_combo)
print("")

best_models = df2.loc[df2[least_featnum+3]>min_model_metric_1]
best_models = best_models.loc[best_models[least_featnum+4]>min_model_metric_2]
print("Final best models:", len(best_models.index))
print(best_models)

best_models.to_csv(outfile, sep='\t', index=False, header=False)





























