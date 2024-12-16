#Program to combine the positive and negative dataset to create the final pocket dataset for training the models

import sys
import csv
import pandas as pd
import numpy as np

data1 = sys.argv[1]
data2 = sys.argv[2]
outfile = sys.argv[3]

df1 = pd.read_csv(data1, sep="\t", header=0, index_col=False)
df2 = pd.read_csv(data2, sep="\t", header=0, index_col=False)

#print(df1.columns)
#print(df2.columns)

#columns_selected = list(df1.columns)[1:]
columns_selected = list(df1.columns)

final_set = []
for i, row in df1.iterrows():
	df_dict = row.to_dict()
	df_dict["class"] = 1
	final_set.append(df_dict)

for j, row in df2.iterrows():
	df_dict = row.to_dict()
	df_dict["Family_name"] = "Neg_"+str(j+1)

	#if(df_dict['Family_name']=="Neg_19"):
		#print(df_dict)

	dict_subset = {k: df_dict[k] for k in columns_selected}
	dict_subset["class"] = 0
	final_set.append(dict_subset)

df_final = pd.DataFrame(final_set)
print("Final dataset size:", len(final_set))

df_final.to_csv(outfile, sep="\t", header=True, index=False)
print("Dataset saved.")


























