#Remove redundancy in positive dataset by reducing each family to a maximum of 3 datapoints (mean, +SD and -SD of features)

import sys
import csv
import pandas as pd
import numpy as np

infile = sys.argv[1]
family = sys.argv[2]
outfile = sys.argv[3]
nomit = 4

input_df = pd.read_csv(infile, sep="\t", header=0)
family_df = pd.read_excel(family, header=0)

families = list(family_df["Unique_binding_regions"])
uniq_families = []
for fam in families:
	if(fam not in uniq_families):
		uniq_families.append(fam)

#print(uniq_families)
print("No. of unique families:", len(uniq_families))

out = open(outfile, "w")
print("Family_name\t"+"\t".join(list(input_df.columns)[nomit:-1]), file=out)

for j, fam in enumerate(uniq_families):
	pdb_list = list(set(family_df[family_df["Unique_binding_regions"]==fam]["id"]))

	target_df = input_df[input_df["PDB_ID"].str.contains("|".join(pdb_list))]

	if(len(list(target_df.index))>1):
		df_subset = target_df[list(input_df.columns)[nomit:-1]]
		feat_matrix = df_subset.to_numpy()

		#print(feat_matrix[0])
		avg_feats = np.mean(feat_matrix, axis=0)
		np.set_printoptions(formatter={'float_kind':'{:f}'.format})
		#print(avg_feats)

		avg_feats_final = [str(x) for x in avg_feats]

		print(fam+"\t"+"\t".join(avg_feats_final), file=out)	
		print(fam+" done.")
	else:
		feat_matrix = target_df[list(input_df.columns)[nomit:-1]].to_numpy()
		if(len(feat_matrix)==1):
			feats_final = [str(x) for x in feat_matrix[0]]
		else:
			feats_final = [str(x) for x in feat_matrix]
		print(fam+"\t"+"\t".join(feats_final), file=out)
		print(fam+" done.")
		
	#break
	



























