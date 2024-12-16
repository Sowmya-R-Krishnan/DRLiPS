#Program to automate fpocket

import os
import sys
import csv
import re
import math
import pandas as pd
import numpy as np

pdb_list = sys.argv[1]
inpath = sys.argv[2]
pocket_info = sys.argv[3]
outpath = sys.argv[4]

pdb_df = list(pd.read_csv(pdb_list, sep="\t", header=None)[0])
pocket_df = pd.read_excel(pocket_info, header=0)

print("No. of structures to be checked:", len(pdb_df))

for i, pdb in enumerate(pdb_df):
	fname = pdb.replace(".pdb","")
	pdb_id = fname.replace("pdb","")

	sub_df = pocket_df[pocket_df["PDB_ID"]==pdb_id]
	
	for i, row in sub_df.iterrows():
		if(int(row["Model_ID"])>0):
			os.system("fpocket -f "+inpath+pdb+" -l "+str(int(row["Model_ID"])))
			os.system("mv "+inpath+fname+"_out/ "+outpath+fname+"_"+str(int(row["Model_ID"]))+"_out/")
		else:
			os.system("fpocket -f "+inpath+pdb)
			os.system("mv "+inpath+fname+"_out/ "+outpath+fname+"_"+str(int(row["Model_ID"]))+"_out/")

	#break
	



































