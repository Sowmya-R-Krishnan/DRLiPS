#All with all RNA-ligand docking (blind docking) with 1000 bootstraps per RNA using NLDock (15 sites sampled per RNA)

import os
import sys
import csv
import subprocess
import pandas as pd
import numpy as np

pdb_path = sys.argv[1]
lig_path = sys.argv[2]
outpath = sys.argv[3]

for prot in os.listdir(pdb_path):
	pdb_id = prot.replace(".mol2","").replace("pdb","")
	dirpath = os.path.join(outpath, pdb_id)
	if(not os.path.exists(dirpath)):
		os.mkdir(dirpath)

	for run in range(1, 1001):
		for lig in os.listdir(lig_path):
			lig_id = lig.replace(".mol2","")

			#Run NLDock global docking for the pair
			out_var = os.system("./NLDock_v1.0/bin/NLDock "+pdb_path+prot+" "+lig_path+lig+" -out "+dirpath+"/"+pdb_id+"_"+lig_id+"_"+str(run)+" -write_max 15")
			#os.system("mv "+pdb_id+"_"+lig_id+".* "+outpath)
			print(run, pdb_id, lig_id)

	break
































