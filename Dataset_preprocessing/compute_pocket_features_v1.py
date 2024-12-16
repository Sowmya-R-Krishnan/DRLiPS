#Program to compute binding site features for a given dataset of binding sites

import os
import sys
import csv
import re
import pandas as pd
import numpy as np
import sympy as sy
import Bio
import warnings
import math
import ResidueDepth_Mod
import subprocess
from Bio import PDB
from Bio.PDB import PDBParser, PDBIO, SASA, ResidueDepth
from Bio.PDB.ResidueDepth import get_surface, residue_depth
from collections import defaultdict, Counter

infile = sys.argv[1]
filetype = sys.argv[2]
pdb_path = sys.argv[3]
msms_path = sys.argv[4]
naccess_path = sys.argv[5]
dssr_path = sys.argv[6]
outpath = sys.argv[7]

if(filetype=="excel"):
	input_df = pd.read_excel(infile, header=0)
elif(filetype=="tsv"):
	input_df = pd.read_csv(infile, sep="\t", header=0)
elif(filetype=="csv"):
	input_df = pd.read_csv(infile, sep=",", header=0)
else:
	print("Only excel, tsv and csv formats supported currently. Please ensure that the input format is one of these.")
	sys.exit()

#--------------------------------------------------------ALL UTILITY FUNCTIONS-------------------------------------------------------------
#XXX: All PDB format definitions for Python taken from: https://cupnet.net/pdb-format/
#Function to extract the atoms corresponding to every model deposited in an NMR structure
def extract_modelwise_atoms(pdb_id):
	model_starts = []
	lines_all = []
	line_index = 0
	with open(pdb_path+pdb_id+".pdb") as f:
		for line in f.readlines():
			line = line.strip()
			lines_all.append(line)
			if(line.startswith("MODEL")):
				model_starts.append(line_index)

			line_index = line_index + 1

	model_dict = {}
	for i, idx in enumerate(model_starts):
		final_atoms = []
		try:
			model_range = (idx, model_starts[i+1]-1)
			lines = lines_all[model_range[0]:model_range[1]]
			final_atoms = [x for x in lines if(x.startswith("ATOM"))]
			model_dict[i] = final_atoms
		except:
			model_range = (idx, len(lines_all))
			lines = lines_all[model_range[0]:model_range[1]]
			final_atoms = [x for x in lines if(x.startswith("ATOM"))]
			model_dict[i] = final_atoms

	return model_dict

#Function to extract atoms corresponding to a residue ID
def extract_atoms_from_resid(pdb_atoms, residue_id):
	atoms = []
	for atom in pdb_atoms:
		chain = atom[21:22]
		residue = atom[22:27].strip()
		resiname = atom[17:20].strip()
		resid = str(resiname)+"_"+str(residue)+"_"+str(chain)

		if(resid==residue_id):
			atoms.append(atom)

	return atoms

#Function to extract the atoms corresponding to the binding site of interest
def extract_binding_site(pdb_atoms, residue_list):
	bs_atoms = []
	for res_id in residue_list:
		res_atoms = extract_atoms_from_resid(pdb_atoms, res_id)
		bs_atoms.extend(res_atoms)
	return bs_atoms

#Function to extract the PATTY atom type mapping for a list of PDB atoms
def extract_patty_types(pdb_atoms, patty_types):
	patty_list = []

	for atom in pdb_atoms:
		resiname = atom[17:20].strip()
		if(resiname in ["A", "C", "G", "U"]):
			atom_name = atom[12:16].strip()
			#To handle hydrogens present in crystal structures
			try:
				patty_type = patty_types[resiname][atom_name]
				patty_list.append(patty_type)
			except:
				continue

	return patty_list

#Function to translate a molecule to its centre of mass so that the principal moments are comparable
def translate_to_cofm(masses, xyz):
	# Position of centre of mass in original coordinates
	cofm = sum(masses[:,np.newaxis] * xyz) / np.sum(masses)
	# Transform to CofM coordinates and return
	xyz -= cofm
	return xyz

#Function to get the inertia matrix from masses and coordinates
def get_inertia_matrix(masses, xyz):
	# Moment of intertia tensor
	xyz = translate_to_cofm(masses, xyz)
	x, y, z = xyz.T
	Ixx = np.sum(masses * (y**2 + z**2))
	Iyy = np.sum(masses * (x**2 + z**2))
	Izz = np.sum(masses * (x**2 + y**2))
	Ixy = -np.sum(masses * x * y)
	Iyz = -np.sum(masses * y * z)
	Ixz = -np.sum(masses * x * z)
	I = np.array([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])
	return I

#Function to get the principal moments of inertia from the inertia matrix through diagonalization
def get_principal_moi(I):
	Ip = np.linalg.eigvals(I)
	# Sort and convert principal moments of inertia to SI (kg.m2)
	Ip.sort()
	return Ip

#Function to calculate euclidean distance between two points in Cartesian space of 3 dimensions
def euclidean_dist(x1, y1, z1, x2, y2, z2):
	return(math.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2))

#Function to calculate roughness of a given point and probe radius
def roughness(x, r):
	val = math.log(x)/math.log(r)
	return val

#------------------------------------------------------------------------------------------------------------------------------------------


#-----------------------------------------------------------MAIN PROGRAM-------------------------------------------------------------------
features_final = []

for j, row in input_df.iterrows():
	feat_dict = {}

	pdb_id = row["PDB_ID"]
	resi_list = row["Binding_site_residues"].split(",")

	feat_dict["PDB_ID"] = pdb_id

	ions = ['NCO', 'NA', 'RHD', 'SR', 'PO4', 'ZN', 'MG', 'AG', 'K', 'AU3', 'NI', 'SE4', 'ACT', 'IRI', 'IR', 'CL', 'LU', 'TB', 'CA', 'CS', 'TL', 'SO4', 'CU', 'CAC', 'AU', 'BR', 'SM', 'CD', 'MN', 'NH4', 'IR3', 'CO', 'HG', '3CO', 'OS', 'CON', 'BA', 'FE2', 'PB', 'SIN', 'MES', 'GOL', 'MPD']
	std_bases = ["A", "C", "G", "U"]  #XXX: Helps omit non-standard bases and ligands, if present
	resi_list_refined =  [x for x in resi_list if x.split("_")[0] not in ions and x.split("_")[0]!="HOH" and x.split("_")[0] in std_bases]

	nmr_flag = 0
	target_atoms = []
	with open(pdb_path+pdb_id+".pdb") as f:
		for line in f.readlines():
			line = line.strip()
			if(re.search("^MODEL", line)):
				nmr_flag = 1
				model_atomdict = extract_modelwise_atoms(pdb_id)
				break

	if(nmr_flag==1):
		target_atoms = model_atomdict[0]
	else:
		with open(pdb_path+pdb_id+".pdb") as f:
			for line in f.readlines():
				line = line.strip()
				if(line.startswith("ATOM")):
					target_atoms.append(line)

	#Binding pocket extraction
	binding_pocket = extract_binding_site(target_atoms, resi_list_refined)

	#--------------------------------------------------------------------------------------------------------------------------------
	#Composition features: A, C, G, U, Purines, Pyrimidines, Total atoms, Total residues, Total heavy atoms
	bases = [x.split("_")[0] for x in resi_list_refined]
	base_count = Counter(bases)
	feat_dict["A"] = base_count["A"]
	feat_dict["C"] = base_count["C"]
	feat_dict["G"] = base_count["G"]
	feat_dict["U"] = base_count["U"]
	feat_dict["purines"] = base_count["A"] + base_count["G"]
	feat_dict["pyrimidines"] = base_count["C"] + base_count["U"]
	feat_dict["total_atoms"] = len(binding_pocket)
	feat_dict["total_residues"] = len(bases)
	#--------------------------------------------------------------------------------------------------------------------------------

	#--------------------------------------------------------------------------------------------------------------------------------
	#Pharmacophore features: HBA (ACC), HBD (DON), Polar (POL), Anion (ANI), Cation (CAT), Hydrophobic (HYD), Others (OTH)
	#XXX: Based on PATTY atom typing in RDKit
	patty_types = {"A": {"P":"OTH", "C5'":"HYD", "O5'":"ACC", "C4'":"HYD", "O4'":"ACC", "C3'":"HYD", "O3'":"POL", "C2'":"HYD", "O2'":"POL", "C1'":"HYD", "N1":"ACC", "C2":"HYD", "N3":"ACC", "C4":"OTH", "C5":"OTH", "C6":"OTH", "N6":"DON", "N7":"POL", "C8":"HYD", "N9":"POL", "OP1":"ANI", "OP2":"ANI"}, "C": {"P":"OTH", "C5'":"HYD", "O5'":"ACC", "C4'":"HYD", "O4'":"ACC", "C3'":"HYD", "O3'":"POL", "C2'":"HYD", "O2'":"POL", "C1'":"HYD", "N1":"OTH", "C2":"OTH", "O2":"ACC", "N3":"ACC", "C4":"OTH", "N4":"DON", "C5":"HYD", "C6":"HYD", "OP1":"ANI", "OP2":"ANI"}, "G": {"P":"OTH", "C5'":"HYD", "O5'":"ACC", "C4'":"HYD", "O4'":"ACC", "C3'":"HYD", "O3'":"POL", "C2'":"HYD", "O2'":"POL", "C1'":"HYD", "N1":"DON", "C2":"OTH", "N2":"DON", "N3":"ACC", "C4":"OTH", "C5":"OTH", "C6":"OTH", "O6":"ACC", "N7":"POL", "C8":"HYD", "N9":"POL", "OP1":"ANI", "OP2":"ANI"}, "U": {"P":"OTH", "C5'":"HYD", "O5'":"ACC", "C4'":"HYD", "O4'":"ACC", "C3'":"HYD", "O3'":"POL", "C2'":"HYD", "O2'":"POL", "C1'":"HYD", "N1":"OTH", "C2":"OTH", "O2":"ACC", "N3":"DON", "C4":"OTH", "O4":"ACC", "C5":"HYD", "C6":"HYD", "OP1":"ANI", "OP2":"ANI"}}

	patty_list = extract_patty_types(binding_pocket, patty_types)
	patty_counter = Counter(patty_list)
	feat_dict["HBA"] = patty_counter["ACC"]
	feat_dict["HBD"] = patty_counter["DON"]
	feat_dict["polar"] = patty_counter["POL"]
	feat_dict["anion"] = patty_counter["ANI"]  #XXX: RNA has only anion count, not cation count since phosphate is negatively charged
	feat_dict["cation"] = patty_counter["CAT"]
	feat_dict["hydrophobic"] = patty_counter["HYD"]
	feat_dict["others"] = patty_counter["OTH"]
	feat_dict["total_heavy_atoms"] = len(patty_list)
	#--------------------------------------------------------------------------------------------------------------------------------

	#--------------------------------------------------------------------------------------------------------------------------------
	#Surface area and depth features: Pocket depth, Molecular SA, Polar SA, Hydrophobic SA, SASA, Polar SASA, Hydrophobic SASA, Relative polar SA, Relative hydrophobic SA, Relative polar SASA, Relative hydrophobic SASA

	with warnings.catch_warnings():
		warnings.simplefilter('ignore')
		pdb = PDBParser().get_structure(pdb_id, pdb_path+pdb_id+".pdb")
		io = PDBIO()
		io.set_structure(pdb)

		relative_sasa = 0
		relative_polar_sasa = 0
		relative_np_sasa = 0
		resi_depths = []

		if(nmr_flag==1):
			sr = SASA.ShrakeRupley()
			sr.compute(pdb[model_id], level="A")

			surface = get_surface(pdb[model_id], MSMS=msms_path)

			for model in pdb:
				if(model.id==model_id):
					for chain in model:
						for resi in chain:
							resid = str(resi.resname)+"_"+str(resi.id[1])+"_"+str(chain.id)
							if(resid in resi_list_refined):
								for atom in resi:
									relative_sasa = relative_sasa + np.round(atom.sasa, 5)
									if(atom.id.startswith("C")):
										relative_np_sasa = relative_np_sasa + np.round(atom.sasa, 5)
									else:
										relative_polar_sasa = relative_polar_sasa + np.round(atom.sasa, 5)

								resi_depths.append(residue_depth(chain[resi.id[1]], surface))
		else:
			sr = SASA.ShrakeRupley()
			sr.compute(pdb[0], level="A")

			surface = get_surface(pdb[0], MSMS=msms_path)

			for model in pdb:
				for chain in model:
					for resi in chain:
						resid = str(resi.resname)+"_"+str(resi.id[1])+"_"+str(chain.id)
						if(resid in resi_list_refined):
							for atom in resi:
								relative_sasa = relative_sasa + np.round(atom.sasa, 5)
								if(atom.id.startswith("C")):
									relative_np_sasa = relative_np_sasa + np.round(atom.sasa, 5)
								else:
									relative_polar_sasa = relative_polar_sasa + np.round(atom.sasa, 5)

							resi_depths.append(residue_depth(chain[resi.id[1]], surface))


		feat_dict["r_sasa"] = relative_sasa
		feat_dict["r_polar_sasa"] = relative_polar_sasa
		feat_dict["r_nonpolar_sasa"] = relative_np_sasa
		feat_dict["pocket_depth"] = np.round(np.mean(resi_depths), 5)

	#XXX: Using NACCESS for rest of the calculations
	os.system("cp "+pdb_path+pdb_id+".pdb "+naccess_path)
	out_var1 = subprocess.check_output(naccess_path+" "+pdb_path+pdb_id+".pdb -c", shell=True, stderr=subprocess.STDOUT)

	mol_sa = 0.0
	polar_sa = 0.0
	non_polar_sa = 0.0
	with open(pdb_id+".asa") as f:
		for line in f.readlines():
			line = line.strip()
			contents = line.split(" ")
			final_contents = [x for x in contents if x!=""]
			atomic_area = float(final_contents[-2])
			mol_sa = mol_sa + atomic_area

			if(final_contents[2][0]=="C"):
				non_polar_sa = non_polar_sa + atomic_area
			else:	
				polar_sa = polar_sa + atomic_area

	os.system("rm "+pdb_id+".*")
	feat_dict["mol_sa"] = np.round(mol_sa, 5)
	feat_dict["polar_sa"] = np.round(polar_sa, 5)
	feat_dict["nonpolar_sa"] = np.round(non_polar_sa, 5)

	#--------------------------------------------------------------------------------------------------------------------------------

	#--------------------------------------------------------------------------------------------------------------------------------
	#Shape descriptors: 3 Principal moments, NPR1, NPR2, Asphericity, Eccentricity, Spherocity index, Inertial shape factor
	#XXX: Formulae for Asphericity, ..., Inertial shape factor taken from DrugPred_RNA SI
	#XXX: Code for PMI calculation adapted from https://scipython.com/book/chapter-6-numpy/problems/p65/the-moment-of-inertia-tensor/
	mass_dict = {"H" :[1  , 1.00794], "He":[2  , 4.002602], "Li":[3  , 6.941], "Be":[4  , 9.012182], "B" :[5  , 10.811], "C" :[6  , 12.0107], "N" :[7  , 14.0067], "O" :[8  , 15.9994], "F" :[9  , 18.9984032], "Ne":[10 , 20.1797], "Na":[11 , 22.98976928], "Mg":[12 , 24.3050], "Al":[13 , 26.9815386], "Si":[14 , 28.0855], "P" :[15 , 30.973762], "S" :[16 , 32.065], "Cl":[17 ,  35.453], "Ar":[18 , 39.948], "K" :[19 , 39.0983], "Ca":[20 , 40.078], "Sc":[21 , 44.955912], "Ti":[22 , 47.867], "V" :[23 , 50.9415], "Cr":[24 , 51.9961], "Mn":[25 , 54.93845], "Fe":[26 ,55.933195], "Co":[27 , 58.933195], "Ni":[28 , 58.6934], "Cu":[29 , 63.546], "Zn":[30 , 65.39], "Ga":[31 , 69.723], "Ge":[32 , 72.64], "As":[33 , 74.92160], "Se":[34 , 78.96], "Br":[35 , 79.904], "Kr":[36 , 83.798], "Rb":[37 , 85.4678], "Sr":[38 , 87.62], "Y" :[39 , 88.90585], "Zr":[40 , 91.224], "Nb":[41 , 92.90638], "Mo":[42 , 95.94], "Tc":[43 , 97.9072], "Ru":[44 , 101.07], "Rh":[45 , 102.90550], "Pd":[46 , 106.42], "Ag":[47 , 107.8682], "Cd":[48 , 112.411], "In":[49 , 114.818], "Sn":[50 , 118.710], "Sb":[51 , 121.760], "Te":[52 , 127.60], "I" :[53 , 126.90447], "Xe":[54 , 131.293], "Cs":[55 , 132.9054519], "Ba":[56 , 137.327], "La":[57 , 138.90547], "Ce":[58 , 140.116], "Pr":[59 , 140.90765], "Nd":[60 , 144.242], "Pm":[61 , 145], "Sm":[62 , 150.36], "Eu":[63 , 151.964], "Gd":[64 , 157.25], "Tb":[65 , 158.92535], "Dy":[66 , 162.500], "Ho":[67 , 164.93032], "Er":[68 , 167.259], "Tm":[69 , 168.93421], "Yb":[70 , 173.04], "Lu":[71 , 174.967], "Hf":[72 , 178.49], "Ta":[73 , 180.94788], "W" :[74 , 183.84], "Re":[75 , 186.207], "Os":[76 , 190.23], "Ir":[77 , 192.217], "Pt":[78 , 195.084], "Au":[79 , 196.966569], "Hg":[80 , 200.59], "Tl":[81 , 204.3833], "Pb":[82 , 207.2], "Bi":[83 , 208.98040], "Po":[84 , 208.9824], "At":[85 , 209.9871], "Rn":[86 , 222.0176], "Fr":[87 , 223], "Ra":[88 , 226], "Ac":[89 , 227], "Th":[90 , 232.03806], "Pa":[91 , 231.03588], "U" :[92 , 238.02891], "Np":[93 , 238.8486], "Pu":[94 , 242.8798], "Am":[95 , 244.8594], "Cm":[96 , 246.911], "Bk":[97 , 248.9266], "Cf":[98 , 252.9578], "Es":[99 , 253.9656], "Fm":[100, 259.0046], "Md":[101, 260.0124], "No":[102, 261.0202], "Lr":[103, 264.0436]}

	masses = []
	coords = []
	for atom in binding_pocket:
		atom_type = atom[76:78].strip()
		masses.append(mass_dict[atom_type][1])
		coords.append([float(atom[30:38].strip()), float(atom[38:46].strip()), float(atom[46:54].strip())])

	masses = np.asarray(masses)
	coords = np.asarray(coords)
	inertia_tensor = get_inertia_matrix(masses, coords)
	pmi = get_principal_moi(inertia_tensor)

	npr1 = np.round(pmi[0]/pmi[2], 5)
	npr2 = np.round(pmi[1]/pmi[2], 5)

	asphericity = np.round(0.5 * ((pmi[2]-pmi[0])**2 + (pmi[2]-pmi[1])**2 + (pmi[1]-pmi[0])**2)/(pmi[0]**2 + pmi[1]**2 + pmi[2]**2), 5)
	eccentricity = np.round(math.sqrt(pmi[2]**2 - pmi[0]**2)/pmi[2]**2, 5)
	spherocity_index = np.round((3*pmi[0])/(pmi[0]+pmi[1]+pmi[2]), 5)
	inertial_shape_factor = np.round(pmi[1]/(pmi[0]*pmi[2]), 5)

	feat_dict["pmi_1"] = pmi[0]
	feat_dict["pmi_2"] = pmi[1]
	feat_dict["pmi_3"] = pmi[2]
	feat_dict["npr1"] = npr1
	feat_dict["npr2"] = npr2
	feat_dict["asphericity"] = asphericity
	feat_dict["eccentricity"] = eccentricity
	feat_dict["spherocity_index"] = spherocity_index
	feat_dict["inertial_shape_factor"] = inertial_shape_factor
	#--------------------------------------------------------------------------------------------------------------------------------

	#--------------------------------------------------------------------------------------------------------------------------------
	#Surface roughness (fractal dimension) of protein - XXX: Formulae taken from: https://doi.org/10.1002/minf.201400090
	radii = [0.4, 0.8, 1.6, 3.2]
	with warnings.catch_warnings():
		warnings.simplefilter('ignore')
		if(nmr_flag==1):
			ses_points = ResidueDepth_Mod.get_surface(pdb[model_id], MSMS=msms_path)
		else:
			ses_points = ResidueDepth_Mod.get_surface(pdb[0], MSMS=msms_path)

	final_surface_points = []
	for atom in binding_pocket:
		dists = []
		for coord in ses_points:
			dist = euclidean_dist(float(atom[30:38].strip()), float(atom[38:46].strip()), float(atom[46:54].strip()), coord[0], coord[1], coord[2])
			dists.append(dist)
		min_dist_point = tuple(ses_points[np.argmin(dists)])   #XXX: Tuple is used to facilitate duplicate removal in next step
		if(min_dist_point not in final_surface_points):
			final_surface_points.append(min_dist_point)

	for r in radii:
		x_sum = 0
		N = len(final_surface_points)  #No. of SES points used in roughness calculation
		for pt1 in final_surface_points:
			for pt2 in final_surface_points:
				dist = euclidean_dist(pt1[0], pt1[1], pt1[2], pt2[0], pt2[1], pt2[2])

				#Heavyside step function implementation
				if(dist>r):
					x_sum = x_sum + 0
				else:
					x_sum = x_sum + 1

		x_val = x_sum/(N*N)
		D = roughness(x_val, r)
		feat_dict["roughness_"+str(r)] = D
	#--------------------------------------------------------------------------------------------------------------------------------

	#--------------------------------------------------------------------------------------------------------------------------------
	#Structural parameters: Base parameters, Sugar torsions, Pseudo torsions - XXX: Based on the x3DNA program stand-alone results
	out_var1 = subprocess.check_output(dssr_path+" --more -i="+pdb_path+pdb_id+".pdb", shell=True, stderr=subprocess.STDOUT)
	main_chain_conf_params = []
	pseudo_torsion_params = []
	sugar_params = []
	suiteness = []
	starts = []
	lines_all = []
	with open("dssr-torsions.txt") as torsions:
		line_index = 0
		for line in torsions.readlines():
			line = line.strip()
			if(line=="nt               alpha    beta   gamma   delta  epsilon   zeta     e-z        chi            phase-angle   sugar-type    ssZp     Dp    splay"):
				starts.append(line_index)
			elif(line=="nt                eta   theta     eta\'  theta\'    eta\"  theta\""):
				starts.append(line_index)
			elif(line=="nt                 v0      v1      v2      v3      v4      tm      P   Puckering"):
				starts.append(line_index)
			elif(line=="nt             bin    cluster   suiteness"):
				starts.append(line_index)
			line_index = line_index + 1
			lines_all.append(line)

	for i, idx in enumerate(starts):
		if(i==0):
			main_chain_conf_params = lines_all[idx:starts[i+1]-1]
		elif(i==1):
			pseudo_torsion_params = lines_all[idx:starts[i+1]-1]	
		elif(i==2):
			sugar_params = lines_all[idx:starts[i+1]-1]
		else:
			suiteness = lines_all[idx:]

	os.system("rm dssr-*")  #Clean up output files from the dssp calculation

	alpha = []
	beta = []
	gamma = []
	delta = []
	epsilon = []
	zeta = []
	e_z = []
	chi = []
	phase_angle = []
	sugar_type = []
	ssZp = []
	Dp = []
	splay = []
	eta = []
	theta = []
	eta_p = []
	theta_p = []
	eta_dp = []
	theta_dp = []
	v0 = []
	v1 = []
	v2 = []
	v3 = []
	v4 = []
	tm = []
	P = []
	puckering = []
	suite = []

	bs_resi_ids = []
	for resi in resi_list_refined:
		resi_parts = resi.split('_')
		bs_resi_ids.append(resi_parts[2]+"."+resi_parts[0]+resi_parts[1])  #Convert to residue ID format as per the DSSR norms

	#[1, G, A.G1, ---, ---, 38.7, 84.9, -136.2, -70.2, -66(BI), -170.4(anti), 358.5(C2'-exo), ~C3'-endo, 4.72, 4.81, 22.29]
	for i, param in enumerate(main_chain_conf_params):
		try:
			if(i>0 and param!="******************************************************************************************"):
				contents = [x for x in param.split(" ") if x!=""]
				if(contents[2] in bs_resi_ids):
					if(contents[3]!="---"):
						alpha.append(float(contents[3]))
					if(contents[4]!="---"):
						beta.append(float(contents[4]))
					if(contents[5]!="---"):
						gamma.append(float(contents[5]))
					if(contents[6]!="---"):
						delta.append(float(contents[6]))
					if(contents[7]!="---"):
						epsilon.append(float(contents[7]))
					if(contents[8]!="---"):
						zeta.append(float(contents[8]))
					if(contents[9]!="---"):
						e_z.append(float(contents[9].split("(")[0]))
					if(contents[10]!="---"):
						chi.append(float(contents[10].split("(")[0]))
					if(contents[11]!="---"):
						phase_angle.append(float(contents[11].split("(")[0]))
					if(contents[12]!="---"):
						sugar_type.append(contents[12])
					if(contents[13]!="---"):
						ssZp.append(float(contents[13]))
					if(contents[14]!="---"):
						Dp.append(float(contents[14]))
					if(contents[15]!="---"):
						splay.append(float(contents[15]))
			elif(param=="******************************************************************************************"):
				break
		except:
			print(contents)

	for i, param in enumerate(pseudo_torsion_params):
		try:
			if(i>0 and param!="******************************************************************************************"):
				contents = [x for x in param.split(" ") if x!=""]
				if(contents[2] in bs_resi_ids):
					if(contents[3]!="---"):
						eta.append(float(contents[3]))
					if(contents[4]!="---"):
						theta.append(float(contents[4]))
					if(contents[5]!="---"):
						eta_p.append(float(contents[5]))
					if(contents[6]!="---"):
						theta_p.append(float(contents[6]))
					if(contents[7]!="---"):
						eta_dp.append(float(contents[7]))
					if(contents[8]!="---"):
						theta_dp.append(float(contents[8]))
			elif(param=="******************************************************************************************"):
				break
		except:
			print(contents)

	for i, param in enumerate(sugar_params):
		try:
			if(i>0 and param!="******************************************************************************************"):
				contents = [x for x in param.split(" ") if x!=""]
				if(contents[2] in bs_resi_ids):
					if(contents[3]!="---"):
						v0.append(float(contents[3]))
					if(contents[4]!="---"):
						v1.append(float(contents[4]))
					if(contents[5]!="---"):
						v2.append(float(contents[5]))
					if(contents[6]!="---"):
						v3.append(float(contents[6]))
					if(contents[7]!="---"):
						v4.append(float(contents[7]))
					if(contents[8]!="---"):
						tm.append(float(contents[8]))
					if(contents[9]!="---"):
						P.append(float(contents[9]))
					if(contents[10]!="---"):
						puckering.append(contents[10])
			elif(param=="******************************************************************************************"):
				break
		except:
			print(contents)

	for i, param in enumerate(suiteness):
		try:
			if(param!="" and suiteness[i+1]!="Concatenated suite string per chain. To avoid confusion of lower case"):
				if(i>0):
					contents = [x for x in param.split(" ") if x!=""]
					if(contents[2] in bs_resi_ids):
						if(contents[5]!="---"):
							suite.append(float(contents[5]))
			else:
				break
		except:
			print(contents)

	alpha = np.round(np.mean(alpha), 5)
	beta = np.round(np.mean(beta), 5)
	gamma = np.round(np.mean(gamma), 5)
	delta = np.round(np.mean(delta), 5)
	epsilon = np.round(np.mean(epsilon), 5)
	zeta = np.round(np.mean(zeta), 5)
	e_z = np.round(np.mean(e_z), 5)
	chi = np.round(np.mean(chi), 5)
	phase_angle = np.round(np.mean(phase_angle), 5)
	sugar_type = Counter(sugar_type)
	ssZp = np.round(np.mean(ssZp), 5)
	Dp = np.round(np.mean(Dp), 5)
	splay = np.round(np.mean(splay), 5)
	eta = np.round(np.mean(eta), 5)
	theta = np.round(np.mean(theta), 5)
	eta_p = np.round(np.mean(eta_p), 5)
	theta_p = np.round(np.mean(theta_p), 5)
	eta_dp = np.round(np.mean(eta_dp), 5)
	theta_dp = np.round(np.mean(theta_dp), 5)
	v0 = np.round(np.mean(v0), 5)
	v1 = np.round(np.mean(v1), 5)
	v2 = np.round(np.mean(v2), 5)
	v3 = np.round(np.mean(v3), 5)
	v4 = np.round(np.mean(v4), 5)
	tm = np.round(np.mean(tm), 5)
	P = np.round(np.mean(P), 5)
	puckering = Counter(puckering)
	suite = np.round(np.mean(suite), 5)

	feat_dict["alpha"] = alpha
	feat_dict["beta"] = beta
	feat_dict["gamma"] = gamma
	feat_dict["delta"] = delta
	feat_dict["epsilon"] = epsilon
	feat_dict["zeta"] = zeta
	feat_dict["e_z"] = e_z
	feat_dict["chi"] = chi
	feat_dict["phase_angle"] = phase_angle
	feat_dict["st_C3'_endo"] = sugar_type["~C3'-endo"]
	feat_dict["st_C2'_endo"] = sugar_type["~C2'-endo"]
	feat_dict["st_C3'_exo"] = sugar_type["~C3'-exo"]
	feat_dict["st_C2'_exo"] = sugar_type["~C2'-exo"]
	feat_dict["ssZp"] = ssZp
	feat_dict["Dp"] = Dp
	feat_dict["splay"] = splay
	feat_dict["eta"] = eta
	feat_dict["eta'"] = eta_p
	feat_dict["theta'"] = theta_p
	feat_dict["eta''"] = eta_dp
	feat_dict["theta''"] = theta_dp
	feat_dict["v0"] = v0
	feat_dict["v1"] = v1
	feat_dict["v2"] = v2
	feat_dict["v3"] = v3
	feat_dict["v4"] = v4
	feat_dict["tm"] = tm
	feat_dict["P"] = P
	feat_dict["pucker_C3'_endo"] = puckering["C3'-endo"]
	feat_dict["pucker_C2'_endo"] = puckering["C2'-endo"]
	feat_dict["pucker_C3'_exo"] = puckering["C3'-exo"]
	feat_dict["pucker_C2'_exo"] = puckering["C2'-exo"]
	feat_dict["suite"] = suite
	#--------------------------------------------------------------------------------------------------------------------------------
	feat_dict["class"] = row["class"]
	
	with open(outpath+pdb_id+"_"+str(j+1)+"_features.out", "w") as out:
		key_str = "\t".join(list(feat_dict.keys()))
		vals = [str(x) for x in list(feat_dict.values())]
		val_str = "\t".join(vals)

		print(key_str, file=out)
		print(val_str, file=out)

	print(pdb_id+"_"+str(j))

































