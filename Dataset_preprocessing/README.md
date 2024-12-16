# Dataset preprocessing stages

## Source of the dataset
* Positive dataset curated from experimental RNA-small molecule complexes (/data/Positive_pockets.csv)
* Negative dataset curated from three strategies (/data/Negative_pockets.csv) - backbone motif search (automate_primos_v1.sh), exhaustive pocket prediction (automate_fpocket_v1.py), all-to-all blind docking (automate_nldock_v1.py)

## Source of the features
* RNA binding site features from custom implementation

## Description
* RNA structure-based features are obtained for each binding site in the positive and negative datasets.
* BioPython along with NACCESS, MSMS and X3DNA-DSSR tools (with appropriate academic license) are used to calculate the features for each binding site as listed in the article.
* Finally, features of both datasets are combined to get the final dataset for model training (/data/Pocket_dataset_final.csv).
* If you want to create your own negative dataset from a set of PDB structures, use the three automation scripts provided first and then perform overlap analysis using a suitable cut-off as discussed in the article.
* To run PRIMOS motif search calculations, a wormbase has to be created for the positive binding sites and all the complete RNA structures.

## Prerequisites
* Anaconda or Miniconda with Python 3.8.
* A conda environment with the following libraries:
	* Python (>v3.8)
	* pandas
	* numpy
	* scipy
	* scikit-learn
	* pickle
	* matplotlib
	* seaborn
	* biopython
	* fpocket
	* Generic libraries: sys, re, csv, os, collections...
* PRIMOS package (https://github.com/pylelab/PRIMOS)
* NLDock package (http://huanglab.phys.hust.edu.cn/software/NLDock/)
* MSMS package (https://ccsb.scripps.edu/msms/)
* NACCESS package (http://www.bioinf.manchester.ac.uk/naccess/)
* X3DNA-DSSR binary executable (https://x3dna.org/)

## Sample commands (Follow the same order to reproduce the files provided in sample_output folder)
```
* python compute_pocket_features_v1.py "./data/Positive_pockets.csv" "tsv" "./sample_input/" "MSMS path" "NACCESS path" "DSSR path" "./sample_output/Pos_feats/"
* python compute_pocket_features_v1.py "./data/Negative_pockets.csv" "tsv" "./sample_input/" "MSMS path" "NACCESS path" "DSSR path" "./sample_output/Neg_feats/"
* python create_final_dataset_v1.py "./data/Positive_pockets.csv" "./data/Family_mapping_v1.xlsx" "./sample_output/Positive_dataset_final.csv"
* python create_final_dataset_v1.py "./data/Negative_pockets.csv" "./data/Family_mapping_v1.xlsx" "./sample_output/Negative_dataset_final.csv"
* python create_final_dataset_v2.py "./sample_output/Positive_dataset_final.csv" "./sample_output/Negative_dataset_final.csv" "./sample_output/Pocket_dataset_final.csv"
```

## Miscellaneous
* The complete dataset of RNA-small molecule complex structures can be downloaded from the <a href="https://www.rcsb.org/" target="_blank">PDB database</a>
* All programs take inputs as command-line arguments only
* The same pipeline is followed to prepare all custom test datasets used in the study, in the format necessary for model evaluation
* The header fields or columns provided in the sample dataset must be present (column order can be random) in your custom dataset file for the programs to work in the intended fashion
