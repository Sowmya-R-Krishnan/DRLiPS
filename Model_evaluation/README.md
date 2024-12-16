# Model evaluation with external test datasets - Regression and classification data

## Source of the datasets
* Studies on genome-wide druggable regions in SARS-CoV-2 - Manfredonia et al., 2020; Sreeramulu et al., 2021
* pre-miR-21 binding site - Costales et al., 2018

## Description
* The dataset provided has been pre-processed as per the pipeline explained in Dataset_preprocessing README file and provided here. They can be directly used for model evaluation.

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
	* Generic libraries: sys, re, csv, os, collections...

## Sample commands (Follow the same order to reproduce the files provided in sample_output folder)
```
* python test_on_external_set.py "./data/External_set_final_feats.csv" "../Feature_selection/data/Pocket_dataset_final.csv" "best_models.log" "test_results.csv" 4
```

## Miscellaneous
* All programs take inputs as command-line arguments only
* The same pipeline is followed to prepare all custom test datasets used in the study, in the format necessary for model evaluation
* The header fields or columns provided in the sample dataset must be present (column order can be random) in your custom dataset file for the programs to work in the intended fashion
