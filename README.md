# A Novel Method for Prediction of Druggable RNA-Small Molecule Binding Pockets Using Machine Learning
**Sowmya Ramaswamy Krishnan, Arijit Roy<sup>#</sup>, Limsoon Wong, M. Michael Gromiha<sup>#</sup>**

**NOTE**: This code should be used for academic purposes only. This code should not be shared or used for commercial purposes without the consent of all the authors involved.

**Description**: Python implementation of the programs used for development of the models hosted in Druggable RNA-Ligand binding Pocket Selector (<a href="https://web.iitm.ac.in/bioinfo2/DRLiPS/" target="_blank">DRLiPS</a>) web server

<a href="https://web.iitm.ac.in/bioinfo2/DRLiPS/" target="_blank"><img align=center src="https://web.iitm.ac.in/bioinfo2/DRLiPS/images/logo/RSS_logo_v1.png"></a>

# Usage restrictions
The code has been provided for academic purposes only.

# Data
The sample datasets used for training and evaluating the models are provided in the `data` folder. In order to use your own data or the complete dataset, you have to prepare the datasets in the format provided for each case to train/test the models. To get the complete datasets used in the study, please visit the <a href="https://www.rcsb.org/" target="_blank">PDB database</a>.

# Dependencies not covered
The feature calculations for binding sites require the <a href="https://x3dna.org/" target="_blank">X3DNA-DSSR</a> binary executable. Since it is a third-party licensed tool, it will not be included here. Kindly ensure that the executable is available in your system before running the code.

# Code usage
Instructions on code usage are available in a README.md file within each folder (links also provided below).
* <a href="https://github.com/Sowmya-R-Krishnan/DRLiPS/blob/main/Dataset_preprocessing/README.md">Dataset pre-processing</a>
* <a href="https://github.com/Sowmya-R-Krishnan/DRLiPS/blob/main/Feature_selection/README.md">Feature selection</a>
* <a href="https://github.com/Sowmya-R-Krishnan/DRLiPS/blob/main/Model_evaluation/README.md">Model evaluation</a>

For further queries related to code usage, please write to us (sowmya.rk1@tcs.com; roy.arijit3@tcs.com; gromiha@iitm.ac.in).

# Citation
Please cite this article if you use the codes in this repository for your research: 
