# CBOS
Correction based on shuffling for 3D virtual screening

Supporting information for paper: Abdulhameed et al. Correction based on shuffling (CBOS) improves 3D ligand-based virtual screening (Under Review).

--------------------------- 

 
1. example_run
This directory contains data and scripts to perform all steps reported in the paper starting from input .sdf file. 

A typical run starts with sdf file of compounds with known activity for a target. The sdf file of compounds are passed through openeye programs (filter, omega, and rocs). 

The output from ROCS is .rpt file. R script is used to convert the .rpt file into a chemical vs chemical similarity matrix. 

Finally, R scripts are used to start from the chemical vs chemical similarity matrix and calcualte CBOS score. In this paper, for virtual screening evaluation, we split the actives into reference compounds and keep the rest of actives along with inactives as the screening set. We evaluated the standard MAX score and CBOS score performance using ROC AUC values and enrichment factor 10% (EF10) using the screening set.

"example_run" has 3 sub-directories:
 - 1_AHR_Input_sdf_other_files

 - 2_AHR_ROCS_output_process

 - 3_AHR_MAX_CBOS_scores
(please refer to readme in each sub-directory)


2. Data_scripts_reproduce_AUC_14targets
This directory has 14 sub-directory for each of the 14 targets analyzed. Each sub-directory has data and R scripts that will run either CBOS or MAX score and provide the AUC and EF10 for each target (as reported in the paper). Please run the R scripts in each sub-directory to obtain the AUC/EF values. The data points together will create Figure 1 in the paper. 


