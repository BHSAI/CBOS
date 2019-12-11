The current directory (./3_AHR_MAX_CBOS_scores) contains the R scripts to perform CBOS analysis 
and also to calculate ROC AUC/Enrichment factor 10% (EF10) values for CBOS and MAX score for arylhydrocarbon receptor (AHR) data.

Script for CBOS score and AUC/EF10:
	- AHR_Script_for_CBOSscore_AUC.R

Script for MAX score and AUC/EF10:
	- AHR_Script_for_MAXscore_AUC.R

Output of each of these two scripts will be a table with AUC and EF10 values. 

Note: 
1. In order to reproduce our results we have provided the reference set and query/test set 
used in our screening performance evaluation as .RData file (AHR_TRAIN_TEST_split.RData).

2. We have also provided the codes used to start from chemical vs chemical similarity matrix 
(obtained by processing ROCS output .rpt file) and get a new reference and query/test set at the end of the R script (AHR_Script_for_CBOSscore_AUC.R).

