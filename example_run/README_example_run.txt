This directory contains worked out example/data for one target – arylhydrocarbon receptor (AHR). The same scripts can be used for other targets by just changing the target name. 

Description:
1)
1_AHR_Input_sdf_other_files - This directory has the input sdf file. The scripts used to process the sdf file through openeye programs (Filter,omega) are provided in the directory. 
(Note: The output final data after all pre-processing/filtering for all targets are provided in: https://github.com/BHSAI/CBOS/tree/master/Input_sdf_files_for_ROCS)

2)
2_AHR_ROCS_output_process - This directory has the R script used to process the .rpt file into a chemical vs chemical similarity matrix is provided in the directory.

3) 
3_AHR_MAX_CBOS_scores - This directory has the R scripts to run CBOS and MAX score (starting from chemical vs chemical similarity matrix) and get AUC/EF10 values as well as associated data required for the run. Note: If interested to run CBOS only and get AUC values, please run this script: AHR_Script_for_CBOSscore_AUC.R.

