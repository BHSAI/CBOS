Use the R script in this directory: AHR_script_process_ROCSoutput_to_chemvschem_table.R
to process the ROCS output (.rpt) files.

#-------------------------
Note-1 : To speed up calculation, we have split the query into two sets and performed two ROCS runs 
(check our ROCS run script in directory: ../1_AHR_Input_sdf_other_files/AHR_ROCS_runscript.txt)
#--------------------------

The current R script converts the .rpt file into one chemical vs chemical similarity matrix table.

We have provided the output file (AHRall_MATRIX_ROCS_merge_FINALout.zip) from our analysis in this directory. 
(Note: The actual .rpt file has size >1GB and hence could not be loaded into github and we have given the output table instead)
