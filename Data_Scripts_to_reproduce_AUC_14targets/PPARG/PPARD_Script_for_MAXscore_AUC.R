#----------------------------
# Script to get AUC for MAX score from ROCS output chem vs chem table
# ---------------------------
# Supporting information: Abdulhameed et al. Correction based on shuffling (CBOS) 
# improves 3D ligand-based virtual screening 
# For any questions, please email: mabdulhameed@bhsai.org
#-----------------------------
# setwd("./Data_Scripts_to_reproduce_AUC_14targets/PPARG")
nTARGET <- 7
TARG <- 'PPARG'
#************************** 
# Load Train/test split for TARG
#************************** 
TAR_ACTANNOT <- paste(TARG,"_agANT_combined_ACTIVITY_forROCplot.txt",sep="")
#######################
FILE_3DCBOS_RDATA <- paste(TARG,"_TRAIN_TEST_split.RData",sep="")
load(FILE_3DCBOS_RDATA)
#######################

#######################
DDqMAT <- TARG3d_DD[row.names(TARG3d_DD) %in% TARG_TR1$CID,colnames(TARG3d_DD) %in% TARG_Q$CID]
tDDqMAT <- t(DDqMAT)
# 
source("./get_max_w3_fn.R")
MAXw3_out1 <- get_max_w3(tDDqMAT,TARG) 
RUN_APPD246 <- MAXw3_out1[[1]] # MAX score
RUN_APPD246_w3 <- MAXw3_out1[[2]]

############# 
# get ROC AUC 
PREP_4ROCplot_1 <- function(TAR_ACTANNOT,TAR_APPD246,TAR_uniprot){
    TARannot <- read.delim(TAR_ACTANNOT, header=TRUE,sep="\t",stringsAsFactors = FALSE)
    rownames(TARannot) <- paste("X",TARannot$CID,sep="")
    TARannot <- TARannot[,-1,drop=FALSE]
    TAR_val1 <- merge(TARannot,TAR_APPD246,by="row.names")
    # TAR_val1_try1 <- TAR_val1[,c(1:2,which(colnames(TAR_val1)==TAR_uniprot))]
    TAR_val1_try1 <- TAR_val1
    row.names(TAR_val1_try1) <- TAR_val1$Row.names
    TAR_val1_try1 <- TAR_val1_try1[,2:3]
    TAR_val1_try1$ACT_AN <- sub("Inactive",0,TAR_val1_try1$ACT_AN)
    TAR_val1_try1$ACT_AN <- sub("Active",1,TAR_val1_try1$ACT_AN)
    TAR_val1_try1 <- TAR_val1_try1[order(-TAR_val1_try1[,2]),]
    TAR_val1_try1$ACT_AN <- as.numeric(TAR_val1_try1$ACT_AN)
    return(TAR_val1_try1)
}
############
TAR_APPD246 <- RUN_APPD246
TAR_val1_try1 <- PREP_4ROCplot_1(TAR_ACTANNOT,RUN_APPD246)
#install.packages('ROCR')
library(ROCR)
pred <- prediction(TAR_val1_try1[,2],TAR_val1_try1$ACT_AN)
pred_test <- performance(pred, "tpr", "fpr")
auc_ROCR <- performance(pred, measure = "auc")
MAXROCAUC <-round(auc_ROCR@y.values[[1]],2)
######################
# EF10%
library('enrichvs')  
MAX_e10 <- round(enrichment_factor(TAR_val1_try1[,2],TAR_val1_try1$ACT_AN,top=0.1,decreasing=TRUE),2)
######################
nACTIV <- nrow(TAR_val1_try1[TAR_val1_try1$ACT_AN==1,]) # number of actives in query/test set
nINACTIV <- nrow(TAR_val1_try1[TAR_val1_try1$ACT_AN==0,]) # number of inactives in query/test set
data.frame(TARG,MAXROCAUC,MAX_e10,nACTIV,nINACTIV)
###############
rm(list=ls())
