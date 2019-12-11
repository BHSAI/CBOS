#----------------------------
# Script to get AUC for CBOS score from chem vs chem similarity matrix table
# ---------------------------
# Supporting information: Abdulhameed et al. Correction based on shuffling (CBOS) 
# improves 3D ligand-based virtual screening 
# For any questions, please email: mabdulhameed@bhsai.org
#-----------------------------
# setwd("./PGR")

nTARGET <- 13
TARG <- 'PGR'
#************************** 
# Load Train/test split for TARG
#************************** 
TAR_ACTANNOT <- paste(TARG,"_agANT_combined_ACTIVITY_forROCplot.txt",sep="")
#######################
FILE_3DCBOS_RDATA <- paste(TARG,"_TRAIN_TEST_split.RData",sep="")
load(FILE_3DCBOS_RDATA)
#######################

TR_COUNT <- nrow(TARG_TR1)
#######################
# uploading the random data used in our run.
#
load(paste("./",nTARGET,"_",TARG,"_all_RANDLIST1000.RData",sep=""))

#######################
RUN_APPD <- TARG3d_DD
# calculate MAX scores between each query vs 1000 random TR sets
RAND_AVGSD1000 <-  vector("list",1) 
#*********************
# used spply twice to get the Max score between each query and 1000 random reference set - output is a matrix (1000 * nCOMPOUNDS)
Q_RANDMAX <- sapply(colnames(RUN_APPD),function(y) sapply(RAND_bb3,function (x) max(RUN_APPD[row.names(RUN_APPD)%in% x,y,drop=FALSE])))  
RAND_MEAN <- round(colMeans(Q_RANDMAX),2)
RAND_SD   <-round(apply(Q_RANDMAX,2,sd),2)
Q_AVGSD <- rbind(Q_RANDMAX,RAND_MEAN,RAND_SD)

TAR1 <- data.frame(t(Q_AVGSD[1001:1002,]))
TAR1$TAR_NAME <- TARG
TAR1$DRU_NAME <- row.names(TAR1)

RAND_AVGSD1000[[1]] <- TAR1
names(RAND_AVGSD1000)[1]<- as.character(unique(RAND_AVGSD1000[[1]]$TAR_NAME))
RUN_RAND <- RAND_AVGSD1000
#*****
# 
DDqMAT <- RUN_APPD
tDDqMAT <- t(DDqMAT)
# 
source("./get_CBOS_fn.R")
CBOS_out1 <- get_CBOS(tDDqMAT,TARG)
RUN_APPD246 <- CBOS_out1[[1]]

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
# Subset data to only query/test set for performance (ROC AUC) evaluation
RUN_APPD246_2 <- RUN_APPD246[row.names(RUN_APPD246) %in% TARG_Q$CID,,drop=FALSE]
TAR_APPD246 <- RUN_APPD246_2
TAR_val1_try1 <- PREP_4ROCplot_1(TAR_ACTANNOT,TAR_APPD246)
#install.packages('ROCR')
library(ROCR)
pred <- prediction(TAR_val1_try1[,2],TAR_val1_try1$ACT_AN)
pred_test <- performance(pred, "tpr", "fpr")
auc_ROCR <- performance(pred, measure = "auc")
CBOSROCAUC <-round(auc_ROCR@y.values[[1]],2)
######################
# EF10%
library('enrichvs')  
CBOS_e10 <- round(enrichment_factor(TAR_val1_try1[,2],TAR_val1_try1$ACT_AN,top=0.1,decreasing=TRUE),2)

nACTIV <- nrow(TAR_val1_try1[TAR_val1_try1$ACT_AN==1,]) # number of actives in query/test set
nINACTIV <- nrow(TAR_val1_try1[TAR_val1_try1$ACT_AN==0,]) # number of inactives in query/test set
data.frame(TARG,CBOSROCAUC,CBOS_e10,nACTIV,nINACTIV)
########################################### 
rm(list=ls())




