#----------------------------
# Script to process output .rpt file from ROCS program and
# generate chem vs chem matrix with ROCS combo scores
# ---------------------------
# Supporting information: Abdulhameed et al. Correction based on shuffling (CBOS) 
# improves 3D ligand-based virtual screening 
#-----------------------------

# install.packages("splitstackshape") # for cSplit function
library(splitstackshape)
# setwd("~/AHR_ROCS_output_process")
# load ROCS output .rpt file
RUNrpt <- read.delim("AHRrocs_1.rpt",header=TRUE,sep="\t", stringsAsFactors = FALSE)
RUNrpt_split <- split(RUNrpt,RUNrpt$ShapeQuery)
# get starting chem-IDs column for merge
testdf1 <- RUNrpt_split[[1]]
testdf1 <- cSplit(testdf1,"Name","_")
testdf2 <- testdf1[,c(17,1,13)]

is.factor(testdf2$ShapeQuery)
# if (is.factor) is true then:
#install.packages("tidyverse")
#library(tidyverse)
# testdf2 <- testdf2 %>% mutate_if(is.factor, as.character)

DD_START <- testdf2[,1]
#--------------------------
#### process ROCS output and get chem vs chem table ###
# split - list - get tables - process - merge
# start chunk-1 (Note: to speed up computation, the ROCS was run using two chunks of input data)
#--------------------------
dim(DD_START)
for (i in 1:length(RUNrpt_split)){
  DFin1 <- RUNrpt_split[[i]]
  DFin1 <- cSplit(DFin1,"Name","_")
  DFsub1 <- DFin1[,c(17,1,13)]
  colnames(DFsub1)[3] <- unique(DFsub1$ShapeQuery)
  DFsub1 <- DFsub1[,c(1,3)]
  DD_START<- merge(DD_START,DFsub1,by="Name_1")
}
write.table(DD_START,file="AHRallvsall_MATRIX_ROCS_1.txt",row.names=FALSE,sep="\t")
##################################### 
rm(list=ls())
##################################### 
#------------
# chunk-2
#------------
# load ROCS output .rpt file
RUNrpt <- read.delim("AHRrocs_2.rpt",header=TRUE,sep="\t", stringsAsFactors = FALSE)
RUNrpt_split <- split(RUNrpt,RUNrpt$ShapeQuery)
# get starting chem-IDs column for merge
testdf1 <- RUNrpt_split[[1]]
testdf1 <- cSplit(testdf1,"Name","_")
testdf2 <- testdf1[,c(17,1,13)]
# 
DD_START <- testdf2[,1]
#--------------------------
#### process ROCS output and get chem vs chem table
# start chunk-2 
#--------------------------
dim(DD_START)
for (i in 1:length(RUNrpt_split)){
  DFin1 <- RUNrpt_split[[i]]
  DFin1 <- cSplit(DFin1,"Name","_")
  DFsub1 <- DFin1[,c(17,1,13)]
  colnames(DFsub1)[3] <- unique(DFsub1$ShapeQuery)
  DFsub1 <- DFsub1[,c(1,3)]
  DD_START<- merge(DD_START,DFsub1,by="Name_1")
}
write.table(DD_START,file="AHRallvsall_MATRIX_ROCS_2.txt",row.names=FALSE,sep="\t")
############################## 
rm(list=ls())
############################## 
#------
# merge and obtain final all chem vs chem matrix with ROCS combo scores

AHR_ROCS1 <- read.delim("AHRallvsall_MATRIX_ROCS_1.txt",header=TRUE,sep="\t", stringsAsFactors = FALSE)
AHR_ROCS2 <- read.delim("AHRallvsall_MATRIX_ROCS_2.txt",header=TRUE,sep="\t", stringsAsFactors = FALSE)

AHRall_MATRIX_ROCS_merge <- merge(AHR_ROCS1,AHR_ROCS2,by="Name_1")
write.table(AHRall_MATRIX_ROCS_merge,"AHRall_MATRIX_ROCS_merge.txt",row.names=FALSE,sep="\t")
###################################### 
rm(list=ls())

