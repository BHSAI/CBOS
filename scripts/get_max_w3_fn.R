get_max_w3 <- function(tDDqMAT,targ1){
  #
  NCID_DTdbout <- data.frame(matrix(0, ncol = 1, nrow = dim(tDDqMAT)[1]))
  row.names(NCID_DTdbout) <- row.names(tDDqMAT)
  colnames(NCID_DTdbout) <- targ1
  NCID_DTdbout_w3 <-NCID_DTdbout
  #
  D_T1_pre <- TARG_TR1
  row.names(D_T1_pre) <- D_T1_pre$CID
  D_T1_pre$ERA <- 1
  colnames(D_T1_pre)[3] <- TARG
  D_T1 <- D_T1_pre[,3,drop=FALSE]
  D_T1 <- D_T1[row.names(D_T1) %in% row.names(DDqMAT),,drop=FALSE]
  
  if (nrow(D_T1) != 0){
    aa <- row.names(D_T1)
    # loop through drugs 
    for (j in 1:dim(tDDqMAT)[1]){
      # subset DDMAT to retain all rows with DBids in 'aa' and find max
      test <- max(DDqMAT[row.names(DDqMAT)%in% aa, j,drop=FALSE])
      qq<- which(row.names(NCID_DTdbout)==colnames(DDqMAT)[j])
      # zz <- which(colnames(DTdbout)== colnames(D_T1))
      NCID_DTdbout[qq,1] <- test
      
      # Instead of MAX score, if we need to explore weighted average of three nearest neighbors
      w3_1 <- DDqMAT[row.names(DDqMAT)%in% aa, j,drop=FALSE]
      w3_1 <- w3_1[order(-w3_1),,drop=FALSE]
      w3_1 <- 2-(w3_1)
      gg <- 0
      for (hh in 1:3) {
        gg<- gg+ round(exp(-((w3_1[hh,1]/0.4)^2)),3)
      }
      w3<- round(gg/3,3)
      qq_w3<- which(row.names(NCID_DTdbout)==colnames(DDqMAT)[j])
      # zz_w3 <- which(colnames(DTdbout)== colnames(D_T1))
      NCID_DTdbout_w3[qq_w3,1] <- w3
    }
  }
  RUN_APPD246 <- NCID_DTdbout
  RUN_APPD246_w3 <- NCID_DTdbout_w3
  ############# 
  return(list(RUN_APPD246,RUN_APPD246_w3))
}