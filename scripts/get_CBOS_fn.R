get_CBOS <- function(tDDqMAT,targ1){
    #
    NCID_DTdbout <- data.frame(matrix(0, ncol = 1, nrow = dim(tDDqMAT)[1]))
    row.names(NCID_DTdbout) <- row.names(tDDqMAT)
    colnames(NCID_DTdbout) <- TARG
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
            COMBO_TRUSCORE <- max(DDqMAT[row.names(DDqMAT)%in% aa, j,drop=FALSE])
            # above combotru for prot i and query j
            # get the randmax and randstd for prot i and query j
            RANDOM_MAX <- vector()
            listnum <- which(names(RUN_RAND) == colnames(D_T1))
            RANDMAXSD <- RUN_RAND[[listnum]][row.names(RUN_RAND[[listnum]]) %in% row.names(tDDqMAT)[j], ]
            RANDMAX_MEAN <- RANDMAXSD[1,1] 
            RANDMAX_SD <- RANDMAXSD[1,2]
            TARG_SCORE <- (COMBO_TRUSCORE - RANDMAX_MEAN)/RANDMAX_SD
            qq<- which(row.names(NCID_DTdbout)==colnames(DDqMAT)[j])
            #zz <- which(colnames(DTdbout)== colnames(D_T1))
            NCID_DTdbout[qq,1] <- TARG_SCORE
        }
    }
RUN_APPD246cbos <- NCID_DTdbout
############# 
return(list(RUN_APPD246cbos))
}