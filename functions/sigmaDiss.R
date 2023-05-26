sigmaDiss <- function(target, refSet, ICVdata, truncThresh = 0){

  # CHECK col dims of data match
  # CHECK refSet > 1
  # CHECK ICVdata > 1

  ICVsd <- apply(ICVdata, 2, sd, na.rm=T)
  refSetS <- sweep(refSet, 2, ICVsd, "/")
  targetS <- sweep(target, 2, ICVsd, "/")
  ICVdataS <- sweep(ICVdata, 2, ICVsd, "/")

  if(sum(complete.cases(ICVdata)) < nrow(ICVdata)){
    warning("ICVdata contains incomplete observations: ",
            sum(!complete.cases(ICVdata)),
            " obs dropped prior to PCA")
    ICVdataSC <- ICVdataS[complete.cases(ICVdataS),]
  }

  ICVpca <- prcomp(ICVdataSC)
  PCs <- ICVpca$sdev / sum(ICVpca$sdev)
  truncPCs <- which(PCs > truncThresh)

  if(length(truncPCs)==0){
    stop(paste0("truncThresh excludes all PC axes: PC1 explains ", round(PCs[1]*100,2), "% of variation"))
  }

  if(sum(PCs[truncPCs] < 1e-10) > 0){
  warning("Truncating additional PC axis with zero variation")
  truncPCs <- which(PCs > truncThresh & PCs[truncPCs] >= 1e-10)
  }

  targetPred <- as.data.frame(predict(ICVpca, targetS))
  refSetPred <- as.data.frame(predict(ICVpca, refSetS))
  ICVPred <- as.data.frame(predict(ICVpca, ICVdataS))

  ICVPredsd <- apply(ICVPred, 2 , sd, na.rm=TRUE)
  targetPredS <- sweep(targetPred, 2, ICVPredsd, "/")
  refSetPredS <- sweep(refSetPred, 2, ICVPredsd, "/")

  neighDists <- as.vector(get.knnx(data=refSetPredS[,truncPCs],
                                   query=targetPredS[,truncPCs],
                                   k=1,
                                   algorithm="brute")[[2]])

  neighChi <- pchi(neighDists, max(truncPCs))

  neighSigma <- qchi(neighChi, 1)

  return(neighSigma)
}
