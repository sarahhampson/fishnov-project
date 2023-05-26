noveldist <- function(x, method, dissWeight = NULL, abundWeight = FALSE){

  require(vegan)

  if(method=="weighted" & is.null(dissWeight)){
    stop('A weight distance matrix must be supplied if method = "weighted". See ?noveldist for details')
  }

  if(!is.null(dissWeight) | method == "weighted"){
    require(picante)
    if(abundWeight){warning("Weighting conducted on individuals. see ?noveldist for details")}
    if(!abundWeight){warning("Weighting conducted on taxa. see ?noveldist for details")}

    return(comdist(x, dis=dissWeight, abundance.weighted = abundWeight))
  }

  # test vegdist compatibility
  vegdistMethod = method %in% c("manhattan", "euclidean", "canberra", "clark",
                                "bray", "kulczynski", "jaccard", "gower",
                                "altGower", "morisita", "horn", "mountford",
                                "raup", "binomial", "chao", "cao",
                                "mahalanobis", "chisq", "chord",
                                "aitchison", "robust.aitchison")

  if(vegdistMethod){
    distMat <- as.matrix(vegdist(x, method=method))
  } else {
    distMat <- as.matrix(designdist(x, method=method))
  }

}
