# x is vector or matrix of target entities
# refSet is matrix-like of ref entities (must have same dimensionality as x)
# method is vegdist method or designdist method
# nSize is number of ranked reference dissimilarities to return

novelty <- function(x, refSet, method, nSize, dissWeight=NULL, abundWeight=FALSE, ...){

  require(vegan)

  # coerce reference set is a matrix
  if(class(refSet)[1] != "matrix"){

    error.trap <- tryCatch(as.matrix(refSet), error=function(e) e)
    if(inherits(error.trap, "error")){
      stop("Reference set must be a matrix or matrix-like object")
    }

    refSet <- as.matrix(refSet)}

  if(nrow(refSet)==1){
    warning("Reference set is of size 1: novelty will be nonsensical")
    }

  if(class(x)[1] != "matrix"){
    x = matrix(x, ncol=length(x))
  }

  if(ncol(x) != ncol(refSet)){
    stop("Target and reference set must have the same number of columns")
  }

  if(!mode(x) %in% c("numeric", "double", "logical") |
     !mode(refSet) %in% c("numeric", "double", "logical")){
  stop("Dissimilarity can only be calculated using binary or numeric values")
  }

  if(nSize > nrow(refSet)){
    warning(paste0("nSize has been corrected to ", nrow(refSet), " (size of reference set)"))
    nSize <- nrow(refSet)
  }

  if(method=="weighted"){
  xMat <- as.matrix(noveldist(rbind(x, refSet),
                              method=method,
                              dissWeight=dissWeight,
                              abundWeight=abundWeight))[-(1:nrow(x)), 1:nrow(x)]
  } else {
  xMat <- as.matrix(noveldist(rbind(x, refSet), method=method))[-(1:nrow(x)), 1:nrow(x)]
  }

  if(class(xMat)[1]=="matrix"){

    return(t(apply(xMat, 1, sort)[1:nSize, 1:nrow(x)]))
  } else {
    return(sort(xMat)[1:nSize])

  }

  return(t(xMat))

  }
