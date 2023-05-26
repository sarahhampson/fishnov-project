# time series is matrix/matrix-like. Assumes first row is oldest and subsequent
# rows are ordered
# method can be in vegdist or a custom design from designmat
# reference can be "baseline", "pairwise" or an integer for time series row

turnover <- function(timeseries, method, reference, dissWeight=NULL, abundWeight=FALSE){

  require(vegan)

  # coerce reference set is a matrix
  if(class(timeseries)[1] != "matrix"){

    error.trap <- tryCatch(as.matrix(timeseries), error=function(e) e)
    if(inherits(error.trap, "error")){
      stop("Time series must be a matrix or matrix-like object")
    }

    timeseries <- as.matrix(timeseries)}

  if(nrow(timeseries)==1){
    stop("Time series must be longer than 1")
    }

  if(!mode(timeseries) %in% c("numeric", "double", "logical")){
  stop("Dissimilarity can only be calculated using binary or numeric values")
  }

  if(is.numeric(reference)){
    if(reference > nrow(timeseries) | reference < 1){
      stop("Selected reference is outside time series dimensions")
    }
  } else {
    if(!reference %in% c("baseline", "pairwise")){
      stop('reference must equal "baseline", "pairwise" or a integer indicating a row')
    }
  }

  if(method=="weighted"){
  xMat <- as.matrix(noveldist(x=timeseries, method=method,
                              dissWeight = dissWeight,
                              abundWeight=abundWeight))
  } else {
    xMat <- as.matrix(noveldist(x=timeseries, method=method))
  }


  if(reference == "baseline"){return(c(NA, xMat[1,-1]))}
  if(reference == "pairwise"){return(c(NA, diag(xMat[-1,-ncol(xMat)])))}

  return(ifelse(xMat[reference,] == 0, NA, xMat[reference,]))

  }
