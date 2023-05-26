novelTsComms <- function(timeseries,
                         timeIds,
                         alpha,
                         method,
                         modelType = "logNormal", # (or "beta")
                         abundWeight=FALSE,
                         dissWeight=NULL,
                         plot=TRUE,
                         gam.max.k = -1){

  require(vegan) # community disimilarity metrics
  require(mgcv) # additive modelling
  require(arm) # invlogit transformation

  if(alpha < 0 | alpha > 1){stop("Alpha must be a probability threshold where 0 < a < 1.")}

  if(class(timeseries)[1] != "matrix"){

    error.trap <- tryCatch(as.matrix(timeseries), error=function(e) e)
    if(inherits(error.trap, "error")){
      stop("Time series must be a matrix or matrix-like object")
    }

  timeseries <- as.matrix(timeseries)}

  # order matrix by timeIds
  timeseries <- timeseries[timeIds,]

  timeIdsNum <- as.numeric(timeIds)

  # generate instantaneous novelty
  if(method=="weighted"){
    instantNov <- suppressWarnings(turnover(timeseries, method, reference="pairwise",
                                            dissWeight = dissWeight,
                                            abundWeight = abundWeight),
                                            classes="warning")
  } else {
    instantNov <- suppressWarnings(turnover(timeseries, method, reference="pairwise"),
                                   classes="warning")
  }

  if(var(instantNov, na.rm=TRUE)==0){stop("There is zero variation in sequential novelty")}

  # generate cumulative novelty
  if(method=="weighted"){

  cumulNov <- sapply(1:nrow(timeseries), function(n){
    if(n==1){return(NA)}

    suppressWarnings(as.vector(novelty(timeseries[n,],
                                       refSet=matrix(timeseries[1:(n-1),],
                                                     ncol=ncol(timeseries),
                                                     dimnames=list(timeIds[1:(n-1)],
                                                                   colnames(timeseries))),
                                       method=method,
                                       dissWeight=dissWeight,
                                       abundWeight=abundWeight,
                                       nSize=1)), classes="warning")
  })

  } else {

    suppressWarnings(as.vector(novelty(timeseries[n,],
                                       refSet=matrix(timeseries[1:(n-1),],
                                                    ncol=ncol(timeseries),
                                                    dimnames=list(timeIds[1:(n-1)],
                                       colnames(timeseries))),
                                       method=method,
                                       dissWeight=dissWeight,
                                       abundWeight=abundWeight,
                                       nSize=1)), classes="warning")

  }
  if(var(cumulNov, na.rm=TRUE)==0){stop("There is zero variation in cumulative novelty")}

  # INSTANTANEOUS NOVELTY MODEL ####

  if(modelType == "beta"){
    # transform to remove 0s and 1s for beta regression
  instantNovtr <- (instantNov * (length(instantNov)-1) + 0.5) / length(instantNov)

  if(var(instantNov, na.rm=TRUE)==0){stop("There is zero variation in sequential novelty")}

  instantGam <- gam(instantNovtr ~ s(timeIdsNum, k = length(timeIds)-2),
                   family=betar(),
                   method="REML")

  instantMu <- c(NA, instantGam$fitted.values)
  instantPhi <- as.numeric(substr(instantGam$family$family,
                           regexpr("\\(", instantGam$family$family)+1,
                           nchar(instantGam$family$family)-1))
  instantA = instantMu * instantPhi
  instantB = instantPhi - instantA

  instantP <- do.call("rbind", lapply(1:length(instantA), function(n){
    data.frame(lwr=qbeta(0.05, shape1=instantA[n], shape2=instantB[n]),
               upr=qbeta(0.95, shape1=instantA[n], shape2=instantB[n]),
               p = pbeta(instantNov[n], shape1=instantA[n], shape2=instantB[n], lower.tail=FALSE))
  }))

  }

  if(sum(instantNov == 0 | cumulNov == 0, na.rm=TRUE) > 0 &
     modelType == "logNormal"){
    warning("There are dissimilarities of 0 which cannot be log-transformed")
  }

  if(modelType == "lognormal"){

    # log transform dissimilarities
    instantNovtr <- log(instantNov+1)

    if(var(instantNov, na.rm=TRUE)==0){stop("There is zero variation in instantaneous novelty")}

    instantGam <- gam(instantNovtr ~ s(timeIdsNum, k = length(timeIds)-2),
                      method="REML")

    instantMu <- instantGam$fitted.values

    predInt <- sqrt(sum(instantGam$residuals^2) / instantGam$df.residual) * sqrt(1 + instantGam$hat)

    tCrit <- abs(qt(0.025, instantGam$df.residual))

    instantP <- do.call("rbind", lapply(1:length(instantMu), function(n){
      data.frame(lwr = exp(instantMu[n] - (tCrit * predInt[n])),
                 upr = exp(instantMu[n] + (tCrit * predInt[n])),
                 p =  pnorm(instantNovtr[n+1], mean=instantMu[n], sd=predInt[n], lower.tail=FALSE))
    }))
    instantP <- rbind(c(NA, NA, NA),
                      instantP)
    instantMu <- c(NA, instantMu)
  }

  # CUMULATIVE DISIMILARITY MODEL ####

  if(modelType == "beta"){
  # transform to remove 0s and 1s for beta regression
  cumulNovTr <- (cumulNov * (length(cumulNov)-1) + 0.5) / length(cumulNov)

  cumulGam <- gam(cumulNovTr ~ s(timeIdsNum, k = length(timeIds)-2),
                    family=betar(),
                    method="REML")

  cumulMu <- c(NA, cumulGam$fitted.values)
  cumulPhi <- as.numeric(substr(cumulGam$family$family,
                                  regexpr("\\(", cumulGam$family$family)+1,
                                  nchar(cumulGam$family$family)-1))
  cumulA = cumulMu * cumulPhi
  cumulB = cumulPhi - cumulA

  cumulP <- do.call("rbind", lapply(1:length(cumulA), function(n){
    data.frame(lwr=qbeta(0.05, shape1=cumulA[n], shape2=cumulB[n]),
               upr=qbeta(0.95, shape1=cumulA[n], shape2=cumulB[n]),
               p = pbeta(cumulNov[n], shape1=cumulA[n], shape2=cumulB[n], lower.tail=FALSE))
  }))

  }

  if(modelType=="lognormal"){

    # log transform dissimilarities
    cumulNovtr <- log(cumulNov+1)

    if(var(cumulNov, na.rm=TRUE)==0){stop("There is zero variation in cumulative novelty")}

    cumulGam <- gam(cumulNovtr ~ s(timeIdsNum, k = length(timeIds)-2),
                      method="REML")

    cumulMu <- cumulGam$fitted.values

    predInt <- sqrt(sum(cumulGam$residuals^2) / cumulGam$df.residual) * sqrt(1 + cumulGam$hat)

    tCrit <- abs(qt(0.025, cumulGam$df.residual))

    cumulP <- do.call("rbind", lapply(1:length(cumulMu), function(n){
      data.frame(lwr = exp(cumulMu[n] - (tCrit * predInt[n])),
                 upr = exp(cumulMu[n] + (tCrit * predInt[n])),
                 p =  pnorm(cumulNovtr[n+1], mean=cumulMu[n], sd=predInt[n], lower.tail=FALSE))
    }))
    cumulP <- rbind(c(NA, NA, NA),
                      cumulP)
    cumulMu <- c(NA, cumulMu)
  }

  # CREATE RETURN DATA-FRAME ####
  return.data <- data.frame(timeIds = timeIds,
                            timeLag = c(NA, abs(diff(timeIdsNum))),
                            instantNov = instantNov,
                            cumulNov = cumulNov,
                            instantP = instantP$p,
                            cumulP = cumulP$p,
                            instantExpected = instantMu,
                            cumulExpected = cumulMu,
                            instantResid = c(NA, resid(instantGam)),
                            cumulResid = c(NA, resid(cumulGam)),
                            isInstantNov = instantP$p <= alpha,
                            isCumulNov = cumulP$p <= alpha,
                            isNovelComm = instantP$p <= alpha & cumulP$p <= alpha)

  return.data$cat = "backComm"
  return.data[return.data$isInstantNov &
                !is.na(return.data$instantNov) &
                !is.na(return.data$cumulNov) &
                !return.data$cumulNov, "cat"] = "instantNov"
  return.data[!return.data$isInstantNov &
                !is.na(return.data$isInstantNov) &
                !is.na(return.data$isCumulNov) &
                return.data$isCumulNov, "cat"] = "cumulNov"
  return.data[return.data$isInstantNov &
                !is.na(return.data$isCumulNov) &
                !is.na(return.data$isInstantNov) &
                return.data$isCumulNov, "cat"] = "novelComm"

  # PLOT TIME SERIES ####

  # this section of code generates the plot seen when running the function.
  if(plot){
    par(mfrow=c(2,1), mar=c(0,0,0,0), oma=c(3,6,0.5,9), las=1)

    plot(instantNov ~ timeIdsNum, type="n",
         ylim=c(max(instantNov, na.rm=TRUE)+0.1,
                min(instantNov, na.rm=TRUE)),
         axes=FALSE, xlab="", ylab="")

    polygon(x=c(timeIdsNum, rev(timeIdsNum)),
            y=c(instantP[,1], rev(instantP[,2])),
            col="grey80", border=NA)

    if(modelType == "beta"){
    lines(plogis(predict(instantGam)) ~ timeIdsNum[-1], col="grey50", lwd=2)
      } else {
        lines(exp(predict(instantGam)) ~ timeIdsNum[-1], col="grey50", lwd=2)
      }
    lines(instantNov ~ timeIdsNum)
    points(y=instantNov[instantP$p< 0.05],
           x=timeIdsNum[instantP$p < 0.05], pch=21, bg="red")
    lims <- par("usr")

    sapply(which(return.data$isNovelComm), function(x){
      segments(x0 = timeIdsNum[x],
               x1 = timeIdsNum[x],
               y0 = instantNov[x] + (0.05 * (par("usr")[3] - par("usr")[4])),
               y1 = par("usr")[3], col="orange", lwd=1)
    })

    segments(x0=par("usr")[1], x1=par("usr")[1], y0=par("usr")[3], y1=par("usr")[4])
    segments(x0=par("usr")[1], x1=par("usr")[2], y0=par("usr")[4], y1=par("usr")[4])
    segments(x0=par("usr")[2], x1=par("usr")[2], y0=par("usr")[3], y1=par("usr")[4])
    axis(side=2)
    mtext(side=2, text = "Sequential\nnovelty", line=3.5, las=0)

    plot(cumulNov ~ timeIdsNum, type="n", axes=FALSE,
         ylim=c(min(cumulNov, na.rm=TRUE),
                max(cumulNov, na.rm=TRUE)+0.1), xlab="", ylab="")
    polygon(x=c(timeIdsNum, rev(timeIdsNum)),
            y=c(cumulP[,1], rev(cumulP[,2])),
            col="grey80", border=NA)

    if(modelType=="beta"){
    lines(plogis(predict(cumulGam)) ~ timeIdsNum[-1], col="grey50", lwd=2)
    } else {
      lines(exp(predict(cumulGam)) ~ timeIdsNum[-1], col="grey50", lwd=2)
    }
    lines(cumulNov ~ timeIdsNum)

    points(y=cumulNov[cumulP$p< 0.05],
           x=timeIds[cumulP$p < 0.05], pch=21, bg="skyblue")

    sapply(which(return.data$isNovelComm), function(x){
      segments(x0 = timeIdsNum[x],
               x1 = timeIdsNum[x],
               y0 = cumulNov[x] + (0.05 * (par("usr")[4] - par("usr")[3])),
               y1 =par("usr")[4], col="orange", lwd=1)
    })

    par(xpd=NA)
    points(y=rep(par("usr")[4], sum(cumulP$p< 0.05 & instantP$p < 0.05, na.rm=TRUE)),
           x=timeIds[cumulP$p< 0.05 & instantP$p < 0.05 & !is.na(cumulP$p)],
           pch=21, bg="orange")
    par(xpd=TRUE)

    segments(x0=par("usr")[1], x1=par("usr")[1], y0=par("usr")[3], y1=par("usr")[4])
    segments(x0=par("usr")[1], x1=par("usr")[2], y0=par("usr")[3], y1=par("usr")[3])
    segments(x0=par("usr")[2], x1=par("usr")[2], y0=par("usr")[3], y1=par("usr")[4])

    axis(side=1)
    mtext(side=1, "Time Ids", line=2)

    axis(side=2)
    mtext(side=2, text = "Cumulative\nnovelty", line=3.5, las=0)

    par(xpd=NA)
    legend(x=par("usr")[2], y= par("usr")[4],
           legend=c("", "", ""),
           pch=c(NA,NA,NA), lwd=c(NA,NA,1),
           pt.bg=c("red","skyblue","orange"), col=c("black","black","orange"),
           yjust=0.5, xjust=0, bty="n", seg.len=1, x.intersp=0.5)

    legend(x=par("usr")[2], y= par("usr")[4],
           legend=c("Instantaneous", "Cumulative", "Novel"),
           pch=c(21,21,21), lwd=c(NA,NA,NA),
           pt.bg=c("red","skyblue","orange"), col=c("black","black","black"),
           yjust=0.5, xjust=0, bty="n", seg.len=1, x.intersp=0.5)
    par(xpd=FALSE)
  }

  return(return.data)

}
