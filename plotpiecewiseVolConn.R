getlinesVals <- function(piecewise, dF, dFgn, normalise=FALSE){
  y <- piecewise$coefficients[[1]]
  piecewise.CIs <- confint(piecewise)
  
  # control linear model
  gnlm <- lm(values ~ aoo, data=dFgn)
  int <- gnlm[["coefficients"]][[1]]
  slope <- gnlm[["coefficients"]][[2]]
  sd.geneNeg <- sd(dFgn[,"values"], na.rm=TRUE)
  mean.geneNeg <- mean(dFgn[,"values"], na.rm=TRUE)
  
  # calculate points to plot
  yMin = min(dF$aoo)*(piecewise$coefficients[[2]]+piecewise$coefficients[[5]])+y+piecewise$coefficients[[3]]
  yMax = max(dF$aoo)*piecewise$coefficients[[2]]+piecewise$coefficients[[1]]
  yEndVal = y+piecewise$coefficients[[3]]
  
  y.lower <- piecewise.CIs[1,1]
  yMin.lower <- min(dF$aoo)*(piecewise.CIs[2,1]+piecewise.CIs[5,1])+y+piecewise.CIs[3,1]
  yMax.lower <- max(dF$aoo)*piecewise.CIs[2,1]+piecewise.CIs[1,1]
  yEndVal.lower <- y+max(piecewise.CIs[3,2],piecewise.CIs[3,1])
  
  y.upper <- piecewise.CIs[1,2]
  yMin.upper <- min(dF$aoo)*(piecewise.CIs[2,2]+piecewise.CIs[5,2])+y+piecewise.CIs[3,2]
  yMax.upper <- max(dF$aoo)*piecewise.CIs[2,2]+piecewise.CIs[1,2]
  yEndVal.upper <- y+min(piecewise.CIs[3,2],piecewise.CIs[3,1])
  
  if(normalise){
    # calculate gene negative values
    y.GeneNeg <- gnlm$coefficients[[1]] 
    yMin.GeneNeg <- min(dF$aoo)*gnlm$coefficients[[2]]+gnlm$coefficients[[1]]
    yMax.GeneNeg <- max(dF$aoo)*gnlm$coefficients[[2]]+gnlm$coefficients[[1]]
    yEndVal.GeneNeg <- y.GeneNeg
    
    yMin = (yMin - yMin.GeneNeg)/sd.geneNeg
    yMax = (yMax - yMax.GeneNeg)/sd.geneNeg
    yEndVal = (yEndVal - yEndVal.GeneNeg)/sd.geneNeg
    yVal = (y-y.GeneNeg)/sd.geneNeg
    
    yMin.lower = (yMin.lower - yMin.GeneNeg) / sd.geneNeg
    yMax.lower = (yMax.lower - yMax.GeneNeg) / sd.geneNeg
    yEndVal.lower = (yEndVal.lower - yEndVal.GeneNeg) / sd.geneNeg
    yVal.lower = (y.lower-y.GeneNeg)/sd.geneNeg
    
    yMin.upper = (yMin.upper - yMin.GeneNeg) / sd.geneNeg
    yMax.upper = (yMax.upper - yMax.GeneNeg) / sd.geneNeg
    yEndVal.upper = (yEndVal.upper - yEndVal.GeneNeg) / sd.geneNeg
    yVal.upper = (y.upper-y.GeneNeg)/sd.geneNeg
  } else {
    yVal = y
    yVal.lower = y.lower
    yVal.upper = y.upper
  }

  return(list(yMin=yMin, yMax = yMax, yVal=yVal, yVal.lower=yVal.lower, yVal.upper=yVal.upper, yEndVal=yEndVal))
}


yMin = lineVals[["yMin"]]
yMax = lineVals[["yMax"]]
yVal = lineVals[["yVal"]]
yVal.lower = lineVals[["yVal.lower"]]
yVal.uppwer = lineVals[["yVal.upper"]]
yEndVal = lineVals[["yEndVal"]]

# p <- p + geom_point(size=ps, shape=16, aes_string(colour="GS"))
p <- p + geom_segment(x=min(dF.wb$aoo),
                      xend=0.,
                      y=yMin,
                      yend=yEndVal,
                      size=1.5,
                      colour="black")

p <- p + geom_segment(x=0.,
                      xend=max(dF.wb$aoo),
                      y=yVal,
                      yend=yMax,
                      size=1.5,
                      colour="black")

# add confidence intervals
p <- p + geom_smooth(method="lm",
                     linetype="blank",
                     data=dF.wb[dF.wb$aoo<=0.,],
                     fill=colList["gene carriers"], alpha=0.5)

p <- p + geom_smooth(method="lm",
                     linetype="blank",
                     data=dF.wb[dF.wb$aoo>0.,],
                     fill=colList["FTD"], alpha=0.5)

