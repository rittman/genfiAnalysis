# Script to plot volume and network values on the same axes to try and demonstrate structural/functional dissociation

getlinesVals <- function(piecewise, dF, resample=TRUE){
  y <- piecewise$coefficients[[1]]
  piecewise.CIs <- confint(piecewise)
  
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
  
  if(resample){
    ygap = yMin - yMax
    yEndVal = (yEndVal - yMax)/ygap
    yVal = (y - yMax)/ygap
    yMin = 1.
    yMax = 0.
    yVal.lower = y.lower
    yVal.upper = y.upper
  } else {
    yVal = y
    yVal.lower = y.lower
    yVal.upper = y.upper
  }

  return(list(yMin=yMin, yMax = yMax, yVal=yVal, yVal.lower=yVal.lower, yVal.upper=yVal.upper, yEndVal=yEndVal))
}

scoreListVols = list(
# 								TIV="TIV",
# 								WBV="WBV",
WBV_C="Whole Brain Volume",  # whole brain value as percentage of total intracranial volume
T_FR_C="Frontal lobe volume",
T_TEM_C="Temporal lobe Volume",
T_PAR_C="Parietal lobe Volume",
T_OCC_C="Occipital lobe Volume",
T_CIN_C="Cingulate volume",
T_INS_C="Insula volume",
T_HIPPO_C="Hippocampal volume",
T_AMYG_C="Amygdala volume",
T_CAUD_C="Caudate volume",
T_PUT_C="Putamen volume",
T_THAL_C="Thalamic volume"
)

plotGMvsVol <- function(cs, metric, weighted, sp){
  ## connectivity vs volume
  # Import data
  dF <- importGraphData(metric, weighted)
  dF <- applySP(dF, sp)
  dF.wb <- stackIt(dF, metric)
  
  csName = scoreListVols[[cs]]
  csdF = read.table("../clinicalScores.csv", sep=",", header = TRUE)
  names(csdF)[which(names(csdF)=="BLINDID")] <- "wbic"
  names(csdF)[which(names(csdF)=="Yrs.from.AV_AAO")] <- "aoo"
  csdF <- csdF[,c("wbic",cs,"aoo")]
  names(csdF)[which(names(csdF)==cs)] = "score"
  # merge network and imaging data
  dF.wb <- merge(dF.wb, csdF,
                  by = c("wbic"),
                  all.x=TRUE,
                  all.y=FALSE)
  
  # split gene positive and gene negative groups
  dFgn <- dF.wb[dF.wb$GS=="gene negative",]
  dF.wb <- dF.wb[dF.wb$GS!="gene negative",]
  dF.wb$GS <- factor(as.character(dF.wb$GS), levels=c("gene carriers", "FTD"))
  
  # calculate z scores
  sd.geneNeg.values <- sd(dFgn[,"values"], na.rm=TRUE)
  mean.geneNeg.values <- mean(dFgn[,"values"], na.rm=TRUE)
  dF.wb <- cbind(dF.wb, values.norm=sapply(dF.wb$values, function(x) (x-mean.geneNeg.values)/sd.geneNeg.values))
  
  sd.geneNeg.score <- sd(dFgn[,"score"], na.rm=TRUE)
  mean.geneNeg.score <- mean(dFgn[,"score"], na.rm=TRUE)
  dF.wb <- cbind(dF.wb, score.norm=sapply(dF.wb$score, function(x) (x-mean.geneNeg.score)/sd.geneNeg.score))
  
  # Do piecewise analysis
  piecewise.vals <- lm(values.norm ~ aoo*(aoo<bkpt) + aoo*(aoo>=bkpt), data=dF.wb)
  piecewise.score <- lm(score.norm ~ aoo*(aoo<bkpt) + aoo*(aoo>=bkpt), data=dF.wb)
  
  lineVals <- getlinesVals(piecewise.vals, dF.wb, resample = TRUE)
  
  yMin.vals = lineVals[["yMin"]]
  yMax.vals = lineVals[["yMax"]]
  yVal.vals = lineVals[["yVal"]]
  yVal.vals.lower = lineVals[["yVal.lower"]]
  yVal.vals.upper = lineVals[["yVal.upper"]]
  yEndVal.vals = lineVals[["yEndVal"]]
  
  p <- ggplot(dF.wb, aes_string(x="aoo", y="values"))
  
  # p <- p + geom_point(size=ps, shape=16, aes_string(colour="GS"))
  p <- p + geom_segment(x=min(dF.wb$aoo),
                        xend=0.,
                        y=yMin.vals,
                        yend=yEndVal.vals,
                        size=1.5,
                        colour="black")
  
  p <- p + geom_segment(x=0.,
                        xend=max(dF.wb$aoo),
                        y=yVal.vals,
                        yend=yMax.vals,
                        size=1.5,
                        colour="black")
  
  # # add confidence intervals
  # p <- p + geom_smooth(method="lm",
  #                      linetype="blank",
  #                      data=dF.wb[dF.wb$aoo<=0.,],
  #                      fill=colList["gene carriers"], alpha=0.5)
  # 
  # p <- p + geom_smooth(method="lm",
  #                      linetype="blank",
  #                      data=dF.wb[dF.wb$aoo>0.,],
  #                      fill=colList["FTD"], alpha=0.5)
  
  lineVals <- getlinesVals(piecewise.score, dF.wb, resample = TRUE)
  
  yMin.score = lineVals[["yMin"]]
  yMax.score = lineVals[["yMax"]]
  yVal.score = lineVals[["yVal"]]
  yVal.score.lower = lineVals[["yVal.lower"]]
  yVal.score.upper = lineVals[["yVal.upper"]]
  yEndVal.score = lineVals[["yEndVal"]]
  
  p <- p + geom_segment(x=min(dF.wb$aoo),
                        xend=0.,
                        y=yMin.score,
                        yend=yEndVal.score,
                        size=1.5,
                        colour="red")
  
  p <- p + geom_segment(x=0.,
                        xend=max(dF.wb$aoo),
                        y=yVal.score,
                        yend=yMax.score,
                        size=1.5,
                        colour="red")
  
  p <- p + geom_blank()
  
  ymin = min(yMin.vals, yMin.score, yMax.vals, yMax.score, yEndVal.vals, yEndVal.score)
  ymax = max(yMin.vals, yMin.score, yMax.vals, yMax.score, yEndVal.vals, yEndVal.score)
  p <- p + scale_y_continuous(limits=c(ymin, ymax))
  
  return(p)
}

# import the spike percentage data
sp = read.table("../all_sm_thld10_SP.csv", header = TRUE, sep=",")

## Connectivity
p.conn <- plotGMvsVol("WBV_C", "degree", weighted=TRUE, sp=sp)

## Global efficiency
p.ge <- plotGMvsVol("WBV_C", "ge", weighted=FALSE, sp=sp)



## plot them all together
dF.deg <- importGraphData("degree", TRUE)
dF.deg <- applySP(dF.deg, sp)
dF.wb.deg <- stackIt(dF.deg, "degree")

dF.ge <- importGraphData("ge", FALSE)
dF.ge <- applySP(dF.ge, sp)
dF.wb.ge <- stackIt(dF.ge, "ge")

csName = scoreListVols[[cs]]
csdF = read.table("../clinicalScores.csv", sep=",", header = TRUE)
names(csdF)[which(names(csdF)=="BLINDID")] <- "wbic"
names(csdF)[which(names(csdF)=="Yrs.from.AV_AAO")] <- "aoo"
csdF <- csdF[,c("wbic",cs,"aoo")]
names(csdF)[which(names(csdF)==cs)] = "score"

# merge network and imaging data
dF.wb.deg <- merge(dF.wb.deg, csdF,
               by = c("wbic"),
               all.x=TRUE,
               all.y=FALSE)

dF.wb.ge <- merge(dF.wb.ge, csdF,
                  by = c("wbic"),
                  all.x=TRUE,
                  all.y=FALSE)

dF.wb <- cbind(dF.wb.deg, ge=dF.wb.ge$values)
names(dF.wb)[which(names(dF.wb)=="values")] <- "degree"

# split gene positive and gene negative groups
dFgn <- dF.wb[dF.wb$GS=="gene negative",]
dF.wb <- dF.wb[dF.wb$GS!="gene negative",]
dF.wb$GS <- factor(as.character(dF.wb$GS), levels=c("gene carriers", "FTD"))

# calculate z scores
sd.geneNeg.degree <- sd(dFgn[,"degree"], na.rm=TRUE)
mean.geneNeg.degree <- mean(dFgn[,"degree"], na.rm=TRUE)
dF.wb <- cbind(dF.wb, degree.norm=sapply(dF.wb$degree, function(x) (x-mean.geneNeg.degree)/sd.geneNeg.degree))

sd.geneNeg.ge <- sd(dFgn[,"ge"], na.rm=TRUE)
mean.geneNeg.ge <- mean(dFgn[,"ge"], na.rm=TRUE)
dF.wb <- cbind(dF.wb, ge.norm=sapply(dF.wb$ge, function(x) (x-mean.geneNeg.ge)/sd.geneNeg.ge))

sd.geneNeg.score <- sd(dFgn[,"score"], na.rm=TRUE)
mean.geneNeg.score <- mean(dFgn[,"score"], na.rm=TRUE)
dF.wb <- cbind(dF.wb, score.norm=sapply(dF.wb$score, function(x) (x-mean.geneNeg.score)/sd.geneNeg.score))

# Do piecewise analysis
piecewise.degree <- lm(degree.norm ~ aoo*(aoo<bkpt) + aoo*(aoo>=bkpt), data=dF.wb)
piecewise.ge <- lm(ge.norm ~ aoo*(aoo<bkpt) + aoo*(aoo>=bkpt), data=dF.wb)
piecewise.score <- lm(score.norm ~ aoo*(aoo<bkpt) + aoo*(aoo>=bkpt), data=dF.wb)

# initiate plot
p <- ggplot(dF.wb, aes_string(x="aoo", y="degree"))

lineVals.degree <- getlinesVals(piecewise.degree, dF.wb, resample = TRUE)

yMin.degree = lineVals.degree[["yMin"]]
yMax.degree = lineVals.degree[["yMax"]]
yVal.degree = lineVals.degree[["yVal"]]
yVal.degree.lower = lineVals.degree[["yVal.lower"]]
yVal.degree.upper = lineVals.degree[["yVal.upper"]]
yEndVal.degree = lineVals.degree[["yEndVal"]]

# p <- p + geom_point(size=ps, shape=16, aes_string(colour="GS"))
p <- p + geom_segment(x=min(dF.wb$aoo),
                      xend=0.,
                      y=yMin.degree,
                      yend=yEndVal.degree,
                      size=1.5,
                      colour="black")

p <- p + geom_segment(x=0.,
                      xend=max(dF.wb$aoo),
                      y=yVal.degree,
                      yend=yMax.degree,
                      size=1.5,
                      colour="black")

# # add confidence interdegree
# p <- p + geom_smooth(method="lm",
#                      linetype="blank",
#                      data=dF.wb[dF.wb$aoo<=0.,],
#                      fill=colList["gene carriers"], alpha=0.5)
# 
# p <- p + geom_smooth(method="lm",
#                      linetype="blank",
#                      data=dF.wb[dF.wb$aoo>0.,],
#                      fill=colList["FTD"], alpha=0.5)

lineVals.ge <- getlinesVals(piecewise.ge, dF.wb, resample = TRUE)

yMin.ge = lineVals.ge[["yMin"]]
yMax.ge = lineVals.ge[["yMax"]]
yVal.ge = lineVals.ge[["yVal"]]
yVal.ge.lower = lineVals.ge[["yVal.lower"]]
yVal.ge.upper = lineVals.ge[["yVal.upper"]]
yEndVal.ge = lineVals.ge[["yEndVal"]]

# p <- p + geom_point(size=ps, shape=16, aes_string(colour="GS"))
p <- p + geom_segment(x=min(dF.wb$aoo),
                      xend=0.,
                      y=yMin.ge,
                      yend=yEndVal.ge,
                      size=1.5,
                      colour="black")

p <- p + geom_segment(x=0.,
                      xend=max(dF.wb$aoo),
                      y=yVal.ge,
                      yend=yMax.ge,
                      size=1.5,
                      colour="black")

lineVals <- getlinesVals(piecewise.score, dF.wb, resample = TRUE)

yMin.score = lineVals[["yMin"]]
yMax.score = lineVals[["yMax"]]
yVal.score = lineVals[["yVal"]]
yVal.score.lower = lineVals[["yVal.lower"]]
yVal.score.upper = lineVals[["yVal.upper"]]
yEndVal.score = lineVals[["yEndVal"]]

p <- p + geom_segment(x=min(dF.wb$aoo),
                      xend=0.,
                      y=yMin.score,
                      yend=yEndVal.score,
                      size=1.5,
                      colour="red")

p <- p + geom_segment(x=0.,
                      xend=max(dF.wb$aoo),
                      y=yVal.score,
                      yend=yMax.score,
                      size=1.5,
                      colour="red")

p <- p + geom_blank()

p <- p + scale_y_continuous(limits=c(0,1))
