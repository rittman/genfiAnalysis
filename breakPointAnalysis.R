library(bfast)
source("wholeBrainMetrics.R")

smoothData <- function(GS, dF, bw=10){
  if(!is.null(GS)){
    dF <- dF[dF$GS==GS,]
  } else {
    GS="all gene positive"
  }
  
  # resample data to enable time series analysis
  aa = approxfun(dF$aoo, dF$geNorm)
  dF.sm <- data.frame(x=seq(-50,50,by=.1), values=aa(seq(-50,50,by=0.1)))
  dF.sm <- dF.sm[complete.cases(dF.sm),]
  
  # smooth using a kernel with specified bandwidth
  dF.sm <- rbind(data.frame(dF.sm, type="orig"),
                 data.frame(x=dF.sm$x,
                            values=as.vector(ksmooth(dF.sm$x, dF.sm$values, bandwidth=bw))$y,
                            type="smooth")
                 )
  dF.sm <- data.frame(dF.sm, GS=GS)
  
  return(dF.sm)
}

dF <- importGraphData("geNorm", FALSE, 3)
# sp <- read.xlsx("../all_sm_thld10_SP.xlsx", sheetIndex = 1)
dF <- importGraphData(metric, weighted, edgePC)
# apply spike percentage threshold
dF <- applySP(dF, sp)
dF <- dF[dF$GS!="gene negative",]
dF$GS <- droplevels(dF$GS)

genfiData <- read.table("/home/tim/GENFI/genfi_Subjects_sjones_1_22_2015_17_47_47_restructure_summary.csv",
sep="\t",
header = TRUE)
dF.aoo <- data.frame(wbic=genfiData$Subject,
aoo=genfiData$Yrs.from.AV_AAO)
dF <- merge(dF, dF.aoo, by="wbic")

dF.sm <- do.call(rbind, lapply(levels(dF$GS), smoothData, dF=dF))
dF.sm <- rbind(dF.sm, smoothData(NULL, dF))

p <- ggplot(dF.sm, aes_string(x="x", y="values", colour="GS"))
p <- p + geom_line(aes_string(linetype="type"))
p <- p + geom_point(data=dF, aes_string(x="aoo", y="geNorm"), linetype="solid")


# new method of analysis
# create model
mod <- lm(geNorm~aoo, data=dF)
# do segmentation analysis assuming the breakpoint is at 0
my.seg <- segmented(mod, seg.Z = ~aoo, psi = c(0))
# print summary
summary(my.seg)

