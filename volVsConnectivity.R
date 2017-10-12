# Script to generate images of volumes vs connectivity fore the whole brain, frontal and temporal regions
library(grid)
library(gridExtra)

options(xtable.comment=FALSE)

# set the GENFI directory
genfiDir = "/media/tim/e87b22f3-4f3a-4dba-8c87-908464748c0d/backup/GENFI/GENFI_camgrid_20150525"

# get some useful functions for whole brain analysis
source("wholeBrainMetrics.R")

# another function used in the spike percentage stuff
getCorr <- function(n,x){
  ct <- cor.test(x[x$spMean<n,"spMean"], x[x$spMean<n,"values"])
  return(list(ct[["statistic"]][[1]],
              ct[["p.value"]][[1]])
  )
}

# function to plot multiple plots with a single legend
grid_arrange_shared_legend <- function(pp, pq, pc){
  # plots <- list(...)
  g <- ggplotGrob(pp + theme(legend.position="right"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)*2
  grid.arrange(
    arrangeGrob(pp + theme(legend.position="none")),
    arrangeGrob(pq + theme(legend.position="none")),
    arrangeGrob(pc + theme(legend.position="none")),
    legend,
    ncol=2,
    #         nrow=1,
    widths = c(1/2,1/2),
    heights = c(1/2,1/2))# c(2/7, 2/7, 2/7, 1/7))   #unit.c(unit(1, "npc") - lheight, lheight))
}

# function to plot multiple plots with a single legend
grid_arrange_shared_legend_generic <- function(...){
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="right"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x) x + theme(legend.position="none"))),
    legend,
    ncol=2,
    #         nrow=1,
    widths = c(3/4, 1/4) #unit.c(unit(1, "npc") - lheight, lheight)
  )
}

# set the colours to use in plots
cols <- list("gene negative"="#E69F00",
             "gene carriers"="#56B4E9",
             "FTD"="#D55E00",
             "estimates"="#000000"
)


dF <- importGraphData("degree", weighted=TRUE)
dF <- stackIt(dF, "degree")

# Summarise the patient data
ptSum = ddply(dF, .(gene), summarise,
              "gene negative" = length(GS[GS=="gene negative"]),
              "gene carriers" = length(GS[GS=="gene carriers"]),
              "FTD"      = length(GS[GS=="FTD"])
)
ptSum <- data.frame(ptSum, Totals=rowSums(ptSum[,-1]))
tmpdF <- data.frame(gene="Totals",
                    as.data.frame(matrix(colSums(ptSum[,-1]), nrow=1))
)
names(tmpdF) <- names(ptSum)
ptSum <- rbind(ptSum, tmpdF)

# print(xtable(ptSum,
#              caption="The total number of subjects with functional imaging",
#              digits = c(0,0,0,0,0,0),
#              display = c("s", "s", "d", "d", "d","d")),
#       # type="html",
#       include.rownames=FALSE)


# import the spike percentage data
sp = read.table("../all_sm_thld10_SP.csv", header = TRUE, sep=",")

sp.sub <- data.frame(wbic=sp$id, spMean=sp$mean)
dF.sp <- merge(dF, sp.sub, by="wbic")

# set cut-offs for the correlation between spike percentage and connection strength
cutOffs <- seq(2,40, by=1)
ctVals <- sapply(cutOffs, getCorr, x=dF.sp)
dF.ctVals <- data.frame(matrix(unlist(ctVals), byrow = TRUE, ncol = length(ctVals[,1])))
names(dF.ctVals) <- c("r", "p")

# plot the results
dF.plot <- data.frame(cutOffs=cutOffs, dF.ctVals)

p <- ggplot(dF.plot, aes_string(x="cutOffs", y="p"))
p <- p + geom_line(colour="red")
p <- p + geom_hline(yintercept=0.05, linetype="dashed", colour="black")
p <- p + theme_bw()
p <- p + labs(y="p-value", x="spike percentage threshold")
p <- p + theme(text=element_text(size=8))
# plot(p)

# set spike percentage threshold
spVal = 10  # would be better to define this using a Bayesian model averaging method


# apply the spike percentage
dF <- applySP(dF, sp, spVal=spVal)

# Summarise the patient data
ptSum.sp = ddply(dF, .(gene), summarise,
                 "gene negative" = length(GS[GS=="gene negative"]),
                 "gene carriers" = length(GS[GS=="gene carriers"]),
                 "FTD"      = length(GS[GS=="FTD"])
)
ptSum.sp <- data.frame(ptSum.sp, Totals=rowSums(ptSum.sp[,-1]))
tmpdF <- data.frame(gene="Totals",
                    as.data.frame(matrix(colSums(ptSum.sp[,-1]), nrow=1))
)
names(tmpdF) <- names(ptSum.sp)
ptSum.sp <- rbind(ptSum.sp, tmpdF)

# print(xtable(ptSum.sp,
#              caption="Subjects included in the analysis",
#              digits = c(0,0,0,0,0,0),
#              display = c("s", "s", "d", "d", "d","d")),
#       # type="html",
#       include.rownames=FALSE)

scanTypeList <- list(site="category",
                     dim1="continuous",
                     dim2="continuous",
                     dim3="continuous",
                     nvols="continuous",
                     tr_to_use="continuous",
                     ET="continuous",
                     scan_duration="continuous",
                     manufacturer="category",
                     tesla="category" )

scanSumList <- names(scanTypeList)
scanSumTab <- lapply(scanSumList, scanSummaries, scanTypeList=scanTypeList, sp=sp)
nameList <- names(scanSumTab[[1]])
scanSumTab <- data.frame(matrix(unlist(scanSumTab), byrow = TRUE, nrow = length(scanSumTab)))
names(scanSumTab) <- nameList
rownames(scanSumTab) <- scanSumList

for(i in seq(6)){ # first 6 columns of dataframe
  scanSumTab[,i] = as.numeric(as.character(scanSumTab[,i]))
}

print(xtable(scanSumTab,
             caption="Scanning parameters",
             digits = c(0,2,2,2,2,2,2,0,0),
             display = c("s", "fg", "fg", "fg", "fg", "fg", "fg", "s", "s")),
      include.rownames=TRUE)

exclDF = cbind(data.frame(gene=ptSum[,1]),ptSum[,-1] - ptSum.sp[,-1])
# print(xtable(exclDF,
#              caption="Subjects excluded in the analysis",
#              digits = c(0,0,0,0,0,0),
#              display = c("s", "s", "d", "d", "d","d")),
#       include.rownames=FALSE)

# print(xtable(ptSum.sp,
#              caption="Remaining subjects included in the analysis",
#              digits = c(0,0,0,0,0,0),
#              display = c("s", "s", "d", "d", "d","d")),
#       include.rownames=FALSE)

metrics = list(#"ge"="global efficiency",
  "geNorm"="global efficiency",
  "leNorm"="local efficiency (normalised)",
  #                "pl"="path length",
  # "sw"="Small worldness",
  #                "eigCentNP"="eigen centrality",
  #                "eigCentNorm"="eigen centrality (normalised)",
  #                "betCent"="betweenness centrality",
  #                "betCentNorm"="betweenness centrality (normalised)",
  #                "closeCent"="closeness centrality",
  "sw"="small worldness",
  "closeCentNorm"="closeness centrality",
  "medianELNorm"="edge distance",
  "degree"="connection strength",
  "plNorm"="path length",
  "clusterCoeffNorm"="cluster coefficient",
  "elnNorm"="Edge length normalised"
)

lobeList = list(Frontal="Frontal",
                Temporal="Temporal",
                Parietal="Parietal",
                Occipital="Occipital",
                Cerebellum="Cerebellum",
                Hippocampus="Hippocampus",
                Cingulate="Cingulate",
                Insula="Insula",
                Subcortical="Subcortical")

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

scoreListCog = list(
  N1="MMSE",
  N2="LogMemoryImmediate",
  N3="LogMemoryDelayed",
  N4="DigitSpanForward",
  N5="DigitSpanBackward",
  N6="TrailsA",
  N7="TrailsB",
  N8="DigitSymbol",
  N9="BostonNaming",
  N10="VerbalFluencyCategory",
  N11="VerbalFluencyLetter",
  N12="BlockDesign"
)

# linear mixed effects model for all subjects
metric = names(metrics)[6]
metricName = metrics[[names(metrics)[6]]]

for(cs in names(scoreListVols)){
  csName = scoreListVols[[cs]]
  
  wbame = clinicalScores(metric, metricName, cs, paste(csName, "\n(mls)"), cols=cols, sp=sp, weighted=TRUE, tsz=10, exclNeg=TRUE)
  
  t1 = wbame[[1]]
  t2 = wbame[[2]]
  t3 = wbame[[3]]
  t4 = wbame[[4]]
  t5 = wbame[[5]]
  t6 = wbame[[6]]
  t7 = wbame[[7]]
  t8 = wbame[[8]]
  t9 = wbame[[9]]
  t10= wbame[[10]]
  t11= wbame[[11]]
  t12= wbame[[12]]
  t13= wbame[[13]]
  t14= wbame[[14]]
  p1 = wbame[[15]]
  p2 = wbame[[16]]
  p3 = wbame[[17]]
  p4 = wbame[[18]]
  
  # plots of random effects variance
  # grid.arrange(p1, p2, nrow=1)
  
  # plot(p4)
  
  # print(t1, include.rownames=FALSE)
  # print(t2, include.rownames=FALSE)
  # print(t3, include.rownames=FALSE)
  # print(t6, include.rownames=FALSE) # table of random effects
  # print(t7, include.rownames=TRUE) # table of fixed effects
  # print(t11, include.rownames=FALSE) # table of random effects - excluding gene interaction term
  # print(t12, include.rownames=TRUE) # table of fixed effects - excluding gene interaction term
  # 
  # model comparison
  # print(t13, include.rownames=TRUE) # table of fixed effects
  
  ggsave(paste("volVsConnectivity/volumesVsConnectivity_",csName,".png",sep=""), plot=p4, width=80, height=80, units="mm")
}
