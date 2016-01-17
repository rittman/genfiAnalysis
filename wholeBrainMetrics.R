# This script will examine whole brain graph metrics in the GENFI dataset by looking for differences
# between affected/carriers/non-carriers, and by examining the relationship between average age at
# onset and the graph measures.

library(plyr)
library(ggplot2)
library(xtable)
library(car)
library(xlsx)
library(lme4)
library(segmented)

genfiDir = "/home/tim/GENFI/GENFI_camgrid_20150525/"

#### helper functions ####
dotTests <- function(gslist, dF){
  tRes <- t.test(dF[dF$GS==gslist[[1]],"values"], dF[dF$GS==gslist[[2]],"values"])
  return(list(paste(gslist, collapse = "/"),
              tRes$statistic,
              tRes$p.value))
}

# formatting function for numbers
fn <- function(x,a=1,b=2){
  if(x>=10){
    return(round(x,digits=a))
  } else {
    return(signif(x, digits=b))
  }
}

importGraphData <- function(metric, weighted, edgePC=3){
  # define input file
  inFile = paste("d2",metric,"local",sep="_")
  
  # import data
  dF = read.table(paste(genfiDir,inFile, sep="/"), header = TRUE, na.strings = "NA")
  
  # select only the desired percentage edge threshold if unweighted metric
  if(!weighted){
    dF <- dF[dF$edgePC==edgePC,]
  }
  
  # convert diagnostic label to a factor and rename
  dF$GS = as.factor(dF$GS)
  dF$GS = revalue(dF$GS, c("0" = "gene negative", "1"="gene positive", "2"="affected"))
  dF$site = as.factor(dF$site)
  
  # return dataframe with graph metrics
  return(dF)
}

applySP <- function(dF, sp, spVal=10){
  # combine data with spike percentage data
  sp.sub <- data.frame(wbic=sp$id, spMean=sp$mean)
  dF <- merge(dF, sp.sub, by="wbic")
  
  # filter by spike percentage
  dF <- dF[dF$spMean < spVal,]
  
  # return filtered data
  return(dF)
}

stackIt <- function(dF, metric){
  # Sort out stacking the nodes up if nodewise measure
  nodeNames = names(dF)[sapply(names(dF), function(x) grepl("X",x))]
  
  if(length(nodeNames)>0){
    # Take the mean values of
    dF.stacked = stack(dF[,nodeNames])
    dF.stacked = data.frame(dF.stacked, GS=dF$GS, gene=dF$gene, wbic=dF$wbic, site=dF$site, Family=dF$Family)
    
    dF.wb = ddply(dF.stacked, .(gene, wbic, GS, site, Family), summarise,
                  values = mean(values, na.rm = TRUE)
    )
  } else {
    dF.wb <- dF
    names(dF.wb)[names(dF.wb)==metric] <- "values"
  }
  
  return(dF.wb)
}

initiateLog <- function(outFile, metricName){
  header = c("\\documentclass[a4paper,10pt]{article}",
             "\\usepackage[utf8]{inputenc}",
             paste("\\title{GENFI data,",metricName,"}"),
             "\\author{Timothy Rittman}",
             "\\date{}",
             "\\pdfinfo{%",
             "/Title    ()",
             "/Author   ()",
             "/Creator  ()",
             "/Producer ()",
             "/Subject  ()",
             "/Keywords ()",
             "}",
             "\\begin{document}",
             "\\maketitle")
  
  write(paste(header, collapse="\n"),
        file=outFile,
        append=FALSE)
}

endLog <- function(outFile){
  write("\\end{document}", file=outFile, append=TRUE)
}
set.seed(1)
SSquadFun <- selfStart(~aq + bq*x + cq*x^2,
                       function(mCall, data, LHS){
                         xy <- sortedXyData(mCall[["x"]], LHS, data)
                         pars <- as.vector(coef(nls(y ~ aq + bq*x + cq*x^2, data=xy)))
                         setNames(c(pars[3], pars[1], pars[2]),
                                  mCall[c("aq", "bq", "cq")])
                         
                       }, c("aq", "bq", "cq"))

SScubicFun <- selfStart(~ac + bc*x + cc*x^2 + dc*x^3,
                       function(mCall, data, LHS){
                         xy <- sortedXyData(mCall[["x"]], LHS, data)
                         pars <- as.vector(coef(nls(y ~ ac + bc*x + cc*x^2 + dc*x^3, data=xy)))
                         setNames(c(pars[3], pars[1], pars[2]),
                                  mCall[c("ac", "bc", "cc", "dc")])
                         
                       }, c("ac", "bc", "cc", "dc"))

cubicFun <- function(x, ac,bc,cc,dd){
  # an attempt at building a quadratic function
  y = sapply(x, function(x) ac + bc*x + cc*x^2 + dc*x^3)
  return(y)
}

addSigBar <- function(p, pVal, d1, d2, xList, sig.value=0.05){
  if(pVal < sig.value){
    ymax = ggplot_build(p)$panel$ranges[[1]]$y.range[2]
    ymin = ggplot_build(p)$panel$ranges[[1]]$y.range[1]
    ygap = ymax - ymin
    
    if(pVal < 0.0001){
      pVal.plot = "****"
    } else if(pVal < 0.001) {
      pVal.plot = "***"
    } else if(pVal < 0.01) {
      pVal.plot = "**"
    } else {
      pVal.plot = "*"
    }
    ymax = ymax + .15*ygap
    
    ypos = ymax - ygap*.05
    xpos1 = xList[d1]
    xpos2 = xList[d2]
    p <- p + geom_segment(x=xpos1, xend=xpos2, y=ypos, yend=ypos, colour="black")
    p <- p + annotate("text", x=mean(c(xpos1, xpos2)), y=ypos, label=pVal.plot, colour="black", size=8)
    p <- p + scale_y_continuous(limits=c(ymin, ymax))
  }
  return(p)
}

addSigBarGenes <- function(p, pVal, d1, d2, xList, ymax, ng, sig.value=0.05){
  if(pVal < sig.value){
    ymin = ggplot_build(p)$panel$ranges[[1]]$y.range[1]
    ygap = ymax - ymin
    
    if(pVal < 0.0001){
      pVal.plot = "****"
    } else if(pVal < 0.001) {
      pVal.plot = "***"
    } else if(pVal < 0.01) {
      pVal.plot = "**"
    } else {
      pVal.plot = "*"
    }
    ymax = ymax + .15*ygap
    
    ypos = ymax - ygap*.05
    xpos1 = xList[d1]/3 + (ng-1) + 1/3 - (xList[d1]-2)*.08
    
    xpos2 = xList[d2]/3 + (ng-1) + 1/3 - (xList[d2]-2)*.08
    
    p <- p + geom_segment(x=xpos1, xend=xpos2, y=ypos, yend=ypos, colour="black")
    p <- p + annotate("text", x=mean(c(xpos1, xpos2)), y=ypos, label=pVal.plot, colour="black", size=8)
  }
  return(list(p=p, ymax=ymax))
}


#### main functions ####
wholeBrainAnalysis <- function(metric,
                               metricName,
                               cols,
                               sp, weighted=TRUE,
                               outDir="wholebrainResults",
                               edgePC=3, ts=12){
  
  if(weighted){
    metric = paste(metric,"wt",sep="_")
  }
  
  # create output directory
  dir.create(outDir, showWarnings = FALSE)
  
  # define log output file
  outFile = paste(outDir,
                  paste(metric,"logFile.tex",sep="_"),
                  sep="/")
  
  # create log file
  initiateLog(outFile, metricName)
  
  # import graph data
  dF <- importGraphData(metric, weighted, edgePC)
  
  # apply spike percentage threshold
  dF <- applySP(dF, sp)
  
  colList = unlist(cols[levels(dF$GS)])
  
  # Summarise the patient data
  ptSum = ddply(dF, .(gene), summarise,
                "gene negative" = length(GS[GS=="gene negative"]),
                "gene positive" = length(GS[GS=="gene positive"]),
                "affected"      = length(GS[GS=="affected"])
  )
  ptSum <- data.frame(ptSum, Totals=rowSums(ptSum[,-1]))
  tmpdF <- data.frame(gene="Totals",
                      as.data.frame(matrix(colSums(ptSum[,-1]), nrow=1))
  )
  names(tmpdF) <- names(ptSum)
  ptSum <- rbind(ptSum, tmpdF)
  
  print(xtable(ptSum,
               caption="Subjects included in the analysis",
               digits = c(0,0,0,0,0,0),
               display = c("s", "s", "d", "d", "d","d")),
        include.rownames=FALSE,
        # file=outFile,
        append=TRUE)
  
  
  
  # stack data and take the mean if it is a node-wise measures
  dF.wb <- stackIt(dF, metric)
  
  dF.wb.summary = ddply(dF.wb, .(), summarise,
                        "gene negative" = paste(sapply(mean(values[GS=="gene negative"], na.rm = TRUE), fn, a=2,b=3),
                                                paste("(", sapply(sd(values[GS=="gene negative"], na.rm = TRUE), fn, a=2,b=3),   ")", sep="")),
                        "gene positive" = paste(sapply(mean(values[GS=="gene positive"], na.rm = TRUE), fn, a=2,b=3),
                                                paste("(", sapply(sd(values[GS=="gene positive"], na.rm = TRUE),fn, a=2,b=3),   ")", sep="")),
                        "affected"      = paste(sapply(mean(values[GS=="affected"], na.rm = TRUE),fn, a=2,b=3),
                                                paste("(", sapply(sd(values[GS=="affected"], na.rm = TRUE),fn, a=2,b=3),   ")", sep=""))
                        )
  
  
  print(xtable(dF.wb.summary[,-1],
               caption=paste("Mean and standard deviations for",metricName,"values in individuals")),
        include.rownames=FALSE,
        file=outFile,
        append=TRUE)
    
  

  # ANOVA of the differences
  mod <- lm(values~GS*gene, data=dF.wb)
  mod.aov <- anova(mod)
  print(xtable(mod.aov,
               digits = c(0,0,2,2,1,2),
               display = c("s","d", "fg", "fg", "f", "fg")),
        include.rownames=TRUE,
        file=outFile,
        append=TRUE)
  
  
  print(xtable(mod.aov,
               digits = c(0,0,2,2,1,2),
               display = c("s","d", "fg", "fg", "f", "fg")),
        include.rownames=TRUE,
        append=TRUE)
  
  
  #### post hoc t-tests ####
  # group pairwise
  pwGS <- combn(levels(dF.wb$GS),2)
  pwTtests <- apply(pwGS, 2, function(x) dotTests(x,dF.wb))
  pwTtests <- data.frame(matrix(unlist(pwTtests), nrow = length(pwTtests), byrow = TRUE))
  names(pwTtests) <- c("comparison", "t", "p")
  pwTtests[,2] <- as.numeric(as.character(pwTtests[,2]))
  pwTtests[,3] <- as.numeric(as.character(pwTtests[,3]))

  print(xtable(pwTtests,
               caption="Pairwise post-hoc t-test results"),
        include.rownames=FALSE,
        file=outFile,
        append=TRUE)
  
  print(xtable(pwTtests,
               caption="Pairwise post-hoc t-test results"),
        include.rownames=FALSE,
        append=TRUE)
  
  
  # t-test for affected vs non-affected (gene negative and gene positive)
  anTtest <- t.test(dF.wb[dF.wb$GS=="affected","values"], dF.wb[dF.wb$GS!="affected","values"])
  
  print(xtable(data.frame("t"=as.numeric(as.character(fn(anTtest[["statistic"]]))),
                          "p"=as.numeric(as.character(fn(anTtest[["p.value"]])))),
               caption="t-test between affected subjects and non-affected (both gene negative and gene positive unaffected)",
               digits=c(0,2,2),
               display=c("s","fg","fg")
               ),
        include.rownames=FALSE,
        file=outFile,
        append=TRUE)
  
  
  
  print(xtable(data.frame("t"=fn(anTtest[["statistic"]]),
                          "p"=fn(anTtest[["p.value"]])),
               caption="t-test between affected subjects and non-affected (both gene negative and gene positive unaffected",
               digits=c(0,2,2),
               display=c("s","fg","fg")
               ),
        include.rownames=FALSE,
        append=TRUE)
        
  # now do a plot for group differences, collapsed across genes
  p <- ggplot(dF.wb, aes_string(x="GS", y="values", fill="GS"))
  p <- p + geom_boxplot()
  p <- p + scale_colour_manual(name="Group",values=as.vector(colList))
  
  xList = seq_along(levels(dF$GS))
  names(xList) = levels(dF$GS)
  for(n in seq(length(pwGS[1,]))){
    pVal = as.numeric(as.character(pwTtests[n, "p"]))
    d1 = pwGS[1,n]
    d2 = pwGS[2,n]
    p <- addSigBar(p, pVal, d1, d2, xList)
  }
  
  p <- p + theme_bw()
  p <- p + labs(y=metricName) + theme(axis.title.x=element_blank())
  p <- p + theme(text=element_text(size=ts))
  ggsave(paste(outDir,
               paste(metric,"allgroups.png",sep="_"),
               sep="/"))
  
  dF.wb.summary = ddply(dF.wb, .(gene), summarise,
                        "gene negative" = paste(sapply(mean(values[GS=="gene negative"], na.rm = TRUE), fn),
                                                paste("(", sapply(sd(values[GS=="gene negative"], na.rm = TRUE),fn),   ")", sep="")),
                        "gene positive" = paste(sapply(mean(values[GS=="gene positive"], na.rm = TRUE), fn),
                                                paste("(", sapply(sd(values[GS=="gene positive"], na.rm = TRUE),fn),   ")", sep="")),
                        "affected"      = paste(sapply(mean(values[GS=="affected"], na.rm = TRUE),fn),
                                                paste("(", sapply(sd(values[GS=="affected"], na.rm = TRUE), fn),   ")", sep=""))
  )
  
  print(xtable(dF.wb.summary,
               caption=paste("Mean and standard deviations for",metricName," values in individuals")),
        file=outFile,
        append=TRUE)
  
  print(xtable(dF.wb.summary,
               caption=paste("Mean and standard deviations for",metricName," values in individuals")),
        append=TRUE)
  
  # # plot this
  # err <- function(x) qnorm(0.975) * sd(x, na.rm = TRUE)/sqrt(length(x)) # function for standard error
  # 
  # dF.stacked.plot = ddply(dF.stacked, .(gene, GS), summarise, valuesMean=mean(values, na.rm = TRUE),
  #                         upper = valuesMean+err(values),
  #                         lower = valuesMean-err(values)
  #                         )
  
  # dF.stacked.plot <- data.frame(dF.stacked.plot, unique=seq(length(dF.stacked.plot[,1])))
  
  # Do t-tests between diagnostic groups within each gene
  dF.tResults <- data.frame()
  for(gene in levels(dF.wb$gene)){
    dF.temp <- dF.wb[dF.wb$gene==gene,]
    pwTtests <- apply(pwGS, 2, function(x) dotTests(x,dF.temp))
    pwTtests <- data.frame(gene=gene, matrix(unlist(pwTtests), nrow = length(pwTtests), byrow = TRUE))
    names(pwTtests) <- c("gene", "comparison", "t", "p")
    dF.tResults <- rbind(dF.tResults,
                         pwTtests)
  }
  
  dF.tResults[,3] <- as.numeric(as.character(dF.tResults[,3]))
  dF.tResults[,4] <- as.numeric(as.character(dF.tResults[,4]))
  
  print(xtable(dF.tResults,
               caption="t-test between diagnostic groups within each gene",
               digits=c(0,0,0,2,2),
               display=c("s","s","s","fg","fg")),
        include.rownames=FALSE,
        file=outFile,
        append=TRUE)
  
  print(xtable(dF.tResults,
               caption="t-test between diagnostic groups within each gene",
               digits=c(0,0,0,2,2),
               display=c("s","s","s","fg","fg")),
        include.rownames=FALSE,
        append=TRUE)
  
  # now do a plot for group differences for separate genes
  # p <- ggplot(dF.stacked.plot, aes_string(x="unique", y="valuesMean", middle="valuesMean", group="unique", fill="gene", upper="upper", lower="lower"))
  pg <- ggplot(dF.wb, aes_string(x="gene", y="values", fill="GS"))
  pg <- pg + geom_boxplot()
  pg <- pg + scale_colour_manual(name="Group",values=as.vector(colList))
  
  xList = seq_along(levels(dF$GS))
  names(xList) = levels(dF$GS)
  ymax = ggplot_build(pg)$panel$ranges[[1]]$y.range[2]
  ymin = ggplot_build(pg)$panel$ranges[[1]]$y.range[1]
  
  ytop.max = ymax
  for(ng in seq_along(levels(dF.wb$gene))){
    ytop = ymax
    gene = levels(dF.wb$gene)[ng]
    for(n in seq(length(pwGS[1,]))){
      pVal = as.numeric(as.character(dF.tResults[dF.tResults$gene==gene,][n, "p"]))
      d1 = pwGS[1,n]
      d2 = pwGS[2,n]
      pL <- addSigBarGenes(pg, pVal, d1, d2, xList, ytop, ng)
      pg = pL[["p"]]
      ytop = pL[["ymax"]]
    }
    if(ytop>ytop.max){ ytop.max=ytop }
  }
  pg <- pg + scale_y_continuous(limits=c(ymin, ytop.max))
  pg <- pg + theme_bw()
  pg <- pg + labs(y=metricName) + theme(axis.title.x=element_blank())
  pg <- pg + theme(text=element_text(size=ts))
  
  ggsave(paste(outDir,
               paste(metric,"byGene.png",sep="_"),
               sep="/"))
  endLog(outFile)
  
  return(list(p,pg))
}

wholeBrainAnalysisMixedEffects <- function(metric,
                                           metricName,
                                           cols,
                                           sp, weighted=TRUE,
                                           outDir="wholebrainResults",
                                           edgePC=3, ts=12,
                                           exclNeg=FALSE){
  
  if(weighted){
    metric = paste(metric,"wt",sep="_")
  }
  
  # create output directory
  dir.create(outDir, showWarnings = FALSE)
  
  # define log output file
  outFile = paste(outDir,
                  paste(metric,"logFile.tex",sep="_"),
                  sep="/")
  
  # import graph data
  dF <- importGraphData(metric, weighted, edgePC)
  
  # if indicated exclude gene negative group
  if(exclNeg){
    dF <- dF[dF$GS!="gene negative",]
  }
  
  # apply spike percentage threshold
  dF <- applySP(dF, sp)
  
  colList = unlist(cols[levels(dF$GS)])
  
  # stack data and take the mean if it is a node-wise measures
  dF.wb <- stackIt(dF, metric)

  # ANOVA of the differences
  mod <- lmer(values ~ GS + (1 | gene) + (1 | site) + (1 | Family),
              data=dF.wb,
              REML=FALSE)
  
  # print random effects of the model
  mod.coef <- ranef(mod)
  
  print(xtable(mod.coef$site,
               digits=c(0,2),
               display=c("s", "fg"),
               caption = "Linear mixed effects model, site coefficients"),
        file=outFile,
        append=TRUE)
  
  print(xtable(mod.coef$gene,
               digits=c(0,2),
               display=c("s", "fg"),
               caption = "Linear mixed effects model, gene coefficients"),
        file=outFile,
        append=TRUE)

  print(xtable(mod.coef$site,
                 digits=c(0,2),
                 display=c("s", "fg"),
               caption = "Linear mixed effects model, site coefficients"))
  
  
  print(xtable(mod.coef$gene,
               digits=c(0,2),
               display=c("s", "fg"),
               caption = "Linear mixed effects model, gene coefficients"))
  
  
  # print variances
  vc <- VarCorr(mod)

  print(xtable(data.frame(vc),
               display=c("s","s","s","s","g","fg"),
               digits=c(0,0,0,0,2,2),
               caption="Variance of random effects"),
        include.rownames=FALSE)

  print(xtable(data.frame(vc),
               display=c("s","s","s","s","g","fg"),
               digits=c(0,0,0,0,2,2),
               caption="Variance of random effects"),
        file=outFile,
        append=TRUE,
        include.rownames=FALSE)
  
  

  # Type II ANOVA
  mod.avo <- Anova(mod, type="II")
  
  print(xtable(mod.avo,
               digits=c(0,1,1,2),
               display=c("s","f","d","fg"),
               caption=paste("ANOVA of linear mixed effects model for", metricName))
        )
  
  print(xtable(mod.avo,
               digits=c(0,1,1,2),
               display=c("s","f","d","fg"),
               caption=paste("ANOVA of linear mixed effects model for", metricName)),
        file=outFile,
        append=TRUE
  )
  
  
  return(mod)
}

graphTimeComparison <- function(metric,
                                metricName,
                                sp, # spike percentage
                                cols,
                                weighted=TRUE, # is this a weighted metric?
                                outDir="wholebrainVsAOOResults",
                                edgePC=3,
                                h=30,w=45,s=4,ts=12,ps=4,
                                exclNeg=TRUE){
  if(weighted){
    metric = paste(metric,"wt",sep="_")
  }
  
  # create output directory
  dir.create(outDir, showWarnings = FALSE)
  
  # define log output file
  outFile = paste(outDir,
                  paste(metric,"logFile.tex",sep="_"),
                  sep="/")
  
  # create log file
  initiateLog(outFile, metricName)
  
  ### function to plot and analyse the relationship between graph metrics and expected time to disease onset
  # import graph metric data
  dF <- importGraphData(metric, weighted, edgePC)
  
  # filter by spike percentage
  dF <- applySP(dF, sp)
  
  # stack the data and take mean if a nodewise measure
  dF <- stackIt(dF, metric)
  
  # get age of onset data
  genfiData <- read.table("/home/tim/GENFI/genfi_Subjects_sjones_1_22_2015_17_47_47_restructure_summary.csv",
                          sep="\t",
                          header = TRUE)
  
  dF.aoo <- data.frame(wbic=genfiData$Subject,
                       aoo=genfiData$Yrs.from.AV_AAO)
  
  # merge age of onset data
  dF <- merge(dF, dF.aoo, by="wbic")
  
  # if indicated exclude gene negative group
  plotOutName = paste(metric,"png",sep=".")
  
  if(exclNeg){
    dF <- dF[dF$GS!="gene negative",]
    plotOutName = paste(paste(metric, "GenePos", sep=""),"png",sep=".")
  }

  lm.all <- lm(values ~ aoo*GS, data=dF)

  # plot linear model NOT MIXED EFFECTS
  dF.plot <- dF
  
  p <- ggplot(dF.plot, aes_string(x="aoo", y="values", colour="GS", group="GS"))
  p <- p + geom_point(size=ps)
  pp <- p + geom_smooth(method="lm")
  colList = unlist(cols[levels(dF.plot$GS)])
  pp <- pp + scale_colour_manual(name="Group",values=as.vector(colList))
  pp <- pp + labs(title=paste(metricName, "linear regression", sep="\n"), y=metricName, x="Estimated age of onset") + theme(axis.title.x=element_blank(), legend.key=element_rect(fill="white", colour="white"))
  pp <- pp + theme_bw() + theme(legend.key = element_rect(colour="#FFFFFF", fill = "#FFFFFF"))
  pp <- pp + theme(text=element_text(size=ts), plot.title=element_text(size=ts))
  
  plot(pp)
  
  ggsave(paste(outDir,
               plotOutName,
               sep="/"),
         scale=s,
         dpi=600,
         height=h, width=w,
         units="mm")
  print(xtable(summary(lm.all),
               caption=paste("summary of linear model for the interaction between estimated age of onset (aoo) and gene status (GS) on", metricName),
               digits=c(0,2,2,2,2),
               display=c("s","fg","fg","fg","g")),
        include.rownames=TRUE,
        file=outFile,
        append=TRUE)
  
  print(xtable(summary(lm.all),
               caption=paste("summary of linear model for the interaction between estimated age of onset (aoo) and gene status (GS) on", metricName),
               digits=c(0,2,2,2,2),
               display=c("s","fg","fg","fg","g")),
        include.rownames=TRUE,
        append=TRUE)

  # now run linear mixed effects model model
  # run model
  mod <- lmer(values ~ GS * aoo + (1 | gene) + (1 | site) + (1 | Family),
              data=dF,
              REML=FALSE)
  
  print(xtable(summary(mod)[["coefficients"]],
               caption="Linear mixed effects model, fixed effects",
               digits=c(0,2,2,2),
               display=c("s","fg","fg","fg")),
        include.rownames=TRUE,
        file=outFile,
        append=TRUE)
  
  print(xtable(summary(mod)[["coefficients"]],
               caption="Linear mixed effects model, fixed effects",
               digits=c(0,2,2,2),
               display=c("s","fg","fg","fg")),
        include.rownames=TRUE,
        append=TRUE)

  print(xtable(data.frame(StdDev = c(attributes(VarCorr(mod)[[1]])[["stddev"]],
                                     attributes(VarCorr(mod)[[2]])[["stddev"]],
                                     attributes(VarCorr(mod))[["sc"]]),
                          row.names = c("Family", "Site", "Residual")),                          
               caption="Linear mixed effects model, random effects",
               digits=c(0,2),
               display=c("s","fg")),
        include.rownames=TRUE,
        file=outFile,
        append=TRUE)
  
  # print the variance and standard deviation of the random effects
  vc <- VarCorr(mod)
  
  print(xtable(data.frame(vc),
               display=c("s","s","s","s","g","fg"),
               digits=c(0,0,0,0,2,2),
               caption="Variance of random effects"),
        include.rownames=FALSE)
  
  print(xtable(data.frame(vc),
               display=c("s","s","s","s","g","fg"),
               digits=c(0,0,0,0,2,2),
               caption="Variance of random effects"),
        file=outFile,
        append=TRUE,
        include.rownames=FALSE)
  
  
  # null model
  nulMod <- lmer(values ~ aoo + (1 | gene) + (1 | site) + (1 | Family),
                 data=dF,
                 REML=FALSE)

  modComparison <- anova(mod,nulMod)
  
  print(xtable(modComparison,
               caption = "ANOVA between model and null model to account for variance explained by the gene status",
               digits = c(0,0,2,2,2,2,2,2,2),
               display = c("s","d","fg","fg","fg","fg","d","fg","fg")),
        include.rownames=TRUE,
        file=outFile,
        append=TRUE)
  
  print(xtable(modComparison,
               caption = "ANOVA between model and null model to account for variance explained by the gene status",
               digits = c(0,0,2,2,2,2,2,2,2),
               display = c("s","d","fg","fg","fg","fg","d","fg","fg")),
        include.rownames=TRUE,
        append=TRUE)
  
  # plot linear model NOT MIXED EFFECTS
  dF.plot <- dF

  p <- ggplot(dF.plot, aes_string(x="aoo", y="values", colour="GS", group="GS"))
  p <- p + geom_point(size=ps)
  pp <- p + geom_smooth(method="lm")
  colList = unlist(cols[levels(dF.plot$GS)])
  pp <- pp + scale_colour_manual(name="Group",values=as.vector(colList))
  pp <- pp + labs(title=paste(metricName, "linear regression", sep="\n"), y=metricName, x="Estimated age of onset") + theme(axis.title.x=element_blank(), legend.key=element_rect(fill="white", colour="white"))
  pp <- pp + theme_bw() + theme(legend.key = element_rect(colour="#FFFFFF", fill = "#FFFFFF"))
  pp <- pp + theme(text=element_text(size=ts), plot.title=element_text(size=ts))

  ggsave(paste(outDir,
               plotOutName,
               sep="/"),
         scale=s,
         dpi=600,
         height=h, width=w,
         units="mm")
  
  return(mod)
}
  
graphTimeComparisonNL <- function(metric,
                                metricName,
                                sp, # spike percentage
                                cols,
                                startvec=NULL, # starting vectors for equation with quadratic term
                                startvecCub=c(dc=0.000000001), # starting vectors for equation with cubic term
                                weighted=TRUE, # is this a weighted metric?
                                outDir="wholebrainVsAOOResults",
                                edgePC=3,
                                h=30,w=45,s=4,ts=12,ps=4,
                                sink=TRUE){
  if(weighted){
    metric = paste(metric,"wt",sep="_")
  }
  
  # create output directory
  dir.create(outDir, showWarnings = FALSE)
  
  # define log output file
  outFile = paste(outDir,
                  paste(metric,"logFile.tex",sep="_"),
                  sep="/")
  
  # create log file
  initiateLog(outFile, metricName)
  
  ### function to plot and analyse the relationship between graph metrics and expected time to disease onset
  # import graph metric data
  dF <- importGraphData(metric, weighted, edgePC)
  
  # filter by spike percentage
  dF <- applySP(dF, sp)
  
  # stack the data and take mean if a nodewise measure
  dF <- stackIt(dF, metric)
  
  # get age of onset data
  genfiData <- read.table("/home/tim/GENFI/genfi_Subjects_sjones_1_22_2015_17_47_47_restructure_summary.csv",
                          sep="\t",
                          header = TRUE)
  
  dF.aoo <- data.frame(wbic=genfiData$Subject,
                       aoo=genfiData$Yrs.from.AV_AAO)
  
  # merge age of onset data
  dF <- merge(dF, dF.aoo, by="wbic")
  
  # now run non-linear mixed effects model with a quadratic term
  nlOutFile = paste(outDir,
                    paste(metric,"NLlogFile.txt",sep="_"),
                    sep="/")
  if(sink){
    sink(nlOutFile)
  }
  # firstly, plot the values to be able to estimate starting parameters
  if(is.null((startvec))){
    affirm = "n"
    while(affirm=="n"){
      print(pp)
      
#       Asym = as.numeric(readline(prompt="Asymptote of the y-axis: "))
#       R0   = as.numeric(readline(prompt="the value of y when x=0: "))
#       lrc  = as.numeric(readline(prompt="natural log of the rate of decline (guess -1 if not sure): "))
      aq = as.numeric(readline(prompt = "Constant term: "))
      bq = as.numeric(readline(prompt = "x multiplier:" ))
      cq = as.numeric(readline(prompt = "quadratic multiplier: "))
      
      # dF.temp <- data.frame(x=dF$aoo, y=SSasymp(dF$aoo, Asym, R0, lrc), GS="Start estimates")
      colList = unlist(cols[levels(dF.plot$GS)])
      p.temp <- p.temp + scale_colour_manual(name="Group",values=as.vector(colList))
      dF.temp <- data.frame(x=dF$aoo, y=SSquadFun(dF$aoo, aq=aq, bq=bq, cq=cq), GS="Start estimates")
      p.temp <- p.temp + geom_point(data=dF.temp, aes_string(x="x", y="y"))
      print(p.temp)
      affirm = readline(prompt = "Are you happy with this?(y/n)")
    }
    # startvec = c(Asym=Asym, R0=R0, lrc=lrc)
    startvec = c(aq=aq, bq=bq, cq=cq)
  }
  
  # save starting vectors for quadratic non-linear regression
  qStartvec = startvec
#   lapply(list("Starting values are as follows.",
#               paste("Constant:", startvec["aq"]),
#               paste("Multiplier of x term:", startvec["bq"]),
#               paste("Multiplier of quadratic term:", startvec["cq"])),
#          function(x) print(x)
#         )
  
  # run non-linear quadratic model
  # dF <- cbind(dF, data.frame(aooNeg=dF$aoo*-1))  # make age of onset negative so that an asymptotic fit can be achieved.

  # The following line calculates the non-linear model
  # the random elements take in to account the difference in graph measures depending on the scanner site (Asym|site) and the age at onset between genes (aooNeg|gene)
  # at present the function is an asymptote.

  nlmq <- nlmer(values ~ SSquadFun(aoo, aq, bq, cq) ~ (aoo|gene) + (aq|site), data=dF, start = startvec)
  
  print(summary(nlmq))

  # plot estimated values
  dF.plot <- dF
  p <- ggplot(dF.plot, aes_string(x="aoo", y="values", colour="GS", group="GS"))
  p <- p + geom_point(size=ps)
  
  x = seq(-45,25,by=1)
  dF.nlmq <- data.frame(x=x, y=SSquadFun(x, nlmq@beta[[1]], nlmq@beta[[2]], nlmq@beta[[3]]), GS="estimates")
  
  pq <- p + geom_line(data=dF.nlmq, aes_string(x="x", y="y"))
  cList = c("affected", "estimates", "gene positive")
  colList = unlist(cols[cList])
  pq <- pq + scale_colour_manual(name="Group",values=as.vector(colList),
                                 breaks=c("gene positive", "affected", "estimates"))
  pq <- pq + theme_bw() + theme(legend.key=element_rect(fill="white", colour="white"))
  pq <- pq + labs(title=paste(metricName, "non-linear regression, quadratic", sep="\n"), y=metricName, x="Estimated time from onset")# + theme(axis.title.x=element_blank())
  pq <- pq + theme(text=element_text(size=ts), plot.title=element_text(size=ts), plot.title=element_text(size=ts))
  
#   print(paste("Saving",paste(outDir,
#                        paste(paste(metric, "nonLinearEstimatesQuadratic",sep="_"),"png",sep="."),
#                        sep="/")))
  ggsave(paste(outDir,
               paste(paste(metric, "nonLinearEstimatesQuadratic",sep="_"),"png",sep="."),
               sep="/"),
         scale=s,
         dpi=600,
         height=h, width=w,
         units="mm")
  
  # now run non-linear mixed effects model with a cubic term
  # firstly, plot the values to be able to estimate starting parameters
  ac = nlmq@beta[[1]]
  bc = nlmq@beta[[2]]
  cc = nlmq@beta[[3]]
  
  if(is.null((startvecCub))){
    affirm = "n"
    while(affirm=="n"){
      print(pp)
      dc = as.numeric(readline(prompt = "cubic multiplier: "))
      
      # dF.temp <- data.frame(x=dF$aoo, y=SSasymp(dF$aoo, Asym, R0, lrc), GS="Start estimates")
      dF.temp <- data.frame(x=dF$aoo, y=SScubicFun(dF$aoo, ac=ac, bc=bc, cc=cc, dc=dc), GS="Start estimates")
      p.temp <- pp + geom_point(data=dF.temp, aes_string(x="x", y="y"))
      print(p.temp)
      affirm = readline(prompt = "Are you happy with this?(y/n)")
    }
    # startvec = c(Asym=Asym, R0=R0, lrc=lrc)
  } else {
    dc=startvecCub[["dc"]]
    print(dc)
  }
  
  startvec = c(ac=ac, bc=bc, cc=cc, dc=dc)
  cStartvec = startvec # save starting values for the cubic non-linear regression model
  
#   lapply(list("Starting values are as follows.",
#               paste("Constant:", startvec["ac"]),
#               paste("Multiplier of x term:", startvec["bc"]),
#               paste("Multiplier of quadratic term:", startvec["cc"]),
#               paste("Multiplier of cubic term:", startvec["dc"])),
#          function(x) print(x)
#   )
  
  # run non-linear quadratic model
  # The following line calculates the non-linear model
  # the random elements take in to account the difference in graph measures depending on the scanner site (aq|site) and the age at onset between genes (aooNeg|gene)
  # this equation includes a cubic term
  nlmc <- nlmer(values ~ SScubicFun(aoo, ac, bc, cc, dc) ~ (aoo|gene) + (ac|site), data=dF, start = startvec)

  # print a summary of the model
  print(summary(nlmc))
  
  # plot estimated values
  x = seq(-45,25,by=1)
  dF.nlmc <- data.frame(x=x,
                        y=SScubicFun(x, nlmc@beta[[1]], nlmc@beta[[2]], nlmc@beta[[3]], nlmc@beta[[4]]),
                        GS="estimates")
  
  # pc <- ggplot()
  pc <- p + geom_line(data=dF.nlmc, aes_string(x="x", y="y"))
  pc <- pc + scale_colour_manual(name="Group",values=as.vector(colList),
                                 breaks=c("gene positive", "affected", "estimates"))
  pc <- pc + theme_bw() + theme(legend.key=element_rect(fill="white", colour="white"))
  pc <- pc + labs(title=paste(metricName, "non-linear regression, cubic", sep="\n"), y=metricName, x="Estimated time from onset")# + theme(axis.title.x=element_blank())
  pc <- pc + theme(text=element_text(size=ts), plot.title=element_text(size=ts))
  
#   print(paste("Saving",paste(outDir,
#                              paste(paste(metric, "nonLinearEstimatesCubic",sep="_"),"png",sep="."),
#                              sep="/")))
  ggsave(paste(outDir,
               paste(paste(metric, "nonLinearEstimatesCubic",sep="_"),"png",sep="."),
               sep="/"),
         scale=s,
         dpi=600,
         height=h, width=w,
         units="mm")
  
  if(sink){
    sink()
  }
  
  # now compare the quadratic and cubic equations
  qvsc.anova = NULL
  try(
    qvsc.anova <- anova(nlmq, nlmc)
  )
  
  if(!is.na(qvsc.anova)){
    print(xtable(qvsc.anova,
                 caption = "ANOVA between models with quadratic and cubic terms to explain rate of change in graph measure with time",
                 digits = c(0,0,2,2,2,2,2,2,2),
                 display = c("s","d","fg","fg","fg","fg","d","fg","fg")),
          include.rownames=TRUE,
          file=outFile,
          append=TRUE)
  }
  
  # finalise log file and return dataframe
  endLog(outFile)
  return(list(list(pp=pp, pq=pq, pc=pc),
              qStartvec,     # starting values for quadratic non-linear regression
              cStartvec,     # starting values for cubic non-linear regression
              summary(nlmq), # summary of non-linear quatric non-linear regression
              summary(nlmc), # summary of non-linear cubic non-linear regression
              qvsc.anova     # ANOVA between models with quadratic and cubic terms
              )
  )
}


# import spike percentage data
# sp <- read.xlsx("../all_sm_thld10_SP.xlsx", sheetIndex = 1)

# metrics = list("degree"="connection strength",
#                "ge"="global efficiency",
#                "le"="local efficiency",
#                "pl"="path length",
#                "eigCentNP"="eigen centrality",
#                "betCent"="betweenness centrality",
#                "closeCent"="closeness centrality")

# lapply(names(metrics), function(x) wholebrainAnalysis(x, metrics[[x]], sp))
# lapply(names(metrics), function(x) graphTimeComparison(x, metrics[[x]], sp))

# do metrics for unweighted graphs
# metrics = list("degree"="connection strength",
#                "degreeNorm"="connection strength (normalised)",
#                "ge"="global efficiency",
#                "geNorm"="global efficiency (normalised)",
#                "le"="local efficiency",
#                "leNorm"="local efficiency (normalised)",
#                "pl"="path length",
#                "plNorm"="path length (normalised)",
#                "eigCentNP"="eigen centrality",
#                "eigCentNorm"="eigen centrality (normalised)",
#                "betCent"="betweenness centrality",
#                "betCentNorm"="betweenness centrality (normalised)",
#                "closeCent"="closeness centrality",
#                "closeCentNorm"="closeness centrality (normalised)")

# lapply(names(metrics), function(x) wholebrainAnalysis(x, metrics[[x]], sp, weighted=FALSE))


startVecListQ = list(degree        = c(aq = 29.8,   bq = 0.000001, cq = 0.000001),
                     degreeNorm    = c(aq = 1.03,   bq = 0.000001, cq = 0.000001),
                     ge            = c(aq = 0.44,   bq = 0.000001, cq = 0.000001),
                     geNorm        = c(aq = 0.89,   bq = 0.000001, cq = 0.000001),
                     le            = c(aq = 20.5,   bq = 0.000001, cq = 0.000001),
                     leNorm        = c(aq = 1.30,   bq = 0.000001, cq = 0.000001),
                     pl            = c(aq = 2.50,   bq = 0.000001, cq = 0.000001),
                     plNorm        = c(aq = 1.18,   bq = 0.000001, cq = 0.000001),
                     eigCentNP     = c(aq = 0.037,  bq = 0.000001, cq = 0.000001),
                     eigCentNorm   = c(aq = 0.85,   bq = 0.000001, cq = 0.000001),
                     betCent       = c(aq = 0.0031, bq = 0.000001, cq = 0.000001),
                     betCentNorm   = c(aq = 1.30,   bq = 0.000001, cq = 0.000001),
                     closeCent     = c(aq = 0.40,   bq = 0.000001, cq = 0.000001),
                     closeCentNorm = c(aq = 0.85,   bq = 0.000001, cq = 0.000001)
                     )

cols <- list("gene negative"="#E69F00",
             "gene positive"="#56B4E9",
             "affected"="#D55E00",
             "estimates"="#000000"
             )

# startVecListC= list(geNorm = c(dc=0.0001),
#                     plNorm = c(dc=0.001))
# 
# dF <- graphTimeComparison("plNorm", "path length (normalised)",
#                           startvec=startVecListQ[["plNorm"]],
#                           startvecCub = startVecListC[["plNorm"]],
#                           sp,
#                           cols=cols,
#                           weighted=FALSE,
#                           sink=FALSE)

# lapply(names(metrics), function(x) graphTimeComparison(x,
#                                                        metrics[[x]],
#                                                        sp, weighted=FALSE,
#                                                        cols=cols,
#                                                        startvec = startVecListQ[[x]]))

breakpoint <- function(metric,
                       metricName,
                       sp, # spike percentage
                       cols,
                       edgePC=3,
                       outDir="wholebrainVsAOOResults",
                       weighted=FALSE,
                       ts=12){
  if(weighted){
    metric = paste(metric,"wt",sep="_")
  }
  
  # create output directory
  dir.create(outDir, showWarnings = FALSE)
  
  # define log output file
  outFile = paste(outDir,
                  paste(metric,"logFile.tex",sep="_"),
                  sep="/")
  
  # create log file
  initiateLog(outFile, metricName)
  
  ### function to plot and analyse the relationship between graph metrics and expected time to disease onset
  # import graph metric data
  dF <- importGraphData(metric, weighted, edgePC)
  
  # filter out gene negative subjects
  dF <- dF[dF$GS!="gene negative",]
  
  # filter by spike percentage
  dF <- applySP(dF, sp)
  
  # stack the data and take mean if a nodewise measure
  dF <- stackIt(dF, metric) # changes the name of the metric to 'values'
  
  # get age of onset data
  genfiData <- read.table("/home/tim/GENFI/genfi_Subjects_sjones_1_22_2015_17_47_47_restructure_summary.csv",
                          sep="\t",
                          header = TRUE)
  
  dF.aoo <- data.frame(wbic=genfiData$Subject,
                       aoo=genfiData$Yrs.from.AV_AAO)
  
  # merge age of onset data
  dF <- merge(dF, dF.aoo, by="wbic")

#   # filter out affected subjects less than t=0 and gene positive greater than t=0
#   dF <- rbind(dF[dF$GS=="affected" && dF$aoo>=0,],
#               dF[dF$GS=="gene positive" && dF$aoo<=0,])
  
  # create model
  mod <- lm(values ~ aoo, data=dF)

  # do segmentation analysis assuming the breakpoint is at 0
  br <- segmented(mod, seg.Z = ~aoo, psi = c(aoo=0.))

  # get summary
  br.summary <- summary(br)
  
  # do statistical testing
  dt = davies.test(mod, seg.Z=~aoo)
  print(dt)
  pVal = dt$p.value
  
  # plot breakpoint data
  y = br$psi[[2]] * br$coefficients[[2]] + br$coefficients[[1]] # get the y value of the breakpoint
  p <- ggplot(dF, aes_string(x="aoo", y="values", colour="GS"))
  p <- p + geom_point()
  p <- p + geom_segment(x=min(dF$aoo),
                        xend=br$psi[[2]],
                        y=min(dF$aoo)*br$coefficients[[2]]+br$coefficients[[1]],
                        yend=y,
                        colour="black")
  
  p <- p + geom_segment(x=br$psi[[2]],
                        xend=max(dF$aoo),
                        y=y,
                        yend=max(dF$aoo)*br$coefficients[[3]]+br$coefficients[[1]],
                        colour="black")
  
  colList = unlist(cols[levels(dF$GS)])
  p <- p + scale_colour_manual(name="Group",values=as.vector(colList))
  p <- p + labs(title=paste(metricName, "breakpoint analysis", sep="\n"), y=metricName, x="Estimated age of onset") + theme(axis.title.x=element_blank(), legend.key=element_rect(fill="white", colour="white"))
  p <- p + theme_bw() + theme(legend.key = element_rect(colour="#FFFFFF", fill = "#FFFFFF"))
  p <- p + theme(text=element_text(size=ts), plot.title=element_text(size=ts))
  
  return(list(br, p, pVal))
}


breakPointDiscontinuous <- function(metric,
                                    metricName,
                                    sp, # spike percentage
                                    cols,
                                    edgePC=3,
                                    outDir="wholebrainVsAOOResults",
                                    weighted=FALSE,
                                    ts=12){
  if(weighted){
    metric = paste(metric,"wt",sep="_")
  }
  
  # create output directory
  dir.create(outDir, showWarnings = FALSE)
  
  # define log output file
  outFile = paste(outDir,
                  paste(metric,"logFile.tex",sep="_"),
                  sep="/")
  
  # create log file
  initiateLog(outFile, metricName)
  
  ### function to plot and analyse the relationship between graph metrics and expected time to disease onset
  # import graph metric data
  dF <- importGraphData(metric, weighted, edgePC)
  
  # filter out gene negative subjects
  dF <- dF[dF$GS!="gene negative",]
  
  # filter by spike percentage
  dF <- applySP(dF, sp)
  
  # stack the data and take mean if a nodewise measure
  dF <- stackIt(dF, metric) # changes the name of the metric to 'values'
  
  # get age of onset data
  genfiData <- read.table("/home/tim/GENFI/genfi_Subjects_sjones_1_22_2015_17_47_47_restructure_summary.csv",
                          sep="\t",
                          header = TRUE)
  
  dF.aoo <- data.frame(wbic=genfiData$Subject,
                       aoo=genfiData$Yrs.from.AV_AAO)
  
  # merge age of onset data
  dF <- merge(dF, dF.aoo, by="wbic")
  
  #   # filter out affected subjects less than t=0 and gene positive greater than t=0
  #   dF <- rbind(dF[dF$GS=="affected" && dF$aoo>=0,],
  #               dF[dF$GS=="gene positive" && dF$aoo<=0,])
  
  # define breakpoint
  bkpt = 0.
  
#   piecewise <- lmer(values ~ aoo*(aoo<bkpt) + aoo*(aoo>bkpt) + (1 | gene) + (1 | site) + (1 | Family),
#                     data=dF,
#                     REML=FALSE)
  piecewise <- lm (values ~ aoo*(aoo<bkpt) + aoo*(aoo>bkpt), data=dF)

  print(summary(piecewise))
  
  # compare with null model
#   nulMod <- lmer(values ~ aoo + (1 | gene) + (1 | site) + (1 | Family),
#                   data=dF,
#                   REML=FALSE)
  nulMod <- lm(values ~ aoo,
               data=dF)
  
  modComparison <- anova(piecewise,nulMod)
  
  print(modComparison)
  
  print(xtable(modComparison,
               caption = "ANOVA between model and null model to assess whether the segmented model fits better",
               digits = c(0,0,2,2,2,2,2),
               display = c("s","d","fg","fg","fg","fg","fg")),
        include.rownames=TRUE,
        file=outFile,
        append=TRUE)
  
  print(xtable(modComparison,
               caption = "ANOVA between model and null model to assess whether the segmented model fits better",
               digits = c(0,0,2,2,2,2,2),
               display = c("s","d","fg","fg","fg","fg","fg")),
        include.rownames=TRUE,
        append=TRUE)
  
  p <- ggplot(dF, aes_string(x="aoo", y="values", colour="GS"))
  p <- p + geom_point()
  
  ### need to finish plot ###
  
  return(piecewise)
  
}
