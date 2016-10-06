# This script will examine whole brain graph metrics in the GENFI dataset by looking for differences
# between affected/carriers/non-carriers, and by examining the relationship between average age at
# onset and the graph measures.

library(plyr)
library(ggplot2)
library(xtable)
library(car)
library(xlsx)
library(lmerTest)
library(segmented)
library(MASS)

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

importGraphData <- function(metric, weighted, edgePC=3, Age=TRUE, lobe=NA, hubT=NA, nonHubs=FALSE){
  # define input file
  if(weighted){
    inFile = paste("d2",metric,"wt","local",sep="_")
  } else {
    inFile = paste("d2",metric,"local",sep="_")
  }
  
  if(!is.na(lobe)){
    inFile = paste(inFile,lobe,sep="_")
  }
  
  # import data
  dF = read.table(paste(genfiDir,inFile, sep="/"), header = TRUE, na.strings = "NA")
  
  # select only the desired percentage edge threshold if unweighted metric
  if(!weighted){
    dF <- dF[dF$edgePC==edgePC,]
  }
  
  # convert diagnostic label to a factor and rename
  dF$GS = as.factor(dF$GS)
  dF$GS = revalue(dF$GS, c("0" = "gene negative", "1"="gene carriers", "2"="FTD"))
  dF$site = as.factor(dF$site)
  
  # if there is a hub threshold, apply it now
  # if a hub threshold is defined, then firstly read the file that contains the hub definitions
  if(!is.na(hubT)){
    lines=readLines("/home/tim/GENFI/GENFI_camgrid_20150525/Control/degree_allSubjects.txt") # read the file as lines
    lineList = lapply(lines, function(x) strsplit(x, " ")) # extract the lines in to a list
    hubs = unlist(lineList[[which(sapply(lineList, function(x) return(x[[1]][[1]]))==as.character(hubT))]])[-1]  # extract the line starting with the required hub threshold (-1 removes the threshold from the final list)
    hubs = sapply(hubs, function(x) paste("X",x,sep = ""))
    
    if(nonHubs){
      dF <- dF[!sapply(names(dF), function(x) x %in% hubs)]
    } else {
      # now filter by the hubs
      X0Col = which(names(dF)=="X0")-1
      dF <- cbind(dF[1:X0Col], dF[hubs])
    }
  }
  
  if(Age){
    # add age
    dF.demog <- read.table("../genfi_Subjects_sjones_1_22_2015_17_47_47_restructure_summary.csv", sep="\t", header = TRUE)
    names(dF.demog)[4] <- "wbic"
    
    #   dF.demog <- dF.demog[dF.demog$wbic[dF$wbic %in% dF.demog$wbic],]
    dF <- merge(dF, dF.demog[,c(4,which(names(dF.demog)=="Age"))],
                by = c("wbic"),
                all.x=TRUE,
                all.y=FALSE)
    dF <- unique(dF)
  }
  
  # return dataframe with graph metrics
  return(dF)
}

applySP <- function(dF, sp, spVal=10){
  # this function filters by spike percentage, but also removes subjects with <100 volumes
  # combine data with spike percentage data
  sp.sub <- data.frame(wbic=sp$id, spMean=sp$mean)
  dF <- merge(dF, sp.sub, by="wbic")
  
  # filter by spike percentage
  dF <- dF[dF$spMean < spVal,]
  
  # filter by length of scan
  dF.demog <- read.table("../genfi_Subjects_sjones_1_22_2015_17_47_47_restructure_summary.csv", sep="\t", header = TRUE)
  names(dF.demog)[4] <- "wbic"
  dF.demog <- dF.demog[dF.demog$scan_duration>100,] # exclude lines with subjects with only short scans
  
  dF <- dF[dF$wbic %in% dF.demog[,"wbic"],]
  
  # return filtered data
  return(dF)
}

stackIt <- function(dF, metric){
  # Sort out stacking the nodes up if nodewise measure
  nodeNames = names(dF)[sapply(names(dF), function(x) grepl("X",x))]
  
  if(length(nodeNames)>0){
    # Take the mean values of
    dF.stacked = stack(dF[,nodeNames])
    dF.stacked = data.frame(dF.stacked, GS=dF$GS, gene=dF$gene, wbic=dF$wbic, site=dF$site, Family=dF$Family, Age=dF$Age)
    
    dF.wb = ddply(dF.stacked, .(gene, wbic, GS, site, Family, Age), summarise,
                  values = mean(values, na.rm = TRUE)
    )
  } else {
    dF.wb <- dF
    names(dF.wb)[names(dF.wb)==metric] <- "values"
  }
  
  return(dF.wb)
}

initiateLog <- function(logFile, metricName){
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
        file=logFile,
        append=FALSE)
}

endLog <- function(logFile){
  write("\\end{document}", file=logFile, append=TRUE)
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

addSigBar <- function(p, pVal, d1, d2, xList, sig.value=0.05, excludeNegs=FALSE){
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

addSigBarGenes <- function(p, pVal, d1, d2, xList, ymax, ng,
                           nLen=3, sig.value=0.05, excludeNegs=FALSE){
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
    xpos1 = ng + (xList[d1]-(nLen/2)-0.5)*(1/(nLen+1))   #
    xpos2 = ng + (xList[d2]-(nLen/2)-0.5)*(1/(nLen+1)) # - 1/nLen - (xList[d2]-nLen+1)*0.7   #
    
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
                               outDir="wholeBrainResults",
                               edgePC=3,
                               h=15,w=15,s=4,tsz=12,ps=4,  # tsz = text size
                               excludeNegs=FALSE, # TRUE to exclude gene negative subjects
                               lobe=NA, # define the lobe of the brain to examine
                               hubT=NA # hub threshold
                               ){
  
  # create output directory
  dir.create(outDir, showWarnings = FALSE)
  
  # define log output file
  if(weighted){
    outFile = paste(outDir,
                    paste(metric,"wt",sep="_"),
                    sep="/")
    
  } else {
    outFile = paste(outDir,
                    paste(metric,sep="_"),
                    sep="/")
    
  }
  
  logFile = paste(outFile, "logFile.tex", sep="_")
  
  # create log file
  initiateLog(logFile, metricName)
  
  # import graph data
  dF <- importGraphData(metric, weighted, edgePC, lobe=lobe, hubT=hubT)
  
  # apply spike percentage threshold
  dF <- applySP(dF, sp)
  
  if(excludeNegs){
    dF <- dF[dF$GS!="gene negative",]
    dF$GS <- factor(as.character(dF$GS), levels=c("gene carriers", "FTD"))
  }
  
  colList = unlist(cols[levels(dF$GS)])
  
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
  
  t1 <- xtable(ptSum,
               caption="Subjects included in the analysis",
               digits = c(0,0,0,0,0,0),
               display = c("s", "s", "d", "d", "d","d"))
  
  # stack data and take the mean if it is a node-wise measures
  dF.wb <- stackIt(dF, metric)
  
  # add lobe information if necessary
  if(!is.na(lobe)){
    dF.wb <- data.frame(dF.wb, lobe=lobe)
  }
  
  if(!excludeNegs){
    dF.wb.summary = ddply(dF.wb, .(), summarise,
                          "gene negative" = paste(sapply(mean(values[GS=="gene negative"], na.rm = TRUE), fn, a=2,b=3),
                                                  paste("(", sapply(sd(values[GS=="gene negative"], na.rm = TRUE), fn, a=2,b=3),   ")", sep="")),
                          "gene carriers" = paste(sapply(mean(values[GS=="gene carriers"], na.rm = TRUE), fn, a=2,b=3),
                                                  paste("(", sapply(sd(values[GS=="gene carriers"], na.rm = TRUE),fn, a=2,b=3),   ")", sep="")),
                          "FTD"      = paste(sapply(mean(values[GS=="FTD"], na.rm = TRUE),fn, a=2,b=3),
                                                  paste("(", sapply(sd(values[GS=="FTD"], na.rm = TRUE),fn, a=2,b=3),   ")", sep=""))
    )
    
  } else {
    dF.wb.summary = ddply(dF.wb, .(), summarise,
                          "gene carriers" = paste(sapply(mean(values[GS=="gene carriers"], na.rm = TRUE), fn, a=2,b=3),
                                                  paste("(", sapply(sd(values[GS=="gene carriers"], na.rm = TRUE),fn, a=2,b=3),   ")", sep="")),
                          "FTD"      = paste(sapply(mean(values[GS=="FTD"], na.rm = TRUE),fn, a=2,b=3),
                                                  paste("(", sapply(sd(values[GS=="FTD"], na.rm = TRUE),fn, a=2,b=3),   ")", sep=""))
    )
    
  }
  
#   print(xtable(dF.wb.summary[,-1],
#                caption=paste("Mean and standard deviations for",metricName,"values in individuals")),
#         include.rownames=FALSE,
#         file=logFile,
#         append=TRUE)
  
  if(excludeNegs){
    dF.wb <- dF.wb[dF.wb$GS!="gene negative",]
    dF.wb$GS <- as.factor(as.character(dF.wb$GS))
  }
  
  # ANOVA of the differences
  mod <- lm(values~GS*gene, data=dF.wb)
  mod.aov <- anova(mod)

  t2 <- xtable(mod.aov,
               digits = c(0,0,2,2,1,2),
               display = c("s","d", "fg", "fg", "f", "fg"),
               caption = paste("ANOVA of the difference in", metricName, "between diagnostic groups"))
  
  print(t2,
        include.rownames=TRUE,
        file=logFile,
        append=TRUE)
  
  
  #### post hoc t-tests ####
  # group pairwise
  pwGS <- combn(levels(dF.wb$GS),2)
  pwTtests <- apply(pwGS, 2, function(x) dotTests(x,dF.wb))
  pwTtests <- data.frame(matrix(unlist(pwTtests), nrow = length(pwTtests), byrow = TRUE))
  names(pwTtests) <- c("comparison", "t", "p")
  pwTtests[,2] <- as.numeric(as.character(pwTtests[,2]))
  pwTtests[,3] <- as.numeric(as.character(pwTtests[,3]))
  
  t3 <- xtable(pwTtests,
               caption="Pairwise post-hoc t-test results",
               digits = c(0,0,2,2),
               display = c("s", "s", "fg", "fg"))
  
  print(t3,
        include.rownames=FALSE,
        file=logFile,
        append=TRUE)
  
  if(!excludeNegs){
    # t-test for affected vs non-affected (gene negative and gene positive)
    anTtest <- t.test(dF.wb[dF.wb$GS=="FTD","values"], dF.wb[dF.wb$GS!="FTD","values"])
    
    t4 <- xtable(data.frame("t"=fn(anTtest[["statistic"]]),
                            "p"=fn(anTtest[["p.value"]])),
                 caption="t-test between affected subjects and non-affected (both gene negative and gene positive unaffected",
                 digits=c(0,2,2),
                 display=c("s","fg","fg")
    )
    
    print(t4,
          include.rownames=FALSE,
          file=logFile,
          append=TRUE)
  } else {
    t4 = NA
  }
  
  if(excludeNegs){
    dF.wb.summary = ddply(dF.wb, .(gene), summarise,
                          "gene carriers" = paste(sapply(mean(values[GS=="gene carriers"], na.rm = TRUE), fn),
                                                  paste("(", sapply(sd(values[GS=="gene carriers"], na.rm = TRUE),fn),   ")", sep="")),
                          "FTD"      = paste(sapply(mean(values[GS=="FTD"], na.rm = TRUE),fn),
                                                  paste("(", sapply(sd(values[GS=="FTD"], na.rm = TRUE), fn),   ")", sep=""))
    )
    
  } else {
    dF.wb.summary = ddply(dF.wb, .(gene), summarise,
                          "gene negative" = paste(sapply(mean(values[GS=="gene negative"], na.rm = TRUE), fn),
                                                  paste("(", sapply(sd(values[GS=="gene negative"], na.rm = TRUE),fn),   ")", sep="")),
                          "gene carriers" = paste(sapply(mean(values[GS=="gene carriers"], na.rm = TRUE), fn),
                                                  paste("(", sapply(sd(values[GS=="gene carriers"], na.rm = TRUE),fn),   ")", sep="")),
                          "FTD"      = paste(sapply(mean(values[GS=="FTD"], na.rm = TRUE),fn),
                                                  paste("(", sapply(sd(values[GS=="FTD"], na.rm = TRUE), fn),   ")", sep=""))
    )
    
  }
  

  t5 <- xtable(dF.wb.summary,
               caption=paste("Mean and standard deviations for",metricName," values in individuals"))
  
  print(t5,
        file=logFile,
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
    pwTtestsG <- apply(pwGS, 2, function(x) dotTests(x,dF.temp))
    pwTtestsG <- data.frame(gene=gene, matrix(unlist(pwTtestsG), nrow = length(pwTtestsG), byrow = TRUE))
    names(pwTtestsG) <- c("gene", "comparison", "t", "p")
    dF.tResults <- rbind(dF.tResults,
                         pwTtestsG)
  }
  
  dF.tResults[,3] <- as.numeric(as.character(dF.tResults[,3]))
  dF.tResults[,4] <- as.numeric(as.character(dF.tResults[,4]))
  
  t6 <- xtable(dF.tResults,
               caption="t-test between diagnostic groups within each gene",
               digits=c(0,0,0,2,2),
               display=c("s","s","s","fg","fg"))
  
  print(t6,
        include.rownames=FALSE,
        file=logFile,
        append=TRUE)
  
  # now do a plot for group differences, collapsed across genes
  # reorder factors if gene negatives removed
  if(excludeNegs){
    dF.wb$GS <- factor(dF$GS, levels=c("gene carriers", "FTD"))
  }
  p <- ggplot(dF.wb, aes_string(x="GS", y="values", fill="GS"))
  p <- p + geom_boxplot()
  
  # add significance bar
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
  p <- p + theme(text=element_text(size=tsz))
  p <- p + scale_fill_manual(name="Group",values=as.vector(colList))
  p <- p + theme(legend.position="none")
  if(!is.na(lobe)){
    p <- p + labs(title=lobe)
  }
  
  if(excludeNegs){
    outw = 2.5
  } else {
    outw = 3
  }
  
  plotFile = paste(paste(outFile,"allgroups",sep="_"), "png", sep=".")
  ggsave(plotFile,
         scale=s,
         dpi=600,
         height=h, width=w,
         units="mm")
  
  # now do a plot for group differences for separate genes
  # p <- ggplot(dF.stacked.plot, aes_string(x="unique", y="valuesMean", middle="valuesMean", group="unique", fill="gene", upper="upper", lower="lower"))
  pg <- ggplot(dF.wb, aes_string(x="gene", y="values", fill="GS"))
  pg <- pg + geom_boxplot()
  pg <- pg + scale_fill_manual(name="Group",values=as.vector(colList))
  pg <- pg + theme(legend.position="none") 
  
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
      pL <- addSigBarGenes(pg, pVal, d1, d2, xList, ytop, ng, nLen=length(levels(dF.wb$GS)), excludeNegs=excludeNegs)
      pg = pL[["p"]]
      ytop = pL[["ymax"]]
    }
    if(ytop>ytop.max){ ytop.max=ytop }
  }
  pg <- pg + scale_y_continuous(limits=c(ymin, ytop.max))
  pg <- pg + theme_bw()
  pg <- pg + theme(axis.title.x=element_blank(), axis.title.y=element_blank())  #
  pg <- pg + theme(text=element_text(size=tsz))
  
  plotFile = paste(paste(outFile,"byGene",sep="_"), "png", sep=".")
  ggsave(plotFile,
         scale=s,
         dpi=600,
         height=h, width=w*2,
         units="mm")
  
  endLog(logFile)
  
  if(is.na(lobe)){
    return(list(t1, t2, t3, t4, t5, t6, p, pg))
  } else {
    return(list(t1, t2, t3, t4, t5, t6, p, pg, dF.wb))
  }
  
}

wholeBrainAnalysisMixedEffects <- function(metric,
                                           metricName,
                                           cols,
                                           sp, weighted=TRUE,
                                           outDir="wholeBrainResults",
                                           edgePC=3, tsz=12,
                                           exclNeg=FALSE,
                                           family=FALSE,
                                           lobe=NA, # define the lobe of the brain to examine
                                           hubT=NA){
  
  # create output directory
  dir.create(outDir, showWarnings = FALSE)
  
  # define log output file
  if(weighted){
    outFile = paste(outDir,
                    paste(metric,"wt",sep="_"),
                    sep="/")
    
  } else {
    outFile = paste(outDir,
                    paste(metric,sep="_"),
                    sep="/")
  }
  logFile = paste(outFile, "logFile.tex", sep="_")
  
  # create label for legends according to lobe
  if(!is.na(lobe)){
    lobeTag=lobe
  } else {
    lobeTag=""
  }
  
  # import graph data
  dF <- importGraphData(metric, weighted, edgePC, lobe=lobe, hubT=hubT)
  
  # apply spike percentage threshold
  dF <- applySP(dF, sp)
  
  colList = unlist(cols[levels(dF$GS)])
  
  # stack data and take the mean if it is a node-wise measures
  dF.wb <- stackIt(dF, metric)
  
  # add lobe information if necessary
  if(!is.na(lobe)){
    dF.wb <- data.frame(dF.wb, lobe=lobe)
  }
  
  # Set contrasts if control group are to be excluded
  # NB this retains the control group in the estimation of the age effect and as a null regressor
  if(exclNeg){
    cMat <- matrix(c(0,-1,1), nrow = 3)
    rownames(cMat) <- levels(dF.wb$GS)
    colnames(cMat) <- c("carriersVsFTD")
    attr(dF.wb$GS, "contrasts") = cMat
  }

  # ANOVA of the differences
  if(family){
    mod <- lmer(values ~ GS + Age + (1 | gene) + (1 | site) + (1 | Family),
                data=dF.wb,
                REML=FALSE)
  } else {
    mod <- lmer(values ~ GS + Age + (1 | gene) + (1 | site),
                data=dF.wb,
                REML=FALSE)
  }
  
  # print random effects of the model
  mod.coef <- ranef(mod)
  
  t1 <- xtable(mod.coef$site,
                 digits=c(0,2),
                 display=c("s", "fg"),
               caption = paste("Linear mixed effects model, site coefficients",lobeTag))

  t2 <- xtable(mod.coef$gene,
               digits=c(0,2),
               display=c("s", "fg"),
               caption = paste("Linear mixed effects model, gene coefficients",lobeTag))
  
  if(family){
    t3 <- xtable(mod.coef$Family,
                 digits=c(0,2),
                 display=c("s", "fg"),
                 caption = paste("Linear mixed effects model, family coefficients",lobeTag))
  } else {
    t3 = NA
  }
  
  print(t1, file=logFile, append=TRUE)
  print(t2, file=logFile, append=TRUE)
  if(family){print(t3, file=logFile, append=TRUE)}
  
  # plot random effects
  dF.plot.site <- data.frame(site = row.names(mod.coef$site),
                             effect = mod.coef$site[[1]])
  p1 <- ggplot(dF.plot.site, aes_string(x="site",y="effect"))
  p1 <- p1 + geom_bar(stat="identity") + theme_bw() + ggtitle("Effect size of scan site")
  
  dF.plot.gene <- data.frame(gene = row.names(mod.coef$gene),
                             effect = mod.coef$gene[[1]])
  p2 <- ggplot(dF.plot.gene, aes_string(x="gene",y="effect"))
  p2 <- p2 + geom_bar(stat="identity") + theme_bw() + ggtitle("Effect size of gene")
  
  if(family){
    dF.plot.Family <- data.frame(Family = row.names(mod.coef$Family),
                                 effect = mod.coef$Family[[1]])
    p3 <- ggplot(dF.plot.Family, aes_string(x="Family",y="effect"))
    p3 <- p3 + geom_bar(stat="identity") + theme_bw() + ggtitle("Effect size of family")
    p3 <- p3 + theme(axis.text.x=element_blank())
  } else {
    p3 = NA
  }
  
  # print variances
  vc <- VarCorr(mod)

  t4 <- xtable(data.frame(vc),
               display=c("s","s","s","s","g","fg"),
               digits=c(0,0,0,0,2,2),
               caption="Variance of random effects")

  print(t4,
        file=logFile,
        append=TRUE,
        include.rownames=FALSE)
  
#   # Type II ANOVA
#   mod.avo <- Anova(mod, type="II")
#   
#   t5 <- xtable(mod.avo,
#                digits=c(0,1,1,2),
#                display=c("s","f","d","fg"),
#                caption=paste("ANOVA of linear mixed effects model for", metricName))
    t5 <- xtable(summary(mod)[[10]],
                 digits=c(0,2,2,1,2,2),
                 display=c("s","fg","f","f","f","g"),
                 caption=paste("Satterthwaite estimates of pvalues of linear mixed effects model for", metricName,lobeTag))
  
  
  print(t5,
        file=logFile,
        append=TRUE)

  return(list(t1, t2, t3, t4, t5, p1, p2, p3))
}

wholeBrainAnalysisMixedEffectsHubComparison <- function(metric,
                                                         metricName,
                                                         cols,
                                                         sp, weighted=TRUE,
                                                         outDir="wholeBrainResults",
                                                         edgePC=3, tsz=12,
                                                         exclNeg=FALSE,
                                                         family=FALSE,
                                                         hubT=1.5,
                                                         h=15,w=15,s=4,ps=4  # tsz = text size
){
  # create output directory
  dir.create(outDir, showWarnings = FALSE)
  
  # define log output file
  if(weighted){
    outFile = paste(outDir,
                    paste(metric,"wt",sep="_"),
                    sep="/")
    
  } else {
    outFile = paste(outDir,
                    paste(metric,sep="_"),
                    sep="/")
  }
  logFile = paste(outFile, "logFile.tex", sep="_")
  
  # import graph data
  dF <- importGraphData(metric, weighted, edgePC, hubT=hubT)
  # apply spike percentage threshold
  dF <- applySP(dF, sp)
  dF.wb <- stackIt(dF, metric)
  dF.wb <- cbind(dF.wb, hub="hub")
  
  dF.nonhubs <- importGraphData(metric, weighted, edgePC, hubT=hubT, nonHubs=TRUE)
  # apply spike percentage threshold
  dF.nonhubs <- applySP(dF.nonhubs, sp)
  dF.wb.nonhubs <- stackIt(dF.nonhubs, metric)
  dF.wb.nonhubs <- cbind(dF.wb.nonhubs, hub="nonhub")
  
  dF.wb <- rbind(dF.wb, dF.wb.nonhubs)
  
  colList = unlist(cols[levels(dF$GS)])
  
  # Set contrasts if control group are to be excluded
  # NB this retains the control group in the estimation of the age effect and as a null regressor
  if(exclNeg){
    cMat <- matrix(c(0,-1,1), nrow = 3)
    rownames(cMat) <- levels(dF.wb$GS)
    colnames(cMat) <- c("carriersVsFTD")
    attr(dF.wb$GS, "contrasts") = cMat
  }
  
  # plot the data
  p4 <- ggplot(dF.wb, aes_string(x="GS", y="values", fill="hub"))
  p4 <- p4 + geom_boxplot()
  
  p4 <- p4 + theme_bw()
  p4 <- p4 + labs(y=metricName) + theme(axis.title.x=element_blank())
  p4 <- p4 + theme(text=element_text(size=tsz), legend.title=element_blank())
#   p4 <- p4 + scale_fill_manual(name="Group",values=as.vector(colList))
#   p4 <- p4 + theme(legend.position="none")
  
  plotFile = paste(paste(outFile,"HubComparison",sep="_"), "png", sep=".")
  ggsave(plotFile,
         scale=s,
         dpi=600,
         height=h, width=w,
         units="mm")
  
  # ANOVA of the differences
  if(family){
    mod <- lmer(values ~ GS*hub + Age + (1 | wbic) + (1 | gene) + (1 | site) + (1 | Family) + (1 | hub:wbic),
                data=dF.wb,
                REML=FALSE)
  } else {
    mod <- lmer(values ~ GS*hub + Age + (1 | wbic) + (1 | hub) + (1 | gene) + (1 | site),
                data=dF.wb,
                REML=FALSE)
  }
  
  # print random effects of the model
  mod.coef <- ranef(mod)

  t1 <- xtable(mod.coef$site,
               digits=c(0,2),
               display=c("s", "fg"),
               caption = "Linear mixed effects model, site coefficients")
  
  t2 <- xtable(mod.coef$gene,
               digits=c(0,2),
               display=c("s", "fg"),
               caption = "Linear mixed effects model, gene coefficients")
  
  if(family){
    t3 <- xtable(mod.coef$Family,
                 digits=c(0,2),
                 display=c("s", "fg"),
                 caption = "Linear mixed effects model, family coefficients")
  } else {
    t3 = NA
  }
  
  print(t1, file=logFile, append=TRUE)
  print(t2, file=logFile, append=TRUE)
  if(family){print(t3, file=logFile, append=TRUE)}
  
  # plot random effects
  dF.plot.site <- data.frame(site = row.names(mod.coef$site),
                             effect = mod.coef$site[[1]])
  p1 <- ggplot(dF.plot.site, aes_string(x="site",y="effect"))
  p1 <- p1 + geom_bar(stat="identity") + theme_bw() + ggtitle("Effect size of scan site")
  
  dF.plot.gene <- data.frame(gene = row.names(mod.coef$gene),
                             effect = mod.coef$gene[[1]])
  p2 <- ggplot(dF.plot.gene, aes_string(x="gene",y="effect"))
  p2 <- p2 + geom_bar(stat="identity") + theme_bw() + ggtitle("Effect size of gene")
  
  if(family){
    dF.plot.Family <- data.frame(Family = row.names(mod.coef$Family),
                                 effect = mod.coef$Family[[1]])
    p3 <- ggplot(dF.plot.Family, aes_string(x="Family",y="effect"))
    p3 <- p3 + geom_bar(stat="identity") + theme_bw() + ggtitle("Effect size of family")
    p3 <- p3 + theme(axis.text.x=element_blank())
  } else {
    p3 = NA
  }
  
  # print variances
  vc <- VarCorr(mod)
  
  t4 <- xtable(data.frame(vc),
               display=c("s","s","s","s","g","fg"),
               digits=c(0,0,0,0,2,2),
               caption="Variance of random effects")
  
  print(t4,
        file=logFile,
        append=TRUE,
        include.rownames=FALSE)
  
  #   # Type II ANOVA
  #   mod.avo <- Anova(mod, type="II")
  #   
  #   t5 <- xtable(mod.avo,
  #                digits=c(0,1,1,2),
  #                display=c("s","f","d","fg"),
  #                caption=paste("ANOVA of linear mixed effects model for", metricName))
  t5 <- xtable(summary(mod)[[10]],
               digits=c(0,2,2,1,2,2),
               display=c("s","fg","f","f","f","g"),
               caption=paste("Satterthwaite estimates of pvalues of linear mixed effects model for", metricName))
  
  
  print(t5,
        file=logFile,
        append=TRUE)
  
  return(list(t1, t2, t3, t4, t5, p1, p2, p3, p4))
}


graphTimeComparison <- function(metric,
                                metricName,
                                sp, # spike percentage
                                cols,
                                weighted=TRUE, # is this a weighted metric?
                                outDir="wholeBrainVsAOOResults",
                                edgePC=3,
                                h=30,w=45,s=4,tsz=12,ps=4,
                                exclNeg=TRUE,
                                family=FALSE,
                                lobe=NA,
                                hubT=NA){
  # create output directory
  dir.create(outDir, showWarnings = FALSE)
  
  # define log output file
  if(weighted){
    outFile = paste(outDir,
                    paste(metric,"wt",sep="_"),
                    sep="/")
    
  } else {
    outFile = paste(outDir,
                    paste(metric,sep="_"),
                    sep="/")
    
  }
  logFile = paste(outFile, "logFile.tex", sep="_")
  
  # create label for legends according to lobe
  if(!is.na(lobe)){
    lobeTag=lobe
  } else {
    lobeTag=""
  }
  
  # create log file
  initiateLog(logFile, metricName)
  
  ### function to plot and analyse the relationship between graph metrics and expected time to disease onset
  # import graph metric data
  dF <- importGraphData(metric, weighted, edgePC, lobe=lobe, hubT=hubT)
  
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
  dF <- unique(dF) # remove duplicates that may have sneaked in
  
  # add lobe information if necessary
  if(!is.na(lobe)){
    dF <- data.frame(dF, lobe=lobe)
  }
  
#   # create dataframe to contrast gene negative with all gene positives
#   dFgngp <- dF
#   dFgngp$GS <- as.factor(sapply(dFgngp$GS, function(x) if(x!="gene negative"){return("gene positive")}else{return(as.character(x))}))
  
  dFgn <- dF[dF$GS=="gene negative",]
  if(exclNeg){
    dF <- dF[dF$GS!="gene negative",]
    outFile = paste(outFile, "GenePos", sep="")
    dF$GS <- factor(as.character(dF$GS), levels=c("gene carriers", "FTD"))
  } #else {
#     dF$GS <- sapply(dF$GS, function(x) if(x!="gene negative"){return("gene carriers")} else {return("gene negative")})
#  }
  # dF$GS <- as.factor(as.character(dF$GS))
  
  dFgn$GS <- factor(as.character(dFgn$GS), levels=c("gene negative"))
  
  # simple linear with all groups combined
  lm.simpLM <- lm(values~aoo, data=dF)
  
  # simple linear 
  lm.all <- lm(values~aoo*GS, data=dF)
  lm.gn <- lm(values~aoo, data=dFgn)
  
  # null model - assumes no interaction between onset and gene status
  lm.null <- lm(values ~ aoo+GS, data=dF)
  
  # plot linear model NOT MIXED EFFECTS
  dF.plot <- dF

  p <- ggplot(dF.plot, aes_string(x="aoo", y="values", colour="GS", group="GS"))
  # p <- p + geom_point(size=ps)
  p <- p + scale_x_continuous(breaks=seq(-100,100,10))
  colList = unlist(cols[levels(dF.plot$GS)])

  # plot simple linear model with all groups combined
  p5 <- ggplot(dF.plot, aes_string(x="aoo", y="values"))
  p5 <- p5 + scale_x_continuous(breaks=seq(-100,100,10))
  p5 <- p5 + geom_smooth(method=lm)
#   int <- lm.simpLM[["coefficients"]][[1]]
#   slope <- lm.simpLM[["coefficients"]][[2]]
#   p5 <- p5 + geom_abline(intercept=int, slope=slope)
  
  # from line below: paste(metricName, "linear regression", sep="\n")
  p5 <- p5 + labs(title=lobeTag, y=metricName, x="Estimated age of onset") + theme(axis.title.x=element_blank()) #, legend.key=element_rect(fill="white", colour="white"))
  p5 <- p5 + theme_bw() # + theme(legend.key = element_rect(colour="#FFFFFF", fill = "#FFFFFF"))
  p5 <- p5 + theme(text=element_text(size=tsz), plot.title=element_text(size=tsz))

  plotOutName = paste(paste(outFile, "simpLM", sep="_"),"png",sep=".")

  ggsave(plotOutName,
         plot=p5,
         scale=s,
         dpi=600,
         height=h, width=w,
         units="mm")
  
  t7 <- xtable(summary(lm.simpLM),
               caption=paste("summary of linear model combining all groups for", metricName, lobeTag),
               digits=c(0,2,2,2,2),
               display=c("s","fg","fg","fg","g"))
  
  p1 <- p + geom_smooth(method="lm")
  coeffs <- summary(lm.all)[["coefficients"]]
  xMin = min(dF[dF$GS=="gene carriers","aoo"])
  xMax = max(dF[dF$GS=="gene carriers","aoo"])
  yMin = xMin*coeffs["aoo","Estimate"] + coeffs["(Intercept)","Estimate"] + coeffs["aoo:GSFTD","Estimate"] 
  yMax = xMax*coeffs["aoo","Estimate"] + coeffs["(Intercept)","Estimate"] + coeffs["aoo:GSFTD","Estimate"]

  xMin.FTD = min(dF[dF$GS=="FTD","aoo"])
  xMax.FTD = max(dF[dF$GS=="FTD","aoo"])

  yMin.FTD = xMin.FTD*coeffs["aoo:GSFTD","Estimate"] + xMin.FTD*coeffs["aoo","Estimate"] + coeffs["(Intercept)","Estimate"] + coeffs["GSFTD","Estimate"]
  yMax.FTD = xMax.FTD*coeffs["aoo:GSFTD","Estimate"] + xMax.FTD*coeffs["aoo","Estimate"] + coeffs["(Intercept)","Estimate"] + coeffs["GSFTD","Estimate"]
 
  p1 <- p1 + geom_segment(x=xMin,
                         xend=xMax,
                         y=yMin,
                         yend=yMax,
                         size=1.5,
                         colour=cols["gene carriers"])
  
  p1 <- p1 + geom_segment(x=xMin.FTD,
                         xend=xMax.FTD,
                         y=yMin.FTD,
                         yend=yMax.FTD,
                         size=1.5,
                         colour=cols["FTD"])

  yMin.scale = min(yMin,yMin.FTD,yMax,yMax.FTD)
  yMax.scale = max(yMax,yMax.FTD,yMin,yMin.FTD)
  p1 <- p1 + scale_y_continuous(limits=c(yMin.scale, yMax.scale))
  
  p1 <- p1 + scale_x_continuous(limits=c(min(dF$aoo), max(dF$aoo)))

  p1 <- p1 + scale_colour_manual(name="Group",values=as.vector(colList))
  # from line below: paste(metricName, "linear regression", sep="\n")
  p1 <- p1 + labs(title=lobeTag, y=metricName, x="Estimated age of onset") + theme(axis.title.x=element_blank(), legend.key=element_rect(fill="white", colour="white"))
  p1 <- p1 + theme_bw() + theme(legend.key = element_rect(colour="#FFFFFF", fill = "#FFFFFF"))
  p1 <- p1 + theme(text=element_text(size=tsz), plot.title=element_text(size=tsz))
  plotOutName = paste(outFile,"png",sep=".")

  ggsave(plotOutName,
         plot=p1,
         scale=s,
         dpi=600,
         height=h, width=w,
         units="mm")
  
  t1 <- xtable(summary(lm.all),
               caption=paste("summary of linear model for the interaction between estimated age of onset (aoo) and gene status (GS) on", metricName, lobeTag),
               digits=c(0,2,2,2,2),
               display=c("s","fg","fg","fg","g"))
  
  print(t1,
        include.rownames=TRUE,
        file=logFile,
        append=TRUE)
  
  # comparison with null model
  aov.lm.null <- anova(lm.all, lm.null)
  
  t2 <- xtable(aov.lm.null,
               display=c("s","d","fg","d","fg","fg","fg"),
               digits = c(0, 0, 2, 0, 2, 2, 2),
               caption = paste("Assessment of whether there is a difference between slopes of gene positive and gene affected groups by comparing models including and excluding and interaction between age at onset and gene status", lobeTag))
  
  print(t2,
        include.rownames=FALSE,
        file=logFile,
        append=TRUE)
  
  # now run linear mixed effects model model
  # run model
  if(exclNeg){
    cMat <- matrix(c(-1,1), nrow = 2)
    rownames(cMat) <- levels(dF$GS)
    colnames(cMat) <- c("carriersVsFTD")
    attr(dF$GS, "contrasts") = cMat
  } else {
    cMat <- matrix(c(0,-1,1), nrow = 3)
    rownames(cMat) <- levels(dF$GS)
    colnames(cMat) <- c("carriersVsFTD")
    attr(dF$GS, "contrasts") = cMat
  }

  if(family){
    mod <- lmer(values ~ aoo * GS + (1 | gene) + (1 | site) + (1 | Family),
                data=dF,
                REML=FALSE)
  } else {
    mod <- lmer(values ~ aoo * GS + (1 | gene) + (1 | site),
                data=dF,
                REML=FALSE)
  }

  ## plot random effects
  mod.coef <- ranef(mod)
  dF.plot.site <- data.frame(site = row.names(mod.coef$site),
                             effect = mod.coef$site[[1]])
  p2 <- ggplot(dF.plot.site, aes_string(x="site",y="effect"))
  p2 <- p2 + geom_bar(stat="identity") + theme_bw() + ggtitle("Effect size of scan site")
  
  dF.plot.gene <- data.frame(gene = row.names(mod.coef$gene),
                             effect = mod.coef$gene[[1]])
  p3 <- ggplot(dF.plot.gene, aes_string(x="gene",y="effect"))
  p3 <- p3 + geom_bar(stat="identity") + theme_bw() + ggtitle("Effect size of gene")
  
  if(family){
    dF.plot.Family <- data.frame(Family = row.names(mod.coef$Family),
                                 effect = mod.coef$Family[[1]])
    p4 <- ggplot(dF.plot.Family, aes_string(x="Family",y="effect"))
    p4 <- p4 + geom_bar(stat="identity") + theme_bw() + ggtitle("Effect size of family")
    p4 <- p4 + theme(axis.text.x=element_blank())
  } else {
    p4 = NA
  }
  
  t3 <- xtable(summary(mod)[["coefficients"]],
               caption=paste("Linear mixed effects model, fixed effects", lobeTag),
               digits=c(0,2,2,2,2,2),
               display=c("s","fg","fg","fg","fg","fg"))

  print(t3,
        include.rownames=TRUE,
        file=logFile,
        append=TRUE)
  
  t4 <- xtable(data.frame(StdDev = c(attributes(VarCorr(mod)[[1]])[["stddev"]],
                                     attributes(VarCorr(mod)[[2]])[["stddev"]],
                                     attributes(VarCorr(mod))[["sc"]]),
                          row.names = c("Family", "Site", "Residual")),                          
               caption=paste("Linear mixed effects model, random effects", lobeTag),
               digits=c(0,2),
               display=c("s","fg"))
  
  print(t4,
        include.rownames=TRUE,
        file=logFile,
        append=TRUE)
  
  # print the variance and standard deviation of the random effects
  vc <- VarCorr(mod)
  
  t5 <- xtable(data.frame(vc),
               display=c("s","s","s","s","g","fg"),
               digits=c(0,0,0,0,2,2),
               caption=paste("Variance of random effects", lobeTag))
  
  print(t5,
        file=logFile,
        append=TRUE,
        include.rownames=FALSE)
  
  # null model
  # this assumes there is no interaction between the gene status and age of onset
  if(family){
    nulMod <- lmer(values ~ aoo + GS + Age + (1 | gene) + (1 | site) + (1 | Family),
                   data=dF,
                   REML=FALSE)
  } else {
    nulMod <- lmer(values ~ aoo + GS + Age + (1 | gene) + (1 | site),
                   data=dF,
                   REML=FALSE)
  }
  
  modComparison <- anova(mod,nulMod)
  
  t6 <- xtable(modComparison,
               caption = paste("Assessment of whether there is a difference between slopes of gene positive and gene affected groups by comparing models including and excluding and interaction between age at onset and gene status", lobeTag),
               digits = c(0,0,2,2,2,2,2,2,2),
               display = c("s","d","fg","fg","fg","fg","d","fg","fg"))
  
  print(t6,
        include.rownames=TRUE,
        file=logFile,
        append=TRUE)
  
  # fixed effects
  t8 <- xtable(summary(mod)[[10]],
                    digits=c(0,2,2,1,2,2),
                    display=c("s","fg","f","f","f","g"),
                    caption=paste("Satterthwaite estimates of pvalues of linear mixed effects model for", metricName, lobeTag))
  
  # now run the mixed effects model comparing gene negative with gene positive group
  if(!exclNeg){
    cMat <- matrix(c(0,1,1), nrow = 3)
    rownames(cMat) <- levels(dF$GS)
    colnames(cMat) <- c("geneNegVsGenePos")
    attr(dF$GS, "contrasts") = cMat

    if(family){
      mod <- lmer(values ~ aoo * GS + (1 | gene) + (1 | site) + (1 | Family),
                  data=dF,
                  REML=FALSE)
    } else {
      mod <- lmer(values ~ aoo * GS + (1 | gene) + (1 | site),
                  data=dF,
                  REML=FALSE)
    }
    
    
    # plot random effects
    mod.coef <- ranef(mod)
    dF.plot.site <- data.frame(site = row.names(mod.coef$site),
                               effect = mod.coef$site[[1]])
    p6 <- ggplot(dF.plot.site, aes_string(x="site",y="effect"))
    p6 <- p6 + geom_bar(stat="identity") + theme_bw() + ggtitle("Effect size of scan site")
    
    dF.plot.gene <- data.frame(gene = row.names(mod.coef$gene),
                               effect = mod.coef$gene[[1]])
    p7 <- ggplot(dF.plot.gene, aes_string(x="gene",y="effect"))
    p7 <- p7 + geom_bar(stat="identity") + theme_bw() + ggtitle("Effect size of gene")
    
    if(family){
      dF.plot.Family <- data.frame(Family = row.names(mod.coef$Family),
                                   effect = mod.coef$Family[[1]])
      p8 <- ggplot(dF.plot.Family, aes_string(x="Family",y="effect"))
      p8 <- p8 + geom_bar(stat="identity") + theme_bw() + ggtitle("Effect size of family")
      p8 <- p8 + theme(axis.text.x=element_blank())
    } else {
      p8 = NA
    }
    
    t9 <- xtable(summary(mod)[["coefficients"]],
                 caption=paste("Linear mixed effects model, fixed effects", lobeTag),
                 digits=c(0,2,2,2,2,2),
                 display=c("s","fg","fg","fg","fg","fg"))
    
    print(t9,
          include.rownames=TRUE,
          file=logFile,
          append=TRUE)
    
    t10 <- xtable(data.frame(StdDev = c(attributes(VarCorr(mod)[[1]])[["stddev"]],
                                       attributes(VarCorr(mod)[[2]])[["stddev"]],
                                       attributes(VarCorr(mod))[["sc"]]),
                            row.names = c("Family", "Site", "Residual")),                          
                 caption=paste("Linear mixed effects model, random effects", lobeTag),
                 digits=c(0,2),
                 display=c("s","fg"))
    
    print(t10,
          include.rownames=TRUE,
          file=logFile,
          append=TRUE)
    
    # print the variance and standard deviation of the random effects
    vc <- VarCorr(mod)
    
    t11 <- xtable(data.frame(vc),
                 display=c("s","s","s","s","g","fg"),
                 digits=c(0,0,0,0,2,2),
                 caption=paste("Variance of random effects", lobeTag))
    
    print(t11,
          file=logFile,
          append=TRUE,
          include.rownames=FALSE)
    
    # null model
    # this assumes there is no interaction between the gene status and age of onset
    if(family){
      nulMod <- lmer(values ~ aoo + GS + Age + (1 | gene) + (1 | site) + (1 | Family),
                     data=dF,
                     REML=FALSE)
    } else {
      nulMod <- lmer(values ~ aoo + GS + Age + (1 | gene) + (1 | site),
                     data=dF,
                     REML=FALSE)
    }
    
    modComparison <- anova(mod,nulMod)
    
    t12 <- xtable(modComparison,
                 caption = paste("Assessment of whether there is a difference between slopes of gene positive and gene affected groups by comparing models including and excluding and interaction between age at onset and gene status", lobeTag),
                 digits = c(0,0,2,2,2,2,2,2,2),
                 display = c("s","d","fg","fg","fg","fg","d","fg","fg"))
    
    print(t12,
          include.rownames=TRUE,
          file=logFile,
          append=TRUE)
    
    # fixed effects
    t13 <- xtable(summary(mod)[[10]],
                 digits=c(0,2,2,1,2,2),
                 display=c("s","fg","f","f","f","g"),
                 caption=paste("Satterthwaite estimates of pvalues of linear mixed effects model for", metricName, lobeTag))
  } else {
    p6=NA
    p7=NA
    p8=NA
    t9=NA
    t10=NA
    t11=NA
    t12=NA
  }

  return(list(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12,
              p1, p2, p3, p4, p5, p6, p7, p8))
}
  
graphTimeComparisonNL <- function(metric,
                                metricName,
                                sp, # spike percentage
                                cols,
                                startvec=NULL, # starting vectors for equation with quadratic term
                                startvecCub=c(dc=0.000000001), # starting vectors for equation with cubic term
                                weighted=TRUE, # is this a weighted metric?
                                outDir="wholeBrainVsAOOResults",
                                edgePC=3,
                                h=30,w=45,s=4,tsz=12,ps=4,
                                sink=TRUE){
  # create output directory
  dir.create(outDir, showWarnings = FALSE)
  
  # define log output file
  if(weighted){
    outFile = paste(outDir,
                    paste(metric,"wt",sep="_"),
                    sep="/")
    
  } else {
    outFile = paste(outDir,
                    paste(metric,sep="_"),
                    sep="/")
    
  }
  logFile = paste(outFile, "logFile.tex", sep="_")
  
  # create log file
  initiateLog(logFile, metricName)
  
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
  dF <- unique(dF) # remove duplicates that may have sneaked in
  
  # now run non-linear mixed effects model with a quadratic term
  nloutFile = paste(outDir,
                    paste(metric,"NLlogFile.txt",sep="_"),
                    sep="/")
  if(sink){
    sink(nllogFile)
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
      # p.temp <- p.temp + geom_point(data=dF.temp, aes_string(x="x", y="y"))
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
  # p <- p + geom_point(size=ps)
  
  x = seq(-45,25,by=1)
  dF.nlmq <- data.frame(x=x, y=SSquadFun(x, nlmq@beta[[1]], nlmq@beta[[2]], nlmq@beta[[3]]), GS="estimates")
  
  pq <- p + geom_line(data=dF.nlmq, aes_string(x="x", y="y"))
  cList = c("FTD", "estimates", "gene carriers")
  colList = unlist(cols[cList])
  pq <- pq + scale_colour_manual(name="Group",values=as.vector(colList),
                                 breaks=c("gene carriers", "FTD", "estimates"))
  pq <- pq + theme_bw() + theme(legend.key=element_rect(fill="white", colour="white"))
  pq <- pq + labs(title=paste(metricName, "non-linear regression, quadratic", sep="\n"), y=metricName, x="Estimated time from onset")# + theme(axis.title.x=element_blank())
  pq <- pq + theme(text=element_text(size=tsz), plot.title=element_text(size=tsz), plot.title=element_text(size=tsz))
  
#   print(paste("Saving",paste(outDir,
#                        paste(paste(metric, "nonLinearEstimatesQuadratic",sep="_"),"png",sep="."),
#                        sep="/")))
  plotFile = paste(paste(outFile, "nonLinearEstimatesQuadratic",sep="_"), "png", sep=".")
  ggsave(plotFile,
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
      # p.temp <- pp + geom_point(data=dF.temp, aes_string(x="x", y="y"))
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
                                 breaks=c("gene carriers", "FTD", "estimates"))
  pc <- pc + theme_bw() + theme(legend.key=element_rect(fill="white", colour="white"))
  pc <- pc + labs(title=paste(metricName, "non-linear regression, cubic", sep="\n"), y=metricName, x="Estimated time from onset")# + theme(axis.title.x=element_blank())
  pc <- pc + theme(text=element_text(size=tsz), plot.title=element_text(size=tsz))
  
#   print(paste("Saving",paste(outDir,
#                              paste(paste(metric, "nonLinearEstimatesCubic",sep="_"),"png",sep="."),
#                              sep="/")))
  plotFile = paste(paste(outFile, "nonLinearEstimatesCubic",sep="_"),"png",sep=".")
  ggsave(plotFile,
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
          file=logFile,
          append=TRUE)
  }
  
  # finalise log file and return dataframe
  endLog(logFile)
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
             "gene carriers"="#56B4E9",
             "FTD"="#D55E00",
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
                       outDir="wholeBrainVsAOOResults",
                       weighted=FALSE,
                       tsz=12,ps=4){
  # define log output file
  if(weighted){
    outFile = paste(outDir,
                    paste(metric,"wt",sep="_"),
                    sep="/")
    
  } else {
    outFile = paste(outDir,
                    paste(metric,sep="_"),
                    sep="/")
    
  }
  
  logFile = paste(outFile, "logFile.tex", sep="_")
  
  # create output directory
  dir.create(outDir, showWarnings = FALSE)
  
  # create log file
  initiateLog(logFile, metricName)
  
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
  dF <- unique(dF) # remove duplicates that may have sneaked in
  
#   # filter out affected subjects less than t=0 and gene positive greater than t=0
#   dF <- rbind(dF[dF$GS=="FTD" && dF$aoo>=0,],
#               dF[dF$GS=="gene carriers" && dF$aoo<=0,])
  
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
  yMax = max(max(dF$aoo)*br$coefficients[[3]]+br$coefficients[[1]])
  yMin = min(dF$aoo)*br$coefficients[[2]]+br$coefficients[[1]]
  p <- ggplot(dF, aes_string(x="aoo", y="values", colour="GS"))
  # p <- p + geom_point(size=ps)
  p <- p + geom_segment(x=min(dF$aoo),
                        xend=br$psi[[2]],
                        y=yMin,
                        yend=y,
                        colour="black")
  
  p <- p + geom_segment(x=br$psi[[2]],
                        xend=max(dF$aoo),
                        y=y,
                        yend=yMax,
                        colour="black")
  
  # specify y axis
  yGap = yMax - yMin
  yMin.axis = yMin - 0.05*yGap
  yMax.axis = yMax + 0.05*yGap
  p <- p + scale_y_continuous(limits=c(yMin.axis, yMax.axis))
  
  colList = unlist(cols[levels(dF$GS)])
  p <- p + scale_colour_manual(name="Group",values=as.vector(colList))
  p <- p + labs(title=paste(metricName, "breakpoint analysis", sep="\n"), y=metricName, x="Estimated age of onset") + theme(axis.title.x=element_blank(), legend.key=element_rect(fill="white", colour="white"))
  p <- p + theme_bw() + theme(legend.key = element_rect(colour="#FFFFFF", fill = "#FFFFFF"))
  p <- p + theme(text=element_text(size=tsz), plot.title=element_text(size=tsz))
  
  return(list(br, p, pVal))
}


breakPointDiscontinuous <- function(metric,
                                    metricName,
                                    sp, # spike percentage
                                    cols,
                                    edgePC=3,
                                    outDir="wholeBrainVsAOOResults",
                                    weighted=FALSE,
                                    exclNeg=FALSE,
                                    h=30,w=45,s=4,tsz=12,ps=4,
                                    normalise=TRUE,
                                    lobe=NA,
                                    hubT=NA){
  # define log output file
  if(weighted){
    outFile = paste(outDir,
                    paste(metric,"wt",sep="_"),
                    sep="/")
    
  } else {
    outFile = paste(outDir,
                    paste(metric,sep="_"),
                    sep="/")
    
  }
  
  logFile = paste(outFile, "logFile.tex", sep="_")
  
  # create label for legends according to lobe
  if(!is.na(lobe)){
    lobeTag=lobe
  } else {
    lobeTag=""
  }
  
  # create output directory
  dir.create(outDir, showWarnings = FALSE)
  
  # create log file
  initiateLog(logFile, metricName)
  
  ### function to plot and analyse the relationship between graph metrics and expected time to disease onset
  # import graph metric data
  dF <- importGraphData(metric, weighted, edgePC, lobe=lobe, hubT=hubT)
  
  # filter by spike percentage
  dF <- applySP(dF, sp)
  
  # stack the data and take mean if a nodewise measure
  if(weighted){
    metric = paste(metric, "wt", sep="_")
  }
  
  dF <- stackIt(dF, metric) # changes the name of the metric to 'values'
  
  # get age of onset data
  genfiData <- read.table("/home/tim/GENFI/genfi_Subjects_sjones_1_22_2015_17_47_47_restructure_summary.csv",
                          sep="\t",
                          header = TRUE)
  
  dF.aoo <- data.frame(wbic=genfiData$Subject,
                       aoo=genfiData$Yrs.from.AV_AAO)
  
  # merge age of onset data
  dF <- merge(dF, dF.aoo, by="wbic")
  dF <- unique(dF) # remove duplicates that may have sneaked in
  
  # make small adjustment for a set of twins in family 172
  dup <- which(duplicated(dF$aoo))
  dF[dup,"aoo"] <- dF[dup,"aoo"] + 0.001
  
  # filter out gene negative subjects
  if(exclNeg==FALSE){
    outFile = paste(outFile, "GenePos", sep="")
  }
  
  dFgn <- dF[dF$GS=="gene negative",]
  dF <- dF[dF$GS!="gene negative",]
  dF$GS <- factor(as.character(dF$GS), levels=c("gene carriers", "FTD"))
  
  # calculate z scores
  sd.geneNeg <- sd(dFgn[,"values"])
  mean.geneNeg <- mean(dFgn[,"values"])
  dF <- cbind(dF, values.norm=sapply(dF$values, function(x) (x-mean.geneNeg)/sd.geneNeg))
  
  #   # filter out affected subjects less than t=0 and gene positive greater than t=0
  #   dF <- rbind(dF[dF$GS=="FTD" && dF$aoo>=0,],
  #               dF[dF$GS=="gene carriers" && dF$aoo<=0,])
  
  # define breakpoint
  bkpt = 0.
  
#   piecewise <- lmer(values ~ aoo*(aoo<bkpt) + aoo*(aoo>bkpt) + (1 | gene) + (1 | site) + (1 | Family),
#                     data=dF,
#                     REML=FALSE)
  piecewise <- lm (values ~ aoo*(aoo<bkpt) + aoo*(aoo>=bkpt), data=dF)

  t1 <- xtable(summary(piecewise),
               caption = paste("Summary of discontinuous piecewise regression analysis",lobeTag),
               digits = c(0,2,2,2,2),
               display = c("s","fg","fg","fg","g"))
  print(t1,
        include.rownames=TRUE,
        file=logFile,
        append=TRUE)
  
  # compare with null model
#   nulMod <- lmer(values ~ aoo + (1 | gene) + (1 | site) + (1 | Family),
#                   data=dF,
#                   REML=FALSE)
  nulMod <- lm(values ~ aoo,
               data=dF)
  
  modComparison <- anova(piecewise,nulMod)
  
  t2 <- xtable(modComparison,
               caption = paste("ANOVA between model and null model to assess whether the segmented model fits better",lobeTag),
               digits = c(0,0,2,2,2,2,2),
               display = c("s","d","fg","fg","fg","fg","fg"))
  
  print(t2,
        include.rownames=TRUE,
        file=logFile,
        append=TRUE)
  
  # control linear model
  gnlm <- lm(values ~ aoo, data=dFgn)
  int <- gnlm[["coefficients"]][[1]]
  slope <- gnlm[["coefficients"]][[2]]
  
  # get confidence intervals
  piecewise.CIs <- confint(piecewise)
  
  # plot breakpoint data
  y <- piecewise$coefficients[[1]] # get the y value of the piecewiseeakpoint
  p <- ggplot(dF, aes_string(x="aoo", y="values.norm"))
  
  if(exclNeg==FALSE){
    yVal = min(dF$aoo)*gnlm$coefficients[[2]] + gnlm$coefficients[[1]]
    yendVal = max(dF$aoo)*gnlm$coefficients[[2]] + gnlm$coefficients[[1]]

    if(!normalise){
      p <- p + geom_segment(x=min(dFgn$aoo),
                            xend=max(dFgn$aoo),
                            y=yVal,
                            yend=yendVal,
                            size=1.5,
                            colour="#E69F00",
                            linetype="longdash")
      xmin = min(dFgn$aoo)
      
    } else {
      xmin = min(dF$aoo)
    }
  }
  
  # set x axis limits
  xmax = max(dF$aoo)
  p <- p + scale_x_continuous(breaks = seq(-100,100,10), limits=c(xmin,xmax))
  
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

  # set colours
  colList = unlist(cols[levels(dF$GS)])
  
#   #significance areas
#   p <- p + geom_ribbon(ymin=seq(yMin.lower, yEndVal.lower,length.out = length(dF[,1])),
#                        ymax=seq(yMin.upper, yEndVal.upper,length.out = length(dF[,1])),
#                        x=seq(min(dF$aoo),0., length.out = length(dF[,1])),
#                        fill=colList["gene carriers"], alpha=0.5)
#   
#   p <- p + geom_ribbon(ymin=seq(yVal.lower, yMax.lower,length.out = length(dF[,1])),
#                        ymax=seq(yVal.upper, yMax.upper,length.out = length(dF[,1])),
#                        x=seq(0., max(dF$aoo), length.out = length(dF[,1])),
#                        fill=colList["FTD"], alpha=0.5)
#   
  # p <- p + geom_point(size=ps, shape=16, aes_string(colour="GS"))
  p <- p + geom_segment(x=min(dF$aoo),
                        xend=0.,
                        y=yMin,
                        yend=yEndVal,
                        size=1.5,
                        colour="black")
  
  p <- p + geom_segment(x=0.,
                        xend=max(dF$aoo),
                        y=yVal,
                        yend=yMax,
                        size=1.5,
                        colour="black")
  
  # add confidence intervals
  p <- p + geom_smooth(method="lm",
                       linetype="blank",
                       data=dF[dF$aoo<=0.,],
                       fill=colList["gene carriers"], alpha=0.5)
  
  p <- p + geom_smooth(method="lm",
                       linetype="blank",
                       data=dF[dF$aoo>0.,],
                       fill=colList["FTD"], alpha=0.5)
  
#   p <- p + geom_segment(x=min(dF$aoo),
#                         xend=0,
#                         y=yMin.lower,
#                         yend=yEndVal.lower,
#                         size=1.5,
#                         colour="red")
#   
#   p <- p + geom_segment(x=0.,
#                         xend=max(dF$aoo),
#                         y=yVal.lower,
#                         yend=yMax.lower,
#                         size=1.5,
#                         colour="red")
#   
#   p <- p + geom_segment(x=min(dF$aoo),
#                         xend=0,
#                         y=yMin.upper,
#                         yend=yEndVal.upper,
#                         size=1.5,
#                         colour="red")
#   
#   p <- p + geom_segment(x=0.,
#                         xend=max(dF$aoo),
#                         y=yVal.upper,
#                         yend=yMax.upper,
#                         size=1.5,
#                         colour="red")
#   

#   # specify y axis
#   yList <- c(yMin,yEndVal,yMax,yVal,
#              yEndVal.lower, yEndVal.upper,
#              yMin.lower,yMin.upper,
#              yMax.lower,yMax.upper,
#              yVal.lower,yVal.upper,
#              yEndVal.lower, yEndVal.upper,
#              yMin.lower,yMin.upper,
#              yMax.lower,yMax.upper,
#              yVal.lower,yVal.upper)
# 
#   yLow <- min(yList)
#   yHigh <- max(yList)
#   
#   yGap = yHigh-yLow
#   yMin.axis = yLow - 0.05*yGap
#   yMax.axis = yHigh + 0.05*yGap
# 
#   yMin.lab <- round(yMin.axis, digits=1)
#   yMax.lab <- round(yMax.axis, digits=1)
#   
#   if(yGap>1){
#     dg = 0
#   } else if(yGap>0.1){
#     dg = 1
#   } else if(yGap>0.01){
#     dg = 2
#   } else if(yGap>0.001){
#     dg=3
#   }
#    
#   byVal = round((yMax.lab - yMin.lab)/8, digits=dg)
#   if(byVal==0){byVal=1}
# 
#   breakList = seq(yMin.lab-byVal, yMax.lab+byVal, by=byVal)
#   
#   # centre on 0
#   breakList = breakList-min(sapply(breakList, function(x) sqrt(x^2)))
#   
#   p <- p + scale_y_continuous(limits=c(yMin.axis, yMax.axis),
#                               breaks = breakList)
#   
  # from line below: title=paste(metricName, "discontinuous breakpoint analysis", sep="\n"), 
  p <- p + labs(title=lobeTag, y=paste(metricName, "(z-score)"), x="Years from estimated age of onset") + theme(axis.title.x=element_blank(), legend.key=element_rect(fill="white", colour="white"))
  p <- p + theme_bw() + theme(legend.key = element_rect(colour="#FFFFFF", fill = "#FFFFFF"))
  # p <- p + scale_fill_manual(name="Group",values=as.vector(colList))
  p <- p + scale_colour_manual(name="Group", values=as.vector(colList))
  p <- p + theme(text=element_text(size=tsz), plot.title=element_text(size=tsz))
  
  plotFile = paste(paste(outFile,"DiscontBkpt",sep="_"),"png",sep=".")
  ggsave(plotFile,
         plot=p,
         scale=s,
         dpi=600,
         height=h, width=w,
         units="mm")
  
  ### need to finish plot ###
  return(list(t1, t2 ,p))
}

countIt <- function(x,n){return(length(x[x==n]))}

demographics <- function(demog, typeList, sp=NA, exclNeg=FALSE){
  # import graph data
  dF <- importGraphData("geNorm", FALSE, Age=FALSE)
  if(exclNeg){
    dF <- dF[dF$GS!="gene negative",]
  }
  
  dF <- unique(dF) # remove duplicates that may have sneaked in
  
  if(is.data.frame(sp)){dF <- applySP(dF, sp)}
  
  type=typeList[demog]

  # merge with demographic data
  dF.demog <- read.table("../genfi_Subjects_sjones_1_22_2015_17_47_47_restructure_summary.csv", sep="\t", header = TRUE)
  names(dF.demog)[4] <- "wbic"
  dF.demog <- dF.demog[dF.demog$duplicate!=1,] # remove duplicates that may have sneaked in
  
  dF.demog <- dF.demog[which(dF.demog$wbic %in% dF$wbic),]
  
  dF <- merge(dF, dF.demog[,-which(names(dF.demog)=="GS")],
              by = c("wbic"),
              all.x=TRUE,
              all.y=FALSE)
  names(dF)[which(names(dF)==demog)] <- "demog"
  
  # create summary table of demographic data
  if(type=="continuous"){
    dF.wb = ddply(dF, .(GS), summarise,
                  mean = mean(demog, na.rm = TRUE),
                  sd = sd(demog, na.rm = TRUE)
    )
    row.names(dF.wb) <- dF.wb[,"GS"]
    dF.wb <- dF.wb[,-1]
  } else if(type=="ordinal"){
    dF.wb = ddply(dF, .(GS), summarise,
                  count0 = countIt(demog,0),
                  count1 = countIt(demog,1),
                  count2 = countIt(demog,2)
    )
    row.names(dF.wb) <- dF.wb[,"GS"]
    dF.wb <- dF.wb[,-1]
  }

  # do statistical test
  if(type=="continuous"){
    # ANOVA
    mod <- lm(demog~GS,data=dF)
    mod.aov <- Anova(mod, type="II")
    
    # post-hoc t-tests
    if(!exclNeg){
      negVsFTD <- t.test(dF[dF$GS=="gene negative","demog"], dF[dF$GS=="FTD","demog"])
      negVsCarrier <- t.test(dF[dF$GS=="gene negative","demog"], dF[dF$GS=="gene carriers","demog"])
    } else {
      negVsFTD = list(p.value=NA,
                      statistic=NA)
      negVsCarrier = list(p.value=NA,
                          statistic=NA)
    }
    FTDVsCarrier <- t.test(dF[dF$GS=="FTD","demog"], dF[dF$GS=="gene carriers","demog"])
    
    if(exclNeg){
      return(list(Demographic=demog,
                  DOF=mod.aov[["Df"]][[1]],
                  pVal=mod.aov[["Pr(>F)"]][[1]],
                  Fval=mod.aov[["F value"]][[1]],
                  ChiSq=NA,
                  geneCarriers=paste(round(dF.wb[["gene carriers", "mean"]], digits=1),
                                     " (",round(dF.wb[["gene carriers", "sd"]], digits=1),")",
                                     sep=""),
                  FTD=paste(round(dF.wb[["FTD", "mean"]], digits=1),
                            " (",round(dF.wb[["FTD", "sd"]], digits=1),")",
                            sep=""),
                  # negVsFTD.stat=negVsFTD[["statistic"]],
                  # negVsFTD.p=negVsFTD[["p.value"]],
                  # negVsCarrier.stat=negVsCarrier[["statistic"]],
                  # negVsCarrier.p=negVsCarrier[["p.value"]],
                  FTDVsCarrier.stat=FTDVsCarrier[["statistic"]],
                  FTDVsCarrier.p=FTDVsCarrier[["p.value"]]
      )
      )
    } else {
      return(list(Demographic=demog,
                  DOF=mod.aov[["Df"]][[1]],
                  pVal=mod.aov[["Pr(>F)"]][[1]],
                  Fval=mod.aov[["F value"]][[1]],
                  ChiSq=NA,
                  geneNeg=paste(round(dF.wb[["gene negative", "mean"]], digits=1),
                                " (",round(dF.wb[["gene negative", "sd"]], digits=1),")",
                                sep=""),
                  geneCarriers=paste(round(dF.wb[["gene carriers", "mean"]], digits=1),
                                     " (",round(dF.wb[["gene carriers", "sd"]], digits=1),")",
                                     sep=""),
                  FTD=paste(round(dF.wb[["FTD", "mean"]], digits=1),
                            " (",round(dF.wb[["FTD", "sd"]], digits=1),")",
                            sep=""),
                  negVsFTD.stat=negVsFTD[["statistic"]],
                  negVsFTD.p=negVsFTD[["p.value"]],
                  negVsCarrier.stat=negVsCarrier[["statistic"]],
                  negVsCarrier.p=negVsCarrier[["p.value"]],
                  FTDVsCarrier.stat=FTDVsCarrier[["statistic"]],
                  FTDVsCarrier.p=FTDVsCarrier[["p.value"]]
                  )
             )
    }
  } else if(type=="ordinal"){
    if(sum(dF.wb[,"count2"]!=0)){ # in the event there are three categories of the demographic
      
      if(!exclNeg){ # if negatives are excluded
        geneNeg=paste(dF.wb[["gene negative", "count0"]],
                      dF.wb[["gene negative", "count1"]],
                      dF.wb[["gene negative", "count2"]],
                      sep="/")
        
        # post-hoc chi-square tests
        negVsFTD <- chisq.test(dF.wb[c("FTD","gene negative"),])
        negVsCarrier <- chisq.test(dF.wb[c("gene carriers","gene negative"),])
        
      } else { # if negatives are included
        geneNeg=NA
        negVsFTD = list(p.value=NA,
                        statistic=NA)
        negVsCarrier = list(p.value=NA,
                            statistic=NA)
      }
      geneCarriers=paste(dF.wb[["gene carriers", "count0"]],
                         dF.wb[["gene carriers", "count1"]],
                         dF.wb[["gene carriers", "count2"]],
                         sep="/")
      
      FTD=paste(dF.wb[["FTD", "count0"]],
                dF.wb[["FTD", "count1"]],
                dF.wb[["FTD", "count2"]],
                sep="/")
      
    } else { # in the event there are only two categorical variables
      dF.wb <- dF.wb[,-3]
      if(!exclNeg){
        geneNeg=paste(dF.wb[["gene negative", "count0"]],
                      dF.wb[["gene negative", "count1"]],
                      sep="/")
        
        # post-hoc chi-square tests
        negVsFTD <- chisq.test(dF.wb[c("FTD","gene negative"),])
        negVsCarrier <- chisq.test(dF.wb[c("gene carriers","gene negative"),])
        
      } else {
        geneNeg=NA
        negVsFTD = list(p.value=NA,
                        statistic=NA)
        negVsCarrier = list(p.value=NA,
                            statistic=NA)
      }
      
      geneCarriers=paste(dF.wb[["gene carriers", "count0"]],
                         dF.wb[["gene carriers", "count1"]],
                         sep="/")
      FTD=paste(dF.wb[["FTD", "count0"]],
                dF.wb[["FTD", "count1"]],
                sep="/")
    }
    dF.cst <- chisq.test(dF.wb) # perform chi-square test

    # post-hoc chi-square
    FTDVsCarrier <- chisq.test(dF.wb[c("FTD","gene carriers"),])
    
    if(exclNeg){
      return(list(Demographic=demog,
                  DOF=dF.cst$parameter,
                  pVal=dF.cst$p.value,
                  Fval=NA,
                  ChiSq=dF.cst$statistic,
                  geneCarriers=geneCarriers,
                  FTD=FTD,
#                   negVsFTD.stat=NA,
#                   negVsFTD.p=NA,
#                   negVsCarrier.stat=NA,
#                   negVsCarrier.p=NA,
                  FTDVsCarrier.stat=FTDVsCarrier[["statistic"]],
                  FTDVsCarrier.p=FTDVsCarrier[["p.value"]]
                  )
             )
    } else {
      return(list(Demographic=demog,
                  DOF=dF.cst$parameter,
                  pVal=dF.cst$p.value,
                  Fval=NA,
                  ChiSq=dF.cst$statistic,
                  geneNeg=geneNeg,
                  geneCarriers=geneCarriers,
                  FTD=FTD,
                  negVsFTD.stat=negVsFTD[["statistic"]],
                  negVsFTD.p=negVsFTD[["p.value"]],
                  negVsCarrier.stat=negVsCarrier[["statistic"]],
                  negVsCarrier.p=negVsCarrier[["p.value"]],
                  FTDVsCarrier.stat=FTDVsCarrier[["statistic"]],
                  FTDVsCarrier.p=FTDVsCarrier[["p.value"]]
                  )
             )
    }
  }
}

scanSummaries <- function(scanField, scanTypeList, sp=NA, exclNeg=FALSE){
  # import graph data
  dF <- importGraphData("geNorm", FALSE)
  if(exclNeg){
    dF <- dF[dF$GS!="gene negative",]
  }
  
  if(is.data.frame(sp)){dF <- applySP(dF, sp, spVal = 10)}
  
  type=scanTypeList[scanField]
  
  # merge with scanFieldraphic data
  dF.scanField <- read.table("../genfi_Subjects_sjones_1_22_2015_17_47_47_restructure_summary.csv", sep="\t", header = TRUE)
  names(dF.scanField)[4] <- "wbic"
  dF.scanField <- dF.scanField[dF.scanField$scan_duration>100,] # exclude lines with subjects with only short scans
  
  dF <- merge(dF, dF.scanField[,-which(names(dF.scanField)=="GS")],
              by = c("wbic"),
              all.x=TRUE,
              all.y=FALSE)
  names(dF)[which(names(dF)==scanField)] <- "scanField"
  
  scanFieldSummary=c()
  if(type=="continuous"){
    scanFieldSummary <- fivenum(dF$scanField, na.rm = TRUE)
    names(scanFieldSummary) <- c("min", "low.hinge", "median", "high.hinge", "max")
    scanFieldSummary$mean = mean(dF$scanField, na.rm=TRUE)
    for(i in c("categories", "n")){
      scanFieldSummary[i] = NA
    }
  } else if(type=="category"){
    for(i in c("min", "low.hinge", "median", "high.hinge", "max", "mean")){
      scanFieldSummary[i] = NA
    }
    tSum <- table(dF$scanField)
    names(tSum) <- sapply(names(tSum), function(x) strsplit(x, " ")[[1]][[1]])
    
    scanFieldSummary$categories = paste(names(tSum), collapse = "/")
    scanFieldSummary$n = paste(tSum, collapse="/")
  }
  return(scanFieldSummary)
}

clinicalScores <- function(metric,
                           metricName,
                           cs,  # clinical score
                           csName,  # clinical score name
                           cols,
                           sp, weighted=TRUE,
                           outDir="wholeBrainResults",
                           edgePC=3,
                           h=15,w=15,s=4,tsz=12,ps=4,  # tsz = text size
                           exclNeg=FALSE, # TRUE to exclude gene negative subjects
                           family=FALSE,
                           lobe=NA, # define the lobe of the brain to examine
                           hubT=NA, # hub threshold
                           csFile="clinicalScores.csv"
                           ){
  
  # create output directory
  dir.create(outDir, showWarnings = FALSE)
  
  # define log output file
  if(weighted){
    outFile = paste(outDir,
                    paste(metric,cs,"wt",sep="_"),
                    sep="/")
    
  } else {
    outFile = paste(outDir,
                    paste(metric,cs,sep="_"),
                    sep="/")
    
  }
  
  logFile = paste(outFile, "logFile.tex", sep="_")
  
  # create log file
  initiateLog(logFile, metricName)
  
  # import graph data
  dF <- importGraphData(metric, weighted, edgePC, lobe=lobe, hubT=hubT)
  
  # apply spike percentage threshold
  dF <- applySP(dF, sp)
  
  if(exclNeg){
    dF <- dF[dF$GS!="gene negative",]
    dF$GS <- factor(as.character(dF$GS), levels=c("gene carriers", "FTD"))
  }
  
  colList = unlist(cols[levels(dF$GS)])
  
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
  
  t1 <- xtable(ptSum,
               caption="Subjects included in the analysis",
               digits = c(0,0,0,0,0,0),
               display = c("s", "s", "d", "d", "d","d"))
  
  # stack data and take the mean if it is a node-wise measures
  dF.wb <- stackIt(dF, metric)
  
  # add lobe information if necessary
  if(!is.na(lobe)){
    dF.wb <- data.frame(dF.wb, lobe=lobe)
    lobeTag=lobe
  } else {
    lobeTag=""
  }

  if(!exclNeg){
    dF.wb.summary = ddply(dF.wb, .(), summarise,
                          "gene negative" = paste(sapply(mean(values[GS=="gene negative"], na.rm = TRUE), fn, a=2,b=3),
                                                  paste("(", sapply(sd(values[GS=="gene negative"], na.rm = TRUE), fn, a=2,b=3),   ")", sep="")),
                          "gene carriers" = paste(sapply(mean(values[GS=="gene carriers"], na.rm = TRUE), fn, a=2,b=3),
                                                  paste("(", sapply(sd(values[GS=="gene carriers"], na.rm = TRUE),fn, a=2,b=3),   ")", sep="")),
                          "FTD"      = paste(sapply(mean(values[GS=="FTD"], na.rm = TRUE),fn, a=2,b=3),
                                             paste("(", sapply(sd(values[GS=="FTD"], na.rm = TRUE),fn, a=2,b=3),   ")", sep=""))
    )
    
  } else {
    dF.wb.summary = ddply(dF.wb, .(), summarise,
                          "gene carriers" = paste(sapply(mean(values[GS=="gene carriers"], na.rm = TRUE), fn, a=2,b=3),
                                                  paste("(", sapply(sd(values[GS=="gene carriers"], na.rm = TRUE),fn, a=2,b=3),   ")", sep="")),
                          "FTD"      = paste(sapply(mean(values[GS=="FTD"], na.rm = TRUE),fn, a=2,b=3),
                                             paste("(", sapply(sd(values[GS=="FTD"], na.rm = TRUE),fn, a=2,b=3),   ")", sep=""))
    )
    
  }
  
  #   print(xtable(dF.wb.summary[,-1],
  #                caption=paste("Mean and standard deviations for",metricName,"values in individuals")),
  #         include.rownames=FALSE,
  #         file=logFile,
  #         append=TRUE)
  
  if(exclNeg){
    dF.wb <- dF.wb[dF.wb$GS!="gene negative",]
    dF.wb$GS <- as.factor(as.character(dF.wb$GS))
  }
  
  # import clinical score data (includes TIV)
  csdF = read.table("/home/tim/GENFI/clinicalScores.csv", sep=",", header = TRUE)
  names(csdF)[which(names(csdF)=="BLINDID")] <- "wbic"
  csdF <- csdF[,c("wbic", cs)]
  names(csdF)[which(names(csdF)==cs)] = "score"
  
  # merge network and imaging data
  dF.wb <- merge(dF.wb, csdF,
              by = c("wbic"),
              all.x=TRUE,
              all.y=FALSE)
  
  # report the number of subjects with/without data
  ptSum = ddply(dF.wb, .(gene), summarise,
                "gene negative" = length(GS[GS=="gene negative"]),
                "gene carriers" = length(GS[GS=="gene carriers"]),
                "FTD"      = length(GS[GS=="FTD"])
  )
  
  dF.wb <- dF.wb[complete.cases(dF.wb),]
  ptSum.included = ddply(dF.wb, .(gene), summarise,
                         "gene negative" = length(GS[GS=="gene negative"]),
                         "gene carriers" = length(GS[GS=="gene carriers"]),
                         "FTD"      = length(GS[GS=="FTD"])
  )
  
  ptSum.excluded = cbind(gene=ptSum[,1], ptSum[,-1] - ptSum.included[,-1])
  

  # perform mixed effects analysis
  # ANOVA of the differences
  # set contrasts
  mat = rbind(c(-1,1))
  cMat = ginv(mat)
  
  if(family){
    mod <- lmer(score ~ values*GS + Age + (1 | gene) + (1 | site) + (1 | Family),
                data=dF.wb,
                REML=FALSE,
                contrasts = list(GS=cMat))
  } else {
    mod <- lmer(score ~ values*GS + Age + (1 | gene) + (1 | site),
                data=dF.wb,
                REML=FALSE,
                contrasts = list(GS=cMat))
  }
  
  # print random effects of the model
  mod.coef <- ranef(mod)
  
  t3 <- xtable(mod.coef$site,
               digits=c(0,2),
               display=c("s", "fg"),
               caption = paste("Linear mixed effects model, site coefficients",lobeTag, csName))
  
  t4 <- xtable(mod.coef$gene,
               digits=c(0,2),
               display=c("s", "fg"),
               caption = paste("Linear mixed effects model, gene coefficients",lobeTag, csName))
  
  if(family){
    t5 <- xtable(mod.coef$Family,
                 digits=c(0,2),
                 display=c("s", "fg"),
                 caption = paste("Linear mixed effects model, family coefficients",lobeTag, csName))
  } else {
    t5 = NA
  }
  
  print(t3, file=logFile, append=TRUE)
  print(t4, file=logFile, append=TRUE)
  if(family){print(t3, file=logFile, append=TRUE)}
  
  # plot random effects
  dF.plot.site <- data.frame(site = row.names(mod.coef$site),
                             effect = mod.coef$site[[1]])
  p1 <- ggplot(dF.plot.site, aes_string(x="site",y="effect"))
  p1 <- p1 + geom_bar(stat="identity") + theme_bw() + ggtitle(paste("Effect size of scan site", csName))
  
  dF.plot.gene <- data.frame(gene = row.names(mod.coef$gene),
                             effect = mod.coef$gene[[1]])
  p2 <- ggplot(dF.plot.gene, aes_string(x="gene",y="effect"))
  p2 <- p2 + geom_bar(stat="identity") + theme_bw() + ggtitle(paste("Effect size of gene", csName))
  
  if(family){
    dF.plot.Family <- data.frame(Family = row.names(mod.coef$Family),
                                 effect = mod.coef$Family[[1]])
    p3 <- ggplot(dF.plot.Family, aes_string(x="Family",y="effect"))
    p3 <- p3 + geom_bar(stat="identity") + theme_bw() + ggtitle(paste("Effect size of family", csName))
    p3 <- p3 + theme(axis.text.x=element_blank())
  } else {
    p3 = NA
  }
  
  # print variances
  vc <- VarCorr(mod)
  
  t6 <- xtable(data.frame(vc),
               display=c("s","s","s","s","g","fg"),
               digits=c(0,0,0,0,2,2),
               caption=paste("Variance of random effects", csName))
  
  print(t6,
        file=logFile,
        append=TRUE,
        include.rownames=FALSE)
  
  t7 <- xtable(summary(mod)[[10]],
               digits=c(0,2,2,1,2,2),
               display=c("s","fg","f","f","f","g"),
               caption=paste("Satterthwaite estimates of pvalues of linear mixed effects model for", csName, metricName,lobeTag))
  
  print(t7,
        file=logFile,
        append=TRUE)
  
  # plot
  p4 <- ggplot(data=dF.wb, aes_string(x="values", y="score", colour="GS"))
  p4 <- p4 + geom_smooth(method="lm")
  
  p4 <- p4 + theme_bw()
  p4 <- p4 + labs(x=csName, y=metricName)
  p4 <- p4 + theme(text=element_text(size=tsz))
  p4 <- p4 + scale_fill_manual(name="Group",values=as.vector(colList))
  p4 <- p4 + theme(legend.position="none")
  p4 <- p4 + labs(title=paste(csName,lobeTag))
  
  if(exclNeg){
    outw = 2.5
  } else {
    outw = 3
  }
  
  plotFile = paste(paste(outFile,"allgroups",sep="_"), "png", sep=".")
  ggsave(plotFile,
         scale=s,
         dpi=600,
         height=h, width=w,
         units="mm")
  
  # return tables and plots
  return(list(t1 = ptSum.included,
              t2 = ptSum.excluded,
              t3 = t3,
              t4 = t4,
              t5 = t5,
              t6 = t6,
              t7 = t7,
              p1 = p1,
              p2 = p2,
              p3 = p3,
              p4 = p4))
}
  