# This script will examine whole brain graph metrics in the GENFI dataset by looking for differences
# between affected/carriers/non-carriers, and by examining the relationship between average age at
# onset and the graph measures.

library(plyr)
library(ggplot2)
library(xtable)
library(car)
library(xlsx)
library(lme4)

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

#### main functions ####
wholeBrainAnalysis <- function(metric,
                               metricName,
                               sp, weighted=TRUE,
                               outDir="wholeBrainResults",
                               edgePC=3){
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
  
  # Summarise the patient data
  ptSum = ddply(dF, .(gene), summarise,
                "gene negative" = length(GS[GS=="gene negative"]),
                "gene positive" = length(GS[GS=="gene positive"]),
                "affected"      = length(GS[GS=="affected"])
  )
  
  print(xtable(ptSum,
               caption="Subjects included in the analysis"),
        include.rownames=FALSE,
        file=outFile,
        append=TRUE)
  
  # apply spike percentage threshold
  dF <- applySP(dF, sp)
  
  # Summarise the patient data
  ptSum = ddply(dF, .(gene), summarise,
                "gene negative" = length(GS[GS=="gene negative"]),
                "gene positive" = length(GS[GS=="gene positive"]),
                "affected"      = length(GS[GS=="affected"])
  )
  
  print(xtable(ptSum,
               caption="Subjects included in the analysis"),
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
  
  #### post hoc t-tests ####
  # group pairwise
  pwGS <- combn(levels(dF.wb$GS),2)
  pwTtests <- apply(pwGS, 2, function(x) dotTests(x,dF.wb))
  pwTtests <- data.frame(matrix(unlist(pwTtests), nrow = length(pwTtests), byrow = TRUE))
  names(pwTtests) <- c("comparison", "t", "p")
  print(xtable(pwTtests,
               caption="Pairwise post-hoc t-test results"),
        include.rownames=FALSE,
        file=outFile,
        append=TRUE)
  
  # t-test for affected vs non-affected (gene negative and gene positive)
  anTtest <- t.test(dF.wb[dF.wb$GS=="affected","values"], dF.wb[dF.wb$GS!="affected","values"])
  print(xtable(data.frame("t"=fn(anTtest[["statistic"]]),
                          "p"=fn(anTtest[["p.value"]])),
               caption="t-test between affected subjects and non-affected (both gene negative and gene positive unaffected",
               digits=c(0,2,2),
               display=c("s","fg","fg")
               ),
        include.rownames=FALSE,
        file=outFile,
        append=TRUE)
        
  
  # now do a plot
  p <- ggplot(dF.wb, aes_string(x="GS", y="values", fill="GS"))
  p <- p + geom_boxplot()
  p <- p + theme_bw()
  p <- p + labs(y=metricName) + theme(axis.title.x=element_blank())
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
  
  # # plot this
  # err <- function(x) qnorm(0.975) * sd(x, na.rm = TRUE)/sqrt(length(x)) # function for standard error
  # 
  # dF.stacked.plot = ddply(dF.stacked, .(gene, GS), summarise, valuesMean=mean(values, na.rm = TRUE),
  #                         upper = valuesMean+err(values),
  #                         lower = valuesMean-err(values)
  #                         )
  
  # dF.stacked.plot <- data.frame(dF.stacked.plot, unique=seq(length(dF.stacked.plot[,1])))
  
  # p <- ggplot(dF.stacked.plot, aes_string(x="unique", y="valuesMean", middle="valuesMean", group="unique", fill="gene", upper="upper", lower="lower"))
  p <- ggplot(dF.wb, aes_string(x="gene", y="values", fill="GS"))
  p <- p + geom_boxplot()
  p <- p + theme_bw()
  p <- p + labs(y=metricName) + theme(axis.title.x=element_blank())
  ggsave(paste(outDir,
               paste(metric,"byGene.png",sep="_"),
               sep="/"))
  
  endLog(outFile)
}


graphTimeComparison <- function(metric,
                               metricName,
                               sp, # spike percentage
                               startvec=NULL, # starting vectors for equation with quadratic term
                               startvecCub=NULL, # starting vectors for equation with cubic term
                               weighted=TRUE, # is this a weighted metric?
                               outDir="wholeBrainVsAOOResults",
                               edgePC=3){
  print(metricName)
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
  dF.plot <- dF
  
  lm.all <- lm(values ~ aoo*GS, data=dF)

  print(xtable(summary(lm.all),
               caption=paste("summary of linear model for the interaction between estimated age of onset (aoo) and gene status (GS) on", metricName),
               digits=c(0,2,2,2,2),
               display=c("s","fg","fg","fg","g")
  ),
  include.rownames=TRUE,
  file=outFile,
  append=TRUE)

  # now run linear mixed effects model model
  # run model
  mod <- lmer(values ~ aoo + GS + (1 | Family) + (1 | site),
              data=dF,
              REML=FALSE)
  print(xtable(summary(mod)[["coefficients"]],
               caption="Linear mixed effects model, fixed effects",
               digits=c(0,2,2,2),
               display=c("s","fg","fg","fg")
  ),
  include.rownames=TRUE,
  file=outFile,
  append=TRUE)

  print(xtable(data.frame(StdDev = c(attributes(VarCorr(mod)[[1]])[["stddev"]],
                                     attributes(VarCorr(mod)[[2]])[["stddev"]],
                                     attributes(VarCorr(mod))[["sc"]]),
                          row.names = c("Family", "Site", "Residual")),                          
               caption="Linear mixed effects model, random effects",
               digits=c(0,2),
               display=c("s","fg")
  ),
  include.rownames=TRUE,
  file=outFile,
  append=TRUE)

  # null model
  nulMod <- lmer(values ~ aoo + (1 | Family) + (1 | site),
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
  
  p <- ggplot(dF.plot, aes_string(x="aoo", y="values", colour="GS", group="GS"))
  p <- p + geom_point()
  p <- p + geom_smooth(method="lm")
  p <- p + theme_bw()
  p <- p + labs(y=metricName, x="Estimated age of onset") + theme(axis.title.x=element_blank())

  ggsave(paste(outDir,
               paste(metric,"png",sep="."),
               sep="/"))
  
  # now run non-linear mixed effects model with a quadratic term
  # firstly, plot the values to be able to estimate starting parameters
  if(is.null((startvec))){

    affirm = "n"
    while(affirm=="n"){
      print(p)
      
#       Asym = as.numeric(readline(prompt="Asymptote of the y-axis: "))
#       R0   = as.numeric(readline(prompt="the value of y when x=0: "))
#       lrc  = as.numeric(readline(prompt="natural log of the rate of decline (guess -1 if not sure): "))
      aq = as.numeric(readline(prompt = "Constant term: "))
      bq = as.numeric(readline(prompt = "x multiplier:" ))
      cq = as.numeric(readline(prompt = "quadratic multiplier: "))
      
      # dF.temp <- data.frame(x=dF$aoo, y=SSasymp(dF$aoo, Asym, R0, lrc), GS="Start estimates")
      dF.temp <- data.frame(x=dF$aoo, y=SSquadFun(dF$aoo, aq=aq, bq=bq, cq=cq), GS="Start estimates")
      p.temp <- p + geom_point(data=dF.temp, aes_string(x="x", y="y"))
      print(p.temp)
      affirm = readline(prompt = "Are you happy with this?(y/n)")
    }
    # startvec = c(Asym=Asym, R0=R0, lrc=lrc)
    startvec = c(aq=aq, bq=bq, cq=cq)
  }
  print(paste("Starting values are as follows.",
              "Constant:", startvec["aq"],
              "Multiplier of x term:", startvec["bq"],
              "Multiplier of quadratic term:", startvec["cq"])
        )
  
  # run non-linear quadratic model
  dF <- dF[dF$GS!="gene negative",]
  # dF <- cbind(dF, data.frame(aooNeg=dF$aoo*-1))  # make age of onset negative so that an asymptotic fit can be achieved.

  # The following line calculates the non-linear model
  # the random elements take in to account the difference in graph measures depending on the scanner site (Asym|site) and the age at onset between genes (aooNeg|gene)
  # at present the function is an asymptote.

  nlmq <- nlmer(values ~ SSquadFun(aoo, aq, bq, cq) ~ (aoo|gene) + (aq|site), data=dF, start = startvec)
  
  print(summary(nlmq))

  # plot estimated values
  x = seq(-45,25,by=1)
  dF.nlmq <- data.frame(x=x, y=SSquadFun(x, nlmq@beta[[1]], nlmq@beta[[2]], nlmq@beta[[3]]), GS="Estimates")
  pq <- ggplot(data=dF.nlmq, aes_string(x="x", y="y"))
  pq <- pq + geom_line()
  pq <- pq + theme_bw()
  pq <- pq + labs(y=metricName, x="Estimated time from onset")# + theme(axis.title.x=element_blank())
  
  print(paste("Saving",paste(outDir,
                       paste(paste(metric, "nonLinearEstimatesQuadratic",sep="_"),"png",sep="."),
                       sep="/")))
  ggsave(paste(outDir,
               paste(paste(metric, "nonLinearEstimatesQuadratic",sep="_"),"png",sep="."),
               sep="/"))
  
  # now run non-linear mixed effects model with a cubic term
  # firstly, plot the values to be able to estimate starting parameters
  ac = nlmq@beta[[1]]
  bc = nlmq@beta[[2]]
  cc = nlmq@beta[[3]]
  
  if(is.null((startvecCub))){
    
    affirm = "n"
    while(affirm=="n"){
      print(p)
      dc = as.numeric(readline(prompt = "cubic multiplier: "))
      
      # dF.temp <- data.frame(x=dF$aoo, y=SSasymp(dF$aoo, Asym, R0, lrc), GS="Start estimates")
      dF.temp <- data.frame(x=dF$aoo, y=SScubicFun(dF$aoo, ac=ac, bc=bc, cc=cc, dc=dc), GS="Start estimates")
      p.temp <- p + geom_point(data=dF.temp, aes_string(x="x", y="y"))
      print(p.temp)
      affirm = readline(prompt = "Are you happy with this?(y/n)")
    }
    # startvec = c(Asym=Asym, R0=R0, lrc=lrc)
    startvec = c(ac=ac, bc=bc, cc=cc, dc=dc)
  }
  print(paste("Starting values are as follows.",
              "Constant:", startvec["ac"],
              "Multiplier of x term:", startvec["bc"],
              "Multiplier of quadratic term:", startvec["cc"],
              "Multiplier of cubic term:", startvec["dc"])
  )
  
  # run non-linear quadratic model
  # The following line calculates the non-linear model
  # the random elements take in to account the difference in graph measures depending on the scanner site (aq|site) and the age at onset between genes (aooNeg|gene)
  # this equation includes a cubic term
  nlmc <- nlmer(values ~ SScubicFun(aoo, ac, bc, cc, dc) ~ (aoo|gene) + (ac|site), data=dF, start = startvec)

  # print a summary of the model
  print(summary(nlmc))
  
  # plot estimated values
  x = seq(-45,25,by=1)
  dF.nlmc <- data.frame(x=x, y=SScubicFun(x, nlmc@beta[[1]], nlmc@beta[[2]], nlmc@beta[[3]], nlmc@beta[[4]]), GS="Estimates")
  pc <- ggplot(data=dF.nlmc, aes_string(x="x", y="y"))
  pc <- pc + geom_line()
  pc <- pc + theme_bw()
  pc <- pc + labs(y=metricName, x="Estimated time from onset")# + theme(axis.title.x=element_blank())
  
  print(paste("Saving",paste(outDir,
                             paste(paste(metric, "nonLinearEstimatesCubic",sep="_"),"png",sep="."),
                             sep="/")))
  ggsave(paste(outDir,
               paste(paste(metric, "nonLinearEstimatesCubic",sep="_"),"png",sep="."),
               sep="/"))
  
  # now compare the quadratic and cubic equations
  qvsc.anova <- anova(nlmq, nlmc)
  print(xtable(qvsc.anova,
               caption = "ANOVA between models with quadratic and cubic terms to explain rate of change in graph measure with time",
               digits = c(0,0,2,2,2,2,2,2,2),
               display = c("s","d","fg","fg","fg","fg","d","fg","fg")),
        include.rownames=TRUE,
        file=outFile,
        append=TRUE)
  
  # finalise log file and return dataframe
  endLog(outFile)
  return(dF)
}


# import spike percentage data
# sp <- read.xlsx("../all_sm_thld10_SP.xlsx", sheetIndex = 1)

metrics = list("degree"="connection strength",
               "ge"="global efficiency",
               "le"="local efficiency",
               "pl"="path length",
               "eigCentNP"="eigen centrality",
               "betCent"="betweenness centrality",
               "closeCent"="closeness centrality")

# lapply(names(metrics), function(x) wholeBrainAnalysis(x, metrics[[x]], sp))
# lapply(names(metrics), function(x) graphTimeComparison(x, metrics[[x]], sp))

# do metrics for unweighted graphs
metrics = list("degree"="connection strength",
               "degreeNorm"="connection strength (normalised)",
               "ge"="global efficiency",
               "geNorm"="global efficiency (normalised)",
               "le"="local efficiency",
               "leNorm"="local efficiency (normalised)",
               "pl"="path length",
               "plNorm"="path length (normalised)",
               "eigCentNP"="eigen centrality",
               "eigCentNorm"="eigen centrality (normalised)",
               "betCent"="betweenness centrality",
               "betCentNorm"="betweenness centrality (normalised)",
               "closeCent"="closeness centrality",
               "closeCentNorm"="closeness centrality (normalised)")

# lapply(names(metrics), function(x) wholeBrainAnalysis(x, metrics[[x]], sp, weighted=FALSE))
# lapply(names(metrics), function(x) graphTimeComparison(x, metrics[[x]], sp, weighted=FALSE))
dF <- graphTimeComparison("plNorm", "path length (normalised)", sp, weighted=FALSE)

startVecList = list(geNorm = c(aq=0.88, bq=-.0001, cq=-.00001)
                    )
