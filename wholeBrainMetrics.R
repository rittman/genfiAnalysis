# This script will examine whole brain graph metrics in the GENFI dataset by looking for differences
# between affected/carriers/non-carriers, and by examining the relationship between average age at
# onset and the graph measures.

library(plyr)
library(ggplot2)
library(xtable)
library(car)
library(xlsx)

genfiDir = "/home/tim/GENFI/GENFI_camgrid_20150525/"
setwd(genfiDir)

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

wholeBrainAnalysis <- function(metric, metricName, sp, weighted=TRUE, outDir="wholeBrainResults"){
  if(weighted){
    metric = paste(metric,"wt",sep="_")
  }
  
  # create output directory
  dir.create(outDir, showWarnings = FALSE)
  
  # define input file
  inFile = paste("d2",metric,"local",sep="_")
  
  # define log output file
  outFile = paste(outDir,
                  paste(metric,"logFile.tex",sep="_"),
                  sep="/")
  
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
  
  # import data
  dF = read.table(inFile, header = TRUE, na.strings = "NA")
  
  # convert diagnostic label to a factor and rename
  dF$GS = as.factor(dF$GS)
  dF$GS = revalue(dF$GS, c("0" = "gene negative", "1"="gene positive", "2"="affected"))
  
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
  
  # set spike percentage threshold
  spVal = 10  # would be better to define this using a Bayesian model averaging method

  # import file with spike percentage data
  sp.sub <- data.frame(wbic=sp$id, spMean=sp$mean)
  dF <- merge(dF, sp.sub, by="wbic")
  
  dF <- dF[dF$spMean < spVal,]
  
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
  
  # Sort out stacking the nodes up if nodewise measure
  nodeNames = names(dF)[sapply(names(dF), function(x) grepl("X",x))]
  
  if(length(nodeNames)>0){
    # Take the mean values of
    dF.stacked = stack(dF[,nodeNames])
    dF.stacked = data.frame(dF.stacked, GS=dF$GS, gene=dF$gene, wbic=dF$wbic)
    
    dF.wb = ddply(dF.stacked, .(gene, wbic, GS), summarise,
                  values = mean(values, na.rm = TRUE)
    )
  } else {
    dF.wb <- dF
    names(dF.wb)[names(dF.wb)==metric] <- "values"
  }
  
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
  
  write("\\end{document}", file=outFile, append=TRUE)
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

lapply(names(metrics), function(x) wholeBrainAnalysis(x, metrics[[x]], sp))

# do metrics for unweighted graphs
metrics = list("degree"="connection strength",
               "degreeNorm"="connection strength",
               "ge"="global efficiency",
               "geNorm"="global efficiency",
               "le"="local efficiency",
               "leNorm"="local efficiency",
               "pl"="path length",
               "plNorm"="path length",
               "eigCentNP"="eigen centrality",
               "eigCentNorm"="eigen centrality",
               "betCent"="betweenness centrality",
               "betCentNorm"="betweenness centrality",
               "closeCent"="closeness centrality",
               "closeCentNorm"="closeness centrality")

lapply(names(metrics), function(x) wholeBrainAnalysis(x, metrics[[x]], sp, weighted=FALSE))




