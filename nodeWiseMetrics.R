# This script will perform t-test analysis between gene negative and gene carrier subjects for node-wise metrics
# The results should be FDR corrected as well

genfiDir = "/home/tim/GENFI/GENFI_camgrid_20150525/"
setwd(genfiDir)

library(xtable)
library(ggplot2)

#### helper functions ####
nodeTtest <- function(vals, g1len, g2len){
  res <- t.test(vals[c(1:g1len)], vals[c((g1len+1):(g1len+g2len))])

  return(list(t=res$statistic,
              p=res$p.value))
}


#### main functions ####
nodeComparison <- function(dF,plotName){
  # function to compare the values of gene carrier subjects and controls
  # control subjects
  conts = data.frame("GS"="Gene negative",
                     dF[dF$GS==0,c(-1,-2,-3,-4)])
  conts.len = length(conts[,1])
  
  # gene carrier values
  gcs = data.frame("GS"="Gene positive",
                   dF[dF$GS==1,c(-1,-2,-3,-4)])
  gcs.len = length(gcs[,1])
  
  # glue together control and gene carriers
  dF <- rbind(conts,gcs)
  
  # remove any columns containing NA values
  dF <- dF[,apply(dF, 2, function(x) !any(is.na(x)))]
  
  # do a t-test for each node
  resVals = apply(dF[,-1], 2, nodeTtest, g1len=conts.len, g2len=gcs.len)
  resVals <- data.frame(matrix(unlist(resVals), byrow=TRUE, nrow=length(resVals)))  # reorganise values
  names(resVals) = c("t", "p") # update column names
  
  # add node names and adjust p-values using Benjamini & Yekutieli (2001) method
  resVals <- data.frame(nodes=names(dF[,-1]),
                        resVals,
                        p.adj=p.adjust(resVals$p, method="BY"))
  
  # select only significant nodes (using corrected p-values)
  res.plot <- resVals[resVals$p<0.001,]

  # carry on and plot if the number of rows is 0
  if(nrow(res.plot)!=0){
    print(res.plot)
    
    # get the control mean and stuff
    res.plot.cont <- data.frame(res.plot,
                           Mean = sapply(res.plot$nodes, function(x) mean(dF[,x][c(1:conts.len)])),
                           upper = sapply(res.plot$nodes, function(x) mean(dF[,x][c(1:conts.len)])) + sapply(res.plot$nodes, function(x) sd(dF[,x][c(1:conts.len)])/sqrt(length(c(1:conts.len)))),
                           lower = sapply(res.plot$nodes, function(x) mean(dF[,x][c(1:conts.len)])) - sapply(res.plot$nodes, function(x) sd(dF[,x][c(1:conts.len)])/sqrt(length(c(1:conts.len)))),
                           group = "Control"
                           )
    
    # now get the same for gene carriers
    res.plot.gcs <- data.frame(res.plot,
                               Mean = sapply(res.plot$nodes, function(x) mean(dF[,x][c((conts.len+1):gcs.len)])),
                               upper = sapply(res.plot$nodes, function(x) mean(dF[,x][c((conts.len+1):gcs.len)])) + sapply(res.plot$nodes, function(x) sd(dF[,x][c((conts.len+1):gcs.len)])/sqrt(length(c((conts.len+1):gcs.len)))),
                               lower = sapply(res.plot$nodes, function(x) mean(dF[,x][c((conts.len+1):gcs.len)])) - sapply(res.plot$nodes, function(x) sd(dF[,x][c((conts.len+1):gcs.len)])/sqrt(length(c((conts.len+1):gcs.len)))),
                               group="Gene carriers"
    )
    
    # stitch the controls to the carriers
    res.plot <- rbind(res.plot.cont, res.plot.gcs)
    
    res.plot <- res.plot[order(res.plot$Mean),] # reorder by the mean values
    res.plot$nodes = factor(res.plot$nodes, levels=unique(res.plot$nodes))
    
    dodge = position_dodge(0.2)
    p <- ggplot(res.plot, aes_string(x="nodes",y="Mean", group="group", colour="group"))
    p <- p + geom_point(position=dodge) + geom_errorbar(aes_string(ymin="lower", ymax="upper"), position=dodge)
    View(res.plot)
    ggsave(paste(plotName,"png",sep="."))
  }
  
  # return the t-test results
  View(resVals)
  return(resVals)
}

#### run functions ####
# import data
gmList = list(#"d2_ccNorm_local",
#               "d2_elnNorm_local",
#               "d2_ccWt_local",
              "d2_degree_wt_local",
              "d2_degreeWt_local",
              "d2_leNorm_local",
              "d2_le_wt_local"
              )

# select only a single percentage edge value
edgePC = 3

for(gm in gmList){
  gm.dF <- read.table(gm, header=TRUE)
  if(!grepl("wt",gm)){ # check whether the graph metric is weighted
    gm.dF <- gm.dF[gm.dF$"edgePC"==edgePC,] # if not weighted, then select only the specified edge percentage
  }
  
#   print(apply(gm.dF, 2, function(x) !any(is.na(x))))
#   gm.dF <- gm.dF[,apply(gm.dF, 2, function(x) !any(is.na(x)))]
  resVals = nodeComparison(gm.dF[,-1], plotName=gm)
}


