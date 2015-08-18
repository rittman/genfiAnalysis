# This script will perform t-test analysis between gene negative and gene carrier subjects for node-wise metrics
# The results should be FDR corrected as well

genfiDir = "/home/tim/GENFI/GENFI_camgrid_20150525/"
setwd(genfiDir)

#### helper functions ####
nodeTtest <- function(vals, g1len, g2len){
  res <- t.test(vals[c(1:g1len)], vals[c((g1len+1):(g1len+g2len))])

  return(list(t=res$statistic,
              p=res$p.value))
}


#### main functions ####
nodeComparison <- function(dF){
  # function to compare the values of gene carrier subjects and controls
  # control subjects
  conts = data.frame("GS"="Gene negative",
                     dF[dF$GS==0,c(-1,-2,-3,-4)])
  conts.len = length(conts[,1])
  
  # gene carrier values
  gcs = data.frame("GS"="Gene positive",
                   dF[dF$GS==1,c(-1,-2,-3,-4)])
  gcs.len = length(gcs[,1])
  
  dF <- rbind(conts,gcs)
  View(dF)
  
  dF <- dF[,apply(dF, 2, function(x) !any(is.na(x)))]
  
  resVals = apply(dF[,-1], 2, nodeTtest, g1len=conts.len, g2len=gcs.len)
  resVals <- data.frame(matrix(unlist(resVals), byrow=TRUE, nrow=length(resVals)))
  names(resVals) = c("t", "p")
  resVals <- data.frame(nodes=names(dF[,-1]),
                        resVals,
                        p.adjust(resVals$p, method="BY"))
  
  View(resVals)
  return(resVals)
}

#### run functions ####
# import data
gmList = list("d2_ccNorm_local",
              "d2_elnNorm_local"
              )

# select only a single percentage edge value
edgePC = 3

for(gm in gmList){
  gm.dF <- read.table(gm, header=TRUE)
  gm.dF <- gm.dF[gm.dF$"edgePC"==edgePC,]
#   print(apply(gm.dF, 2, function(x) !any(is.na(x))))
#   gm.dF <- gm.dF[,apply(gm.dF, 2, function(x) !any(is.na(x)))]
  resVals = nodeComparison(gm.dF[,-1])
}


