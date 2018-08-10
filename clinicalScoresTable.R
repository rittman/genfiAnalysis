# Script to tabulate clinical scores as per reviewer request
# Timothy Rittman 10/8/2018

library(plyr)
library(xtable)

# import clinical scores
csdF = read.table("/home/tim/GENFI/clinicalScores.csv", sep=",", header = TRUE)
names(csdF)[which(names(csdF)=="BLINDID")] <- "wbic"

# factorise the gene status
csdF$GS <- factor(csdF$GS, labels = c("Negative", "Carrier", "FTD"))

# summarise clinical scores
cs.summary <- ddply(csdF, .(GS), summarise,
                    mmse.mean=mean(N1, na.rm = TRUE),
                    mmse.sd=sd(N1, na.rm = TRUE),
                    LogMemoryImmediate.mean=mean(N2, na.rm = TRUE),
                    LogMemoryImmediate.sd=sd(N2, na.rm = TRUE),
                    LogMemoryDelayed.mean=mean(N3, na.rm = TRUE),
                    LogMemoryDelayed.sd=sd(N3, na.rm = TRUE),
                    DigitSpanForward.mean=mean(N4, na.rm = TRUE),
                    DigitSpanForward.sd=sd(N4, na.rm = TRUE),
                    DigitSpanBackward.mean=mean(N5, na.rm = TRUE),
                    DigitSpanBackward.sd=sd(N5, na.rm = TRUE),
                    TrailsA.mean=mean(N6, na.rm = TRUE),
                    TrailsA.sd=sd(N6, na.rm = TRUE),
                    TrailsB.mean=mean(N7, na.rm = TRUE),
                    TrailsB.sd=sd(N7, na.rm = TRUE),
                    DigitSymbol.mean=mean(N8, na.rm = TRUE),
                    DigitSymbol.sd=sd(N8, na.rm = TRUE),
                    BostonNaming.mean=mean(N9, na.rm = TRUE),
                    BostonNaming.sd=sd(N9, na.rm = TRUE),
                    VerbalFluencyCategory.mean=mean(N10, na.rm = TRUE),
                    VerbalFluencyCategory.sd=sd(N10, na.rm = TRUE),
                    VerbalFluencyLetter.mean=mean(N11, na.rm = TRUE),
                    VerbalFluencyLetter.sd=sd(N11, na.rm = TRUE),
                    BlockDesign.mean=mean(N12, na.rm = TRUE),
                    BlockDesign.sd=sd(N12, na.rm = TRUE))

cs.summary.tab <- ddply(csdF, .(GS), summarise,
                    mmse=paste(round(mean(N1, na.rm = TRUE), digits = 1),
                               " (",round(sd(N1, na.rm = TRUE), digits=1),")",
                               sep=""),
                    LogMemoryImmediate=paste(round(mean(N2, na.rm = TRUE), digits = 2),
                                                  " (",round(sd(N2, na.rm = TRUE), digits=2),")",
                                                  sep=""),
                    
                    LogMemoryDelayed=paste(round(mean(N3, na.rm = TRUE), digits = 2),
                                                " (",round(sd(N3, na.rm = TRUE), digits=2),")",
                                                sep=""),
                    DigitSpanForward=paste(round(mean(N4, na.rm = TRUE), digits = 2),
                                                " (",round(sd(N4, na.rm = TRUE), digits=2),")",
                                                sep=""),
                    DigitSpanBackward=paste(round(mean(N5, na.rm = TRUE), digits = 2),
                                                 " (",round(sd(N5, na.rm = TRUE), digits=2),")",
                                                 sep=""),
                    TrailsA=paste(round(mean(N6, na.rm = TRUE), digits = 2),
                                       " (",round(sd(N6, na.rm = TRUE), digits=2),")",
                                       sep=""),
                    TrailsB=paste(round(mean(N7, na.rm = TRUE), digits = 2),
                                       " (",round(sd(N7, na.rm = TRUE), digits=2),")",
                                       sep=""),
                    DigitSymbol=paste(round(mean(N8, na.rm = TRUE), digits = 2),
                                      " (",round(sd(N8, na.rm = TRUE), digits = 2),")",
                                      sep=""),
                    BostonNaming=paste(round(mean(N9, na.rm = TRUE), digits = 2),
                                       " (",round(sd(N9, na.rm = TRUE), digits = 2),")",
                                       sep=""),
                    VerbalFluencyCategory=paste(round(mean(N10, na.rm = TRUE), digits = 2),
                                                " (",round(sd(N10, na.rm = TRUE), digits = 2),")",
                                                sep=""),
                    VerbalFluencyLetter=paste(round(mean(N11, na.rm = TRUE), digits = 2),
                                              " (",round(sd(N11, na.rm = TRUE), digits = 2),")",
                                              sep=""),
                    BlockDesign.mean=paste(round(mean(N12, na.rm = TRUE), digits = 2),
                                           " (",round(sd(N12, na.rm = TRUE), digits = 2),")",
                                           sep="")
)


# the following is for reference only
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


# write table of clinical scores

# rotate for output
cs.summary.tab.out <- t(cs.summary.tab)

coglist =list(
  GS="",
  N1="MMSE",
  N2="Log Immediate Memory",
  N3="Log Delayed Memory",
  N4="Forward Digit Span",
  N5="Backwards Digit Span",
  N6="Trails A",
  N7="Trails B",
  N8="Digit Symbol Task",
  N9="Boston Naming Task",
  N10="Verbal Fluency (Category)",
  N11="Verbal Fluency (Letter)",
  N12="Block Design Task"
)
cs.summary.tab.out <- cbind(coglist, cs.summary.tab.out)

outtable <- xtable(cs.summary.tab.out)

print(outtable, include.rownames=FALSE, include.colnames = FALSE)

write.table(cs.summary.tab.out, "clinicalScoresTable.csv", col.names = FALSE, row.names = FALSE,
            sep=",")
