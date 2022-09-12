wd="C:/Users/dave/OneDrive/sharedseq/UKBB/annot"
scoFile="UKBB.HL.gerp.20210225.LDLR.sco"

library(ggplot2)
library(reshape)

setwd(wd)
scores=data.frame(read.table(scoFile,header=FALSE))
colnames(scores)=c("ID","CC","Score")
drops=c("ID")
toPlot=scores
toPlot=toPlot[,!(names(toPlot) %in% drops)]
dfm=melt(toPlot,id.vars="Score")
ggplot(dfm, aes(x=Score,y=value, color=variable)) + geom_smooth() 
