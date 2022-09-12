wd="C:/Users/dave/OneDrive/sharedseq/UKBB/annot"
scoFile="UKBB.HL.gerpLin.20210225.LDLR.sco"

library(ggplot2)
setwd(wd)
scores=data.frame(read.table(scoFile,header=FALSE))
colnames(scores)=c("ID","CC","Score")

scores$CC=as.factor(scores$CC)
p <- ggplot(scores, aes(x=CC, y=Score)) + geom_violin()
p


