#!/share/apps/R-3.6.1/bin/Rscript

TestNames=c("HL.rec.20230823","HT.rec.20230823","T2D.rec.20230823")

library(ggplot2)


qqSLPs<-function(SLPs) {
        rankSLP <- rank(SLPs, na.last=TRUE, ties.method="first")
        nelem=length(SLPs)
        midrank <- (nelem+1)/2
        rMinusMr <- rankSLP-midrank
        absDiff <- abs(rMinusMr/midrank)
        pVal <- 1-absDiff
        logP <- log10(pVal)
        eSLPs <- sign(rankSLP - midrank) * -logP
        return(eSLPs)
}

for (TestName in TestNames) {

SummResFile=sprintf("AllRec.%s.summ.txt",TestName)
SummRes=data.frame(read.table(SummResFile,header=TRUE,sep="\t",stringsAsFactors=FALSE))

options(bitmapType='cairo')
for (MinCarrier in c(1,20,50)) {
for (s in 1:3) {
Narrowed=SummRes[SummRes$NumCarriers>=MinCarrier,]
SLP=na.omit(Narrowed[,2+s])
top=max(SLP)
bottom=min(SLP)
top=5
bottom=-5
SLP=SLP[order(-SLP)]
eSLP=qqSLPs(SLP)
filename=sprintf("QQ.%s.%s.%02d.png",TestName,colnames(SummRes)[2+s],MinCarrier)
ppi=600
png(filename,width=6*ppi, height=6*ppi, res=ppi)
toPlot=data.frame(matrix(ncol=2,nrow=length(SLP)))
toPlot[,1]=eSLP
toPlot[,2]=SLP
colnames(toPlot)=c("eSLP","SLP")
myplot=ggplot(toPlot,aes_q(x=as.name("eSLP"),y=as.name("SLP")))+geom_point(size=1)+ theme_bw() +
                geom_hline(yintercept=0,size=1.0) +
                geom_vline(xintercept=0,size=1.0) +
                theme(panel.grid.major=element_line(colour = "black",size=0.25)) +
                scale_x_continuous(breaks = seq(2*floor(bottom/2),2*ceiling(top/2),by =2),minor_breaks=NULL,limits=c(2*floor(bottom/2),2*ceiling(top/2))) +
                scale_y_continuous(breaks = seq(2*floor(bottom/2),2*ceiling(top/2),by =2),minor_breaks=NULL,limits=c(2*floor(bottom/2),2*ceiling(top/2)))
print(myplot)
dev.off()
pos=toPlot[which(SLP>0),]
pos=pos[101:nrow(pos),] # remove top 100 genes
l=lm(pos[,1] ~ pos[,2])
print(l)
neg=toPlot[which(SLP<0),]
neg=neg[1:(nrow(pos)-100),]
l=lm(neg[,1] ~ neg[,2])
print(sprintf("%s %s %02d",TestName,colnames(SummRes)[2+s],MinCarrier))
print(l)
}
}
}