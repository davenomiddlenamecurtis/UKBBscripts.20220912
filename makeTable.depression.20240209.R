#!/share/apps/R-3.6.1/bin/Rscript

# script to take make table of MaxMLPs and MLPs from 270K

BestModel="UKBB.depression.forAnnot.20240208"
TopTestsFile="depression.topTests.240208.txt"
NewModel="UKBB.depression.best.annot.20240208"
ResFile="depression.allResults.240208.txt"

wd="/home/rejudcu/UKBB/depression.20240208"
setwd(wd)

TopTests=data.frame(read.table(TopTestsFile,header=TRUE,sep="\t",stringsAsFactors=FALSE))
for (r in 1:nrow(TopTests)) {
	gene=TopTests[r,1]
	NewRes=data.frame(read.table(sprintf("summ.%s.%s.txt",NewModel,gene),header=TRUE,sep="\t",stringsAsFactors=FALSE))
	TopTests[r,4]=NewRes[1,2]
}
colnames(TopTests)[4]="MLP270K"
write.table(TopTests,ResFile,col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

library(ggplot2)
options(bitmapType='cairo')

qqMLPs<-function(MLPs) {
	rankMLP <- rank(MLPs, na.last=TRUE, ties.method="first")
	nelem=length(MLPs)
	pVal=1-(rankMLP-1)/nelem
	eMLPs=-log10(pVal)
	return(eMLPs)
}

results=data.frame(read.table(ResFile,header=TRUE,sep="\t",stringsAsFactors=FALSE))

for (c in 4:4)
  {
	MLP=results[,c]
	top=max(MLP)
	bottom=0
	MLP=MLP[order(-MLP)]
	eMLP=qqMLPs(MLP)
	filename=sprintf("QQ.%s.png",colnames(results)[c])
	ppi=600
	png(filename,width=6*ppi, height=6*ppi, res=ppi)

	toPlot=data.frame(matrix(ncol=2,nrow=length(MLP)))
	toPlot[,1]=MLP
	toPlot[,2]=eMLP
#	colnames(toPlot)=c(colnames(results)[c],eMLPname)
	colnames(toPlot)=c("MLP","eMLP")
	myplot=ggplot(toPlot,aes_q(x=as.name("eMLP"),y=as.name("MLP")))+geom_point(size=1)+ theme_bw() + 
		geom_hline(yintercept=0,size=1.0) +
		geom_vline(xintercept=0,size=1.0) +
		theme(panel.grid.major=element_line(colour = "black",size=0.25)) +
		scale_x_continuous(breaks = seq(0,ceiling(top),by =1),minor_breaks=NULL,limits=c(0,ceiling(top))) +
		scale_y_continuous(breaks = seq(0,ceiling(top),by =1),minor_breaks=NULL,limits=c(0,ceiling(top))) 
	print(myplot)
	dev.off()
}



