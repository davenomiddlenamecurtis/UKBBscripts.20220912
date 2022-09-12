wd="C:/Users/dave/OneDrive/sharedseq/UKBB/corrected"
tgsoFile="UKBB.BMI.corrected.tgso"

library(ggplot2)

qqMLPs<-function(MLPs) {
	rankMLP <- rank(MLPs, na.last=TRUE, ties.method="first")
	nelem=length(MLPs)
	pVal=1-(rankMLP-1)/nelem
	eMLPs=-log10(pVal)
	return(eMLPs)
}

setwd(wd)
results=data.frame(read.table(tgsoFile,header=TRUE))
total=rowSums(results[,2:ncol(results)])
results=results[which(total != 0),] # remove uninformative genes
colnames(results)[2]="MLPwithPCs"
colnames(results)[3]="MLPwithoutPCs"
for (c in 2:2)
# for (c in 2:ncol(results))
  {
	MLP=na.omit(results[,c])
	top=max(MLP)
	bottom=min(MLP)
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
	pos=toPlot[which(MLP>=0),]
	pos=pos[20:nrow(pos),] # remove top 20 gene sets
	l=lm(pos[,1] ~ pos[,2])
	print(l)
}
nrow(results)
