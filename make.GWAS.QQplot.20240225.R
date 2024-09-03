#!/share/apps/R-3.6.1/bin/Rscript

# script to make a QQ plot from GWAS results
# just provide the results file name as an argument

args = commandArgs(trailingOnly=TRUE)
if (length(args)<1) {
	print("Run with GWAS results file as first argument")
	quit()
}

ResFile=args[1]
results=data.frame(read.table(ResFile,header=TRUE,stringsAsFactors=FALSE,fill=TRUE))

library(ggplot2)
options(bitmapType='cairo')

qqMLPs<-function(MLPs) {
	rankMLP <- rank(MLPs, na.last=TRUE, ties.method="first")
	nelem=length(MLPs)
	pVal=1-(rankMLP-1)/nelem
	eMLPs=-log10(pVal)
	return(eMLPs)
}

results$MLP=-log10(as.numeric(results$P))

	MLP=results$MLP
	top=max(MLP)
	bottom=0
	MLP=MLP[order(-MLP)]
	eMLP=qqMLPs(MLP)
	filename=sprintf("%s.QQ.png",ResFile)
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
		scale_x_continuous(breaks = seq(0,ceiling(max(eMLP)),by =1),minor_breaks=NULL,limits=c(0,ceiling(max(eMLP)))) +
		scale_y_continuous(breaks = seq(0,ceiling(top),by =1),minor_breaks=NULL,limits=c(0,ceiling(top))) 
	print(myplot)
	dev.off()



