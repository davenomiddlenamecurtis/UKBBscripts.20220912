#!/share/apps/R-3.6.1/bin/Rscript
library(ggplot2)
library(corrplot)

wd="/home/rejudcu/UKBB/annot.2025"
setwd(wd)
geneModelsFile="geneModels.20231109.txt"
ver="forAnnot.20250106"
catsFn="categories.txt"
extraWeightsFn="extraWeights.20210413.txt"
protectiveGenes=c("NPC1L1","PCSK9","ANGPTL3","APOC3","INPPL1","DBH")

options(bitmapType='cairo')

geneModels=data.frame(read.table(geneModelsFile,header=FALSE,sep="",stringsAsFactors=FALSE))
cats=data.frame(read.table(catsFn,header=FALSE,sep="",stringsAsFactors=FALSE))
vwNames=cats[cats[,1]!="Unused",1]
extraWeights=data.frame(read.table(extraWeightsFn,header=TRUE,sep="",stringsAsFactors=FALSE))
ewNames=c("AlphaMissenseCategory","AlphaMissenseScore","MINUSGPN",extraWeights[,1])
pheno=geneModels[1,1]
gene=geneModels[1,2]
resultsFile=sprintf("results.%s.%s.%s.txt",ver,pheno,gene)
firstResult=data.frame(read.table(resultsFile,header=TRUE,sep="\t",stringsAsFactors=FALSE))
firstResult[9:10,1]=c("AlphaMissense_Category","AlphaMissense_Score")

allResults=data.frame(matrix(ncol=nrow(geneModels)+1,nrow=nrow(firstResult)))
colnames(allResults)[1]="Annotation"
allResults[,1]=firstResult[,1]
colnames(allResults)[1]=colnames(firstResult)[1]
allMultResults=data.frame(matrix(ncol=nrow(geneModels)+1,nrow=nrow(firstResult)))
colnames(allMultResults)[1]="Annotation"
allMultResults[,1]=firstResult[,1]
colnames(allMultResults)[1]=colnames(firstResult)[1]
for (r in 1:nrow(geneModels)) {
	pheno=geneModels[r,1]
	gene=geneModels[r,2]
	resultsFile=sprintf("results.%s.%s.%s.txt",ver,pheno,gene)
	results=data.frame(read.table(resultsFile,header=TRUE,sep="\t"))
	allResults[,r+1]=results$SLP
	colnames(allResults)[r+1]=gene
	allMultResults[,r+1]=results$multSLP
	colnames(allMultResults)[r+1]=gene
	if (gene %in% protectiveGenes) {
		allResults[,r+1]=allResults[,r+1]*-1
		allMultResults[,r+1]=allMultResults[,r+1]*-1
	}
}
allResults
write.table(allResults,sprintf("allAnnotSingleSLPs.%s.txt",ver),row.names=FALSE,quote=FALSE,sep="\t")
write.table(allMultResults,sprintf("allAnnotMultSLPs.%s.txt",ver),row.names=FALSE,quote=FALSE,sep="\t")

vwToPlot=data.matrix(allResults[1:length(vwNames),-1])
rownames(vwToPlot)=allResults[1:length(vwNames),1]
ewToPlot=data.matrix(allResults[(1+length(vwNames)):(length(vwNames)+length(ewNames)),-1])
rownames(ewToPlot)=allResults[(1+length(vwNames)):(length(vwNames)+length(ewNames)),1]
for (c in 1:ncol(vwToPlot)) {
	vwToPlot[,c]=vwToPlot[,c]/max(abs(vwToPlot[,c]))
	ewToPlot[,c]=ewToPlot[,c]/max(abs(ewToPlot[,c]))
}
vwToPlot
whiteblack <- c("white", "black")
png(height=600, width=600, pointsize=25, file=sprintf("vwSLPs.%s.png",ver))
corrplot(vwToPlot,is.corr=FALSE, col = whiteblack, bg = "gray",tl.col="black",cl.pos="n")
dev.off()
png(height=2400, width=1200, pointsize=25, file=sprintf("ewSLPs.%s.png",ver))
corrplot(ewToPlot,is.corr=FALSE, col = whiteblack, bg = "gray",tl.col="black",cl.pos="n")
dev.off()

means=rowMeans(allResults[,2:ncol(allResults)])
allResults$Average=means
write.table(allResults,sprintf("allAnnotSingleSLPs.withMean.%s.txt",ver),row.names=FALSE,quote=FALSE,sep="\t")

# toExclude=c("PCSK1", "ANGPLT3","GIGYF1" )
# ewToPlot=ewToPlot[,! colnames(ewToPlot) %in% toExclude]
# png(height=600, width=600, pointsize=25, file="ewCorrs.png")
# corrplot(ewToPlot,is.corr=FALSE, col = whiteblack, bg = "gray",tl.col="black",cl.pos="n")
# dev.off()
topScorers=vector()
for (gg in 1:ncol(ewToPlot)) {
	topScorers[gg]=rownames(ewToPlot)[abs(ewToPlot[,gg])==1]
}
topScorers
