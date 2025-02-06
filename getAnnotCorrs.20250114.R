#!/share/apps/R-3.6.1/bin/Rscript
library(ggplot2)

# wd="C:/Users/dave/OneDrive/sharedSeq/UKBB/annot"
# setwd(wd)

geneModelsFile="geneModels.20231109.txt"
ver="forAnnot.20250114"
catsFn="categories.txt"
extraWeightsFn="extraWeights.20210413.txt"

geneModels=data.frame(read.table(geneModelsFile,header=FALSE,sep=""))
cats=data.frame(read.table(catsFn,header=FALSE,sep="",stringsAsFactors=FALSE))
vwNames=cats[cats[,1]!="Unused",1]
extraWeights=data.frame(read.table(extraWeightsFn,header=TRUE,sep="",stringsAsFactors=FALSE))
ewNames=c("GPN_MSA","AlphaMissense_Category","AlphaMissense_Score",extraWeights[,1])
ewScores=data.frame(matrix(ncol=length(ewNames),nrow=0))
colnames(ewScores)=ewNames
for (r in 1:nrow(geneModels)) {
	pheno=geneModels[r,1]
	gene=geneModels[r,2]
	firstFields=15
	nf=firstFields+2 # vepWeight + comment
	nf=nf+length(vwNames)
	nf=nf+length(ewNames)
	saoFn=sprintf("UKBB.%s.%s.%s.sao",pheno,ver,gene)
	tableFn=sprintf("UKBB.%s.%s.%s.tab",pheno,ver,gene)
	commLine=sprintf("if [ ! -e %s ] ; then tail -n +2 %s | awk ' NF==%d {print ;}' > %s; fi",tableFn,saoFn,nf,tableFn)
	system(commLine)
	tab=data.frame(read.table(tableFn,header=FALSE,sep="",stringsAsFactors=FALSE))
	tab=tab[grepl("missense",tab[,nf],fixed=TRUE),]
	tab=tab[,(firstFields+length(vwNames)+2):(firstFields+length(vwNames)+2+length(ewNames)-1)]
	colnames(tab)=ewNames
	ewScores=rbind(ewScores,tab)
}
ewScores=ewScores[rowSums(ewScores)!=0,]
write.table(ewScores,sprintf("ewScores.%s.txt",ver),row.names=FALSE,quote=FALSE,sep="\t")
ewScores=data.frame(read.table(sprintf("ewScores.%s.txt",ver),header=TRUE,sep="\t"))

library(corrplot)

options(bitmapType='cairo')

png(height=2400, width=2400, pointsize=25, file=sprintf("ewScoreCorrs.%s.png",ver))
whiteblack <- c("white", "black")
M=cor(ewScores)
corrplot(M, col = whiteblack, bg = "gray",tl.col="black",cl.pos="n")
dev.off()
corrMatrix=data.frame(matrix(ncol=length(ewNames)+1,nrow=length(ewNames)))
colnames(corrMatrix)=c("Annotation",ewNames)
corrMatrix[,1]=ewNames
corrMatrix[,2:(length(ewNames)+1)]=M
write.table(corrMatrix,sprintf("corrMatrix.%s.txt",ver),row.names=FALSE,quote=FALSE,sep="\t")


