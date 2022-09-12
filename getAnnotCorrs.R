#!/share/apps/R-3.6.1/bin/Rscript
library(ggplot2)

wd="C:/Users/dave/OneDrive/sharedSeq/UKBB/annot"
geneModelsFile="geneModels.txt"
ver="forAnnot.20210330"
catsFn="categories.txt"
extraWeightsFn="extraWeights.20210413.txt"
ver="forAnnot.20210330"

setwd(wd)

geneModels=data.frame(read.table(geneModelsFile,header=FALSE,sep=""))
cats=data.frame(read.table(catsFn,header=FALSE,sep="",stringsAsFactors=FALSE))
vwNames=cats[cats[,1]!="Unused",1]
extraWeights=data.frame(read.table(extraWeightsFn,header=TRUE,sep="",stringsAsFactors=FALSE))
ewNames=extraWeights[,1]
ewScores=data.frame(matrix(ncol=length(ewNames),nrow=0))
colnames(ewScores)=ewNames
for (r in 1:nrow(geneModels)) {
	pheno=geneModels[r,1]
	gene=geneModels[r,2]
	if (pheno=="BMI") {
		firstFields=13
	} else {
		firstFields=15
	}
	nf=firstFields+2 # vepWeight and comment
	nf=nf+length(vwNames)
	nf=nf+length(ewNames)
	saoFn=sprintf("%s.%s.%s.sao",pheno,ver,gene)
	tableFn=sprintf("%s.%s.%s.tab",pheno,ver,gene)
	commLine=sprintf("if [ ! -e %s ] ; then tail -n+3 %s | awk ' NF==%d {print ;}' > %s; fi",tableFn,saoFn,nf,tableFn)
	system(commLine)
	tab=data.frame(read.table(tableFn,header=FALSE,sep=""))
	tab=tab[,(firstFields+length(vwNames)+2):(firstFields+length(vwNames)+2+length(ewNames)-1)]
	colnames(tab)=ewNames
	ewScores=rbind(ewScores,tab)
}
ewScores=ewScores[rowSums(ewScores)!=0,]
write.table(ewScores,sprintf("ewScores.%s.txt",ver),row.names=FALSE,quote=FALSE,sep="\t")
ewScores=data.frame(read.table(sprintf("ewScores.%s.txt",ver),header=TRUE,sep="\t"))

library(corrplot)
M=cor(ewScores)
png(height=4800, width=4800, pointsize=25, res=200, file="ewScoreCorrs.png")
whiteblack <- c("white", "black")
corrplot(M, col = whiteblack, bg = "gray",tl.col="black",cl.pos="n")
dev.off()
