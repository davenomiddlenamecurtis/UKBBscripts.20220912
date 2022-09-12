#!/share/apps/R-3.6.1/bin/Rscript
library(ggplot2)

wd="C:/Users/dave/OneDrive/sharedSeq/UKBB/annot"
geneModelsFile="geneModels.txt"
ver="forAnnot.20210330"
catsFn="categories.txt"
extraWeightsFn="extraWeights.20210413.txt"

weightsToPlot=c("LOF","ProteinAltering")

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
#    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
#                    ncol = cols, nrow = ceiling(numPlots/cols))
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols),byrow=TRUE)
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

setwd(wd)

geneModels=data.frame(read.table(geneModelsFile,header=FALSE,sep=""))
cats=data.frame(read.table(catsFn,header=FALSE,sep="",stringsAsFactors=FALSE))
vwNames=cats[cats[,1]!="Unused",1]
extraWeights=data.frame(read.table(extraWeightsFn,header=TRUE,sep="",stringsAsFactors=FALSE))
ewNames=extraWeights[,1]

pheno=geneModels[1,1]
gene=geneModels[1,2]
resultsFile=sprintf("results.%s.%s.%s.1.00.txt",ver,pheno,gene)
firstResult=data.frame(read.table(resultsFile,header=TRUE,sep="\t"))
weightsToPlot=vwNames
allResults=data.frame(matrix(0,ncol=2+length(weightsToPlot),nrow=nrow(geneModels)*length(seq(0,2,0.2))))
colnames(allResults)=c("Gene","WF",weightsToPlot)
catTables=list() 
for (t in 1:length(vwNames)) {
	catTables[[t]]=data.frame(matrix(ncol=1+nrow(geneModels),nrow=length(seq(0,2,0.2))))
	catTables[[t]][,1]=sprintf("%.1f",seq(0,2,0.2))
	colnames(catTables[[t]])[1]="log10(WF)"
	}
rr=1
for (r in 1:nrow(geneModels)) {
	pheno=geneModels[r,1]
	gene=geneModels[r,2]
for (lw in seq(0,2,0.2)) {
wf=sprintf("%.2f",10^lw)
	resultsFile=sprintf("results.%s.%s.%s.%s.txt",ver,pheno,gene,wf)
	results=data.frame(read.table(resultsFile,header=TRUE,sep="\t"))
	allResults$Gene[rr]=gene
	allResults$WF[rr]=wf
	for (w in 1:length(vwNames)) {
		allResults[rr,2+w]=results$SLP[results$Weight==weightsToPlot[w]]
		if (gene=="PCSK9" || gene=="ANGPTL3") {
			allResults[rr,2+w]=allResults[rr,2+w]*-1	
		}
	}
	rr=rr+1
}
	for (t in 1:length(vwNames)) {
		colnames(catTables[[t]])[r+1]=gene
		catTables[[t]][,r+1]=allResults[allResults$Gene==gene,t+2]
	}
}
allResults


write.table(allResults,sprintf("allWeightEffects.%s.txt",ver),row.names=FALSE,quote=FALSE,sep="\t")
for (t in 1:length(vwNames)) {
		catTables[[t]][,r+1]=allResults[allResults$Gene==gene,t+2]
		write.table(catTables[[t]],sprintf("effects.%s.%s.txt",vwNames[t],ver),row.names=FALSE,quote=FALSE,sep="\t")
	}

allResults$WF=as.numeric(allResults$WF)

plotList=list()
for (w in 1:length(weightsToPlot)) {
	plotList[[w]]=ggplot(allResults,aes_string(x="WF",y=colnames(allResults)[2+w],group="Gene",color="Gene")) +geom_line() +geom_point(aes(shape=Gene))+ theme_bw()
}
multiplot(plotlist=plotList,cols=2)

genesToPlot=allResults[allResults$Gene %in% c("LDLR","MC4R","PCSK1") ,c(1,2,7)]
colnames(genesToPlot)[3]="SLP"
ggplot(genesToPlot,aes_string(x="WF",y="SLP",group="Gene")) +geom_line() +geom_point(aes(shape=Gene))+ theme_bw()

png(height=2400, width=2400, pointsize=25, res=600, file="threeGenesSLPsByWF.png")
ggplot(genesToPlot,aes_string(x="WF",y="SLP",group="Gene")) +geom_line() +geom_point(size=2,aes(shape=Gene))+ theme_bw()
dev.off()

