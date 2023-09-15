#!/share/apps/R-3.6.1/bin/Rscript

SLPFile="/cluster/project9/bipolargenomes/UKBB/UKBB.BMI.rec.20230108/UKBB.BMI.rec.20230108.summ.txt"
TestName="BMI.rec.best.20230108"
ResultsFolder="/cluster/project9/bipolargenomes/UKBB/UKBB.BMI.rec.20230108/results"
ResultsTemplate="UKBB.BMI.rec.20230108.%s.recsco"

InterestingGenes=c("MC4R","PCSK1")
quantitative=TRUE

library(ggplot2)

is_try_error <- function(x) inherits(x, "try-error")

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

SLPs=na.omit(data.frame(read.table(SLPFile,header=TRUE,stringsAsFactors=FALSE)))
print(sprintf("Number of genes with SLPs = %d",nrow(SLPs)))

SLPs=SLPs[SLPs$linrSLP!=0,]
print(sprintf("Number of genes with non-zero SLPs = %d",nrow(SLPs)))

options(bitmapType='cairo')
SLP=SLPs$linrSLP
top=max(SLP)
bottom=min(SLP)
SLP=SLP[order(-SLP)]
eSLP=qqSLPs(SLP)
filename=sprintf("QQ.%s.png",TestName)
ppi=600
png(filename,width=6*ppi, height=6*ppi, res=ppi)
toPlot=data.frame(matrix(ncol=2,nrow=nrow(SLPs)))
toPlot[,1]=SLP
toPlot[,2]=eSLP
colnames(toPlot)=c("SLP","eSLP")
myplot=ggplot(toPlot,aes_q(x=eSLP,y=SLP))+geom_point(size=1)+ theme_bw() +
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
print(l)

PCsFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb23155.common.all.eigenvec"
SexFile="/home/rejudcu/UKBB/UKBB.sex.20201111.txt"
PCsTable=data.frame(read.table(PCsFile,header=TRUE,sep="\t"))
SexTable=data.frame(read.table(SexFile,header=TRUE,sep="\t"))

colnames(PCsTable)[1:2]=c("FID","IID")
formulaString=("Pheno ~")
for (p in 1:20) {
  formulaString=sprintf("%s PC%d +",formulaString,p)
  colnames(PCsTable)[2+p]=sprintf("PC%d",p)
}

Covars=merge(PCsTable,SexTable,by="IID",all=FALSE)

Best=SLPs[abs(SLPs$linrSLP)>3,]$Gene
InterestingGenes=c(Best,InterestingGenes)
SummRes=data.frame(matrix(NA, ncol = 6, nrow = length(InterestingGenes)),stringsAsFactors=FALSE)
colnames(SummRes)=c("Gene","NumCarriers","ScoreSLP","CarrierSLP","PccSLP","KccSLP")
SummRes[,1]=InterestingGenes
r=0
for (Gene in InterestingGenes) {
	r=r+1
	print(Gene)
	ResFile=sprintf(ResultsTemplate,Gene)
	CommStr=sprintf("cp %s/%s %s",ResultsFolder,ResFile,ResFile)
	system(CommStr)
	Results=data.frame(read.table(ResFile,header=FALSE,sep="",stringsAsFactors=FALSE))
	colnames(Results)[1:3]=c("IID","Pheno","Score")
	Results$Carrier=(Results$Score!=0)*1
	AllData=merge(Covars,Results,by="IID",all=FALSE)
	Form=sprintf("%s Sex + Score",formulaString)
	if (quantitative) {
		m=try(glm(as.formula(Form), data = AllData))
	} else {
		m=try(glm(as.formula(Form), data = AllData, family = "binomial"))
	}
	print(summary(m))
	SummRes$ScoreSLP[r]=log10(coef(summary(m))[23,4])
	if (coef(summary(m))[23,1]>0) {
		SummRes$ScoreSLP[r]=SummRes$ScoreSLP[r]*-1
	}
	Form=sprintf("%s Sex + Carrier",formulaString)
	if (quantitative) {
		m=try(glm(as.formula(Form), data = AllData))
	} else {
		m=try(glm(as.formula(Form), data = AllData, family = "binomial"))
	}
	print(summary(m))
	SummRes$CarrierSLP[r]=log10(coef(summary(m))[23,4])
	if (coef(summary(m))[23,1]>0) {
		SummRes$CarrierSLP[r]=SummRes$CarrierSLP[r]*-1
	}
	Carriers=Results[Results$Score!=0,]
	SummRes$NumCarriers[r]=nrow(Carriers)
	if (nrow(Carriers)<=1) {
		print(sprintf("Only one valid carrier for gene %s", Gene))
		next
	}
	png(sprintf("%s.png",ResFile),width=6*ppi, height=6*ppi, res=ppi)
	myplot=ggplot(Carriers,aes(x=Score,y=Pheno))+geom_point(size=1)+ theme_bw() +
                scale_x_continuous(breaks = seq(0,1500,by =250),minor_breaks=NULL,limits=c(0,1500)) +
                scale_y_continuous(breaks = seq(10,70,by =10),minor_breaks=NULL,limits=c(10,70))

	print(myplot)
	dev.off()
	if (nrow(Carriers)==2) {
		print(sprintf("Only two valid carriers for gene %s", Gene))
		next
	}
	if (sd(Carriers$Score)==0) {
		print(sprintf("All carriers have the same score for gene %s", Gene))
		next
	}
	Pcc=cor.test(Carriers$Pheno,Carriers$Score,method="pearson")
	SummRes$PccSLP[r]=log10(Pcc$p.value)
	if (Pcc$estimate>0) {
		SummRes$PccSLP[r]=SummRes$PccSLP[r]*-1
	}
	Kcc=cor.test(Carriers$Pheno,Carriers$Score,method="kendall")
	SummRes$KccSLP[r]=log10(Kcc$p.value)
	if (Kcc$estimate>0) {
		SummRes$KccSLP[r]=SummRes$KccSLP[r]*-1
	}
}
write.table(SummRes,sprintf("best.%s.summ.txt",TestName),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
