#!/share/apps/R-3.6.1/bin/Rscript

IDphenotypefile="/home/rejudcu/UKBB/lipids/UKBB.HL.20201103.txt"
TestName="HL.rec.20230802"
IDphenotypefile="/home/rejudcu/UKBB/T2D.20210104/UKBB.T2D.txt"
TestName="T2D.rec.20230802"
IDphenotypefile="/home/rejudcu/UKBB/HT/UKBB.HT.txt"
TestName="HT.rec..20230802"

SLPFile="/cluster/project9/bipolargenomes/UKBB/UKBB.HT.rec.20230728/UKBB.HT.rec.20230728.summ.txt"
ResultsFolder="/cluster/project9/bipolargenomes/UKBB/UKBB.HT.rec.20230728/results"
ResultsTemplate="UKBB.HT.rec.20230728.%s.recsco"
SaoTemplate="UKBB.HT.rec.20230728.%s.sao" # only to find if on X chromosome

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

args = commandArgs(trailingOnly=TRUE)
if (length(args)>1) {
	TestName=args[1]
	IDphenotypefile=args[2]
}

SummResFile=sprintf("AllRec.%s.summ.txt",TestName)

SLPs=na.omit(data.frame(read.table(SLPFile,header=TRUE,stringsAsFactors=FALSE)))
print(sprintf("Number of genes with SLPs = %d",nrow(SLPs)))

SLPs=SLPs[SLPs$lrSLP!=0,]
print(sprintf("Number of genes with non-zero SLPs = %d",nrow(SLPs)))

PCsFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb23155.common.all.eigenvec"
SexFile="/home/rejudcu/UKBB/UKBB.sex.20201111.txt"
PCsTable=data.frame(read.table(PCsFile,header=TRUE,sep="\t"))
SexTable=data.frame(read.table(SexFile,header=TRUE,sep="\t"))

colnames(PCsTable)[1:2]=c("FID","IID")
FormulaString=c("Pheno ~")
for (p in 1:19) {
	FormulaString=sprintf("%s PC%d +",FormulaString,p)
	colnames(PCsTable)[2+p]=sprintf("PC%d",p)
}
FormulaString=sprintf("%s PC20",FormulaString)

colnames(PCsTable)[2+20]="PC20"

Covars=merge(PCsTable,SexTable,by="IID",all=FALSE)
phenoTypes=data.frame(read.table(IDphenotypefile,header=FALSE,stringsAsFactors=FALSE,sep="",fill=TRUE))
if (phenoTypes[1,1]=="IID") {
  phenoTypes=phenoTypes[2:nrow(phenoTypes),]
}
colnames(phenoTypes)=c("IID","Pheno")
phenoTypes=phenoTypes[complete.cases(phenoTypes),]
phenoTypes$Pheno=as.numeric(phenoTypes$Pheno)
phenoTypes=phenoTypes[which(phenoTypes$Pheno==0 | phenoTypes$Pheno==1),]

PhenosAndCovars=merge(phenoTypes,Covars,by="IID",all=FALSE)

if (file.exists(SummResFile)) {
	SummRes=data.frame(read.table(SummResFile,header=TRUE,sep="\t",stringsAsFactors=FALSE))
} else {
	SummRes=data.frame(matrix(NA, ncol = 5, nrow = nrow(SLPs)),stringsAsFactors=FALSE)
	colnames(SummRes)=c("Gene","NumCarriers","ScorePhenoSLP","RankScorePhenoSLP","CarrierPhenoSLP")
	SummRes$Gene=SLPs$Gene
}

Gene=SummRes$Gene[1] # some subjects do not have a score, not sure why - maybe no BMI
ResFile=sprintf(ResultsTemplate,Gene)
Results=data.frame(read.table(sprintf("%s/%s",ResultsFolder,ResFile),header=FALSE,sep="",stringsAsFactors=FALSE))
colnames(Results)[1:3]=c("IID","OriginalPheno","Score")
AllData=merge(PhenosAndCovars,Results,by="IID",all=FALSE)
Form=sprintf("%s + Sex",FormulaString)
m=try(glm(as.formula(Form), data = AllData, family = "binomial"))
AllLL0=as.numeric(logLik(m))
FemaleData=AllData[AllData$Sex==0,]
Form=FormulaString
m=try(glm(as.formula(FormulaString), data = FemaleData, family = "binomial"))
FemaleLL0=as.numeric(logLik(m))
for (r in 1:nrow(SummRes)) {
	if (!is.na(SummRes$RankScorePhenoSLP[r])) {
		next # done this gene already
	}
	Gene=SummRes$Gene[r]
	print(Gene)
	SaoFile=sprintf(SaoTemplate,Gene)
	SaoLines=readLines(sprintf("%s/%s",ResultsFolder,SaoFile),n=3)
	Coord=strsplit(SaoLines[3],":")[[1]][1]
	if (Coord=="23") {
		XGene=TRUE
	} else {
		XGene=FALSE
	}
	ResFile=sprintf(ResultsTemplate,Gene)
	Results=data.frame(read.table(sprintf("%s/%s",ResultsFolder,ResFile),header=FALSE,sep="",stringsAsFactors=FALSE))
	colnames(Results)[1:3]=c("IID","OriginalPheno","Score")
	Results$Carrier=(Results$Score>0)*1
	SummRes$NumCarriers[r]=sum(Results$Carrier)
	AllData=merge(PhenosAndCovars,Results,by="IID",all=FALSE)
	if (XGene) {
		AllData=AllData[AllData$Sex==0,] # only females
	}
	AllData$RankScore=0
	AllData$RankScore[AllData$Score!=0]=rank(AllData$Score[AllData$Score!=0])
	Predictor=c("Score","RankScore","Carrier")
	for (p in 1:3) {
			if (XGene) {
				LL0=FemaleLL0
				Form=sprintf("%s + %s",FormulaString,Predictor[p])
			} else {
				LL0=AllLL0
				Form=sprintf("%s + %s + Sex",FormulaString,Predictor[p])
			}
			m=try(glm(as.formula(Form), data = AllData, family = "binomial"))
			if (is_try_error(m)) {
				next
			}
			LL1=as.numeric(logLik(m))
			print(summary(m)$coefficients)
			ch2=2*(LL1-LL0)
			pVal=pchisq(ch2,1,lower.tail=FALSE)
			SLP=-log10(pVal)*sign(summary(m)$coefficients[22,1])
			SummRes[r,2+p]=SLP
	}

write.table(SummRes,SummResFile,sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
}
write.table(SummRes,SummResFile,sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

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
}
}