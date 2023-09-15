#!/share/apps/R-3.6.1/bin/Rscript

SLPFile="/cluster/project9/bipolargenomes/UKBB/UKBB.BMI.rec.20230108/UKBB.BMI.rec.20230108.summ.txt"
TestName="BMI.AllRec.20230122"
ResultsFolder="/cluster/project9/bipolargenomes/UKBB/UKBB.BMI.rec.20230108/results"
ResultsTemplate="UKBB.BMI.rec.20230108.%s.recsco"
SaoTemplate="UKBB.BMI.rec.20230108.%s.sao"
SummResFile=sprintf("AllRec.%s.summ.txt",TestName)

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

PCsFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb23155.common.all.eigenvec"
SexFile="/home/rejudcu/UKBB/UKBB.sex.20201111.txt"
PCsTable=data.frame(read.table(PCsFile,header=TRUE,sep="\t"))
SexTable=data.frame(read.table(SexFile,header=TRUE,sep="\t"))

colnames(PCsTable)[1:2]=c("FID","IID")
FormulaString=c("Pheno ~","Obese ~")
for (p in 1:19) {
	for (f in 1:2) {
		FormulaString[f]=sprintf("%s PC%d +",FormulaString[f],p)
	}
	colnames(PCsTable)[2+p]=sprintf("PC%d",p)
}
for (f in 1:2) {
	FormulaString[f]=sprintf("%s PC20",FormulaString[f])
}
colnames(PCsTable)[2+20]="PC20"

Covars=merge(PCsTable,SexTable,by="IID",all=FALSE)

if (file.exists(SummResFile)) {
	SummRes=data.frame(read.table(SummResFile,header=TRUE,sep="\t",stringsAsFactors=FALSE))
} else {
	SummRes=data.frame(matrix(NA, ncol = 7, nrow = nrow(SLPs)),stringsAsFactors=FALSE)
	colnames(SummRes)=c("Gene","NumCarriers","ScoreBMISLP","RankScoreBMISLP","CarrierBMISLP","RankScoreObeseSLP","CarrierObeseSLP")
	SummRes$Gene=SLPs$Gene
	SummRes$ScoreBMISLP=SLPs$linrSLP
}

# get BMIs for LL0 from first score file
Gene=SummRes$Gene[1]
ResFile=sprintf(ResultsTemplate,Gene)
Results=data.frame(read.table(sprintf("%s/%s",ResultsFolder,ResFile),header=FALSE,sep="",stringsAsFactors=FALSE))
colnames(Results)[1:3]=c("IID","Pheno","Score")
Results$Obese=(Results$Pheno>=30)*1
AllData=merge(Covars,Results,by="IID",all=FALSE)
AllLL0=numeric(length=2)
FemaleLL0=numeric(length=2)
for (f in 1:2) {
	if (f==2) {
		quantitative=FALSE
	} else{
		quantitative=TRUE
	}
	Form=sprintf("%s + Sex",FormulaString[f])
	if (quantitative) {
		m=try(glm(as.formula(Form), data = AllData))
	} else {
		m=try(glm(as.formula(Form), data = AllData, family = "binomial"))
	}
	AllLL0[f]=as.numeric(logLik(m))
	FemaleData=AllData[AllData$Sex==0,]
	Form=FormulaString[f]
	if (quantitative) {
		m=try(glm(as.formula(Form), data = FemaleData))
	} else {
		m=try(glm(as.formula(Form), data = FemaleData, family = "binomial"))
	}
	FemaleLL0[f]=as.numeric(logLik(m))
}

for (r in 1:nrow(SummRes)) {
	if (!is.na(SummRes$RankScoreBMISLP[r])) {
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
	colnames(Results)[1:3]=c("IID","Pheno","Score")
	Results$Obese=(Results$Pheno>=30)*1
	Results$Carrier=(Results$Score>0)*1
	AllData=merge(Covars,Results,by="IID",all=FALSE)
	if (XGene) {
		AllData=AllData[AllData$Sex==0,]
	}
	AllData$RankScore=0
	AllData$RankScore[AllData$Score!=0]=rank(AllData$Score[AllData$Score!=0])
	for (f in 1:2) {
		if (f==2) {
			quantitative=FALSE
		} else {
			quantitative=TRUE # obese
		}
		Predictor=c("RankScore","Carrier")
		for (p in 1:2) {
			if (XGene) {
				LL0=FemaleLL0[f]
				Form=sprintf("%s + %s",FormulaString[f],Predictor[p])
			} else {
				LL0=AllLL0[f]
				Form=sprintf("%s + %s + Sex",FormulaString[f],Predictor[p])
			}
			if (quantitative) {
				m=try(glm(as.formula(Form), data = AllData))
			} else {
				m=try(glm(as.formula(Form), data = AllData, family = "binomial"))
			}
			LL1=as.numeric(logLik(m))
			print(summary(m)$coefficients)
			ch2=2*(LL1-LL0)
			pVal=pchisq(ch2,1,lower.tail=FALSE)
			SLP=-log10(pVal)*sign(summary(m)$coefficients[22,1])
			SummRes[r,1+2*f+p]=SLP
		}
	}

write.table(SummRes,SummResFile,sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
}
write.table(SummRes,SummResFile,sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

options(bitmapType='cairo')
for (MinCarrier in c(1,20,50)) {
for (s in 1:5) {
Narrowed=SummRes[SummRes$NumCarriers>=MinCarrier,]
SLP=na.omit(Narrowed[,2+s])
top=max(SLP)
bottom=min(SLP)
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