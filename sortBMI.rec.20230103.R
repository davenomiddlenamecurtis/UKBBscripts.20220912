#!/share/apps/R-3.6.1/bin/Rscript

SLPFile="/cluster/project9/bipolargenomes/UKBB/UKBB.BMI.rec.20221221/UKBB.BMI.rec.20221221.summ.txt"
NumCompHetFile="/cluster/project9/bipolargenomes/UKBB/UKBB.BMI.rec.20221221/UKBB.BMI.rec.20221221.nCompHets.txt"
TestName="BMI.rec.20221221"
ResultsFolder="/cluster/project9/bipolargenomes/UKBB/UKBB.BMI.rec.20221221/results"
ResultsTemplate="UKBB.BMI.rec.20221221.%s.recsco"

quantitative=TRUE

library(ggplot2)

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


SLPs=data.frame(read.table(SLPFile,header=TRUE,sep="\t",stringsAsFactors=FALSE))
NumCompHets=data.frame(read.table(NumCompHetFile,header=TRUE,sep="\t",stringsAsFactors=FALSE))

print(sprintf("Number of genes with valid single variants = %d",nrow(SLPs)))
AllSLPs=merge(SLPs,NumCompHets,by="Gene",all=FALSE)
AllSLPs=AllSLPs[AllSLPs$numInformative>=200,]
print(sprintf("Number of genes with >=200 compound heterozygotes %d",nrow(AllSLPs)))

options(bitmapType='cairo')
SLP=AllSLPs$linrSLP
top=max(SLP)
bottom=min(SLP)
SLP=SLP[order(-SLP)]
eSLP=qqSLPs(SLP)
filename=sprintf("QQ.%s.png",TestName)
ppi=600
png(filename,width=6*ppi, height=6*ppi, res=ppi)
toPlot=data.frame(matrix(ncol=2,nrow=nrow(AllSLPs)))
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

Best=AllSLPs[abs(AllSLPs$linrSLP)>4,]
for (r in 1:nrow(Best)) {
	Gene=Best$Gene[r]
	ResFile=sprintf(ResultsTemplate,Gene)
	CommStr=sprintf("cp %s/%s %s",ResultsFolder,ResFile,ResFile)
	system(CommStr)
	Results=data.frame(read.table(ResFile,header=FALSE,sep="",stringsAsFactors=FALSE))
	colnames(Results)[1:3]=c("IID","Pheno","Score")
	Results$Carrier=(Results$Score!=0)*1
	AllData=merge(Covars,Results,by="IID",all=FALSE)
	Form=sprintf("%s Sex + Score",formulaString)
	if (quantitative) {
		m=glm(as.formula(Form), data = AllData)
	} else {
		m=glm(as.formula(Form), data = AllData, family = "binomial")
	}
	print(summary(m))
	Form=sprintf("%s Sex + Carrier",formulaString)
	if (quantitative) {
		m=glm(as.formula(Form), data = AllData)
	} else {
		m=glm(as.formula(Form), data = AllData, family = "binomial")
	}
	print(summary(m))
	Carriers=Results[Results$Score!=0,]
	print(cor.test(Carriers$Pheno,Carriers$Score,method="pearson"))
	print(cor.test(Carriers$Pheno,Carriers$Score,method="kendall"))
	png(sprintf("%s.png",ResFile),width=6*ppi, height=6*ppi, res=ppi)
	myplot=ggplot(Carriers,aes(x=Score,y=Pheno))+geom_point(size=1)+ theme_bw()
	print(myplot)
	dev.off()
}
