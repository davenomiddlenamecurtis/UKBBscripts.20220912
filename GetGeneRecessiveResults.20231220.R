#!/share/apps/R-3.6.1/bin/Rscript

IDphenotypefile="/home/rejudcu/UKBB/lipids/UKBB.HL.20201103.txt"
TestName="HL.rec.20231220"
IDphenotypefile="/home/rejudcu/UKBB/T2D.20210104/UKBB.T2D.txt"
TestName="T2D.rec.20231220"
IDphenotypefile="/home/rejudcu/UKBB/HT/UKBB.HT.txt"
TestName="HT.rec.20231220"

SLPFile="/cluster/project9/bipolargenomes/UKBB/UKBB.HT.rec.20230728/UKBB.HT.rec.20230728.summ.txt"
ResultsFolder="/cluster/project9/bipolargenomes/UKBB/UKBB.HT.rec.20230728/results"
ResultsTemplate="UKBB.HT.rec.20230728.%s.recsco"
SaoTemplate="UKBB.HT.rec.20230728.%s.sao" # only to find if on X chromosome

is_try_error <- function(x) inherits(x, "try-error")

args = commandArgs(trailingOnly=TRUE)
if (length(args)>2) {
	TestName=args[1]
	IDphenotypefile=args[2]
	Gene=args[3]
} else {
	print("Usage: RScript GetGeneRecessiveResults.20231220.R TestName IDphenotypefile Gene")
	quit()
}

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
	AllData=merge(PhenosAndCovars,Results,by="IID",all=FALSE)
	if (XGene) {
		AllData=AllData[AllData$Sex==0,] # only females
	}
	AllData$RankScore=0
	AllData$RankScore[AllData$Score!=0]=rank(AllData$Score[AllData$Score!=0])
	Predictor=c("Score","RankScore","Carrier")
	OutputSaoFile=sprintf("Rec.%s.%s.recsao",TestName,Gene)
	sink(OutputSaoFile)
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
			writeLines(sprintf("\nLL0 = %f\nLL1 = %f\nCh2=%f\np = %f\nSLP = %f\n",LL0,LL1,ch2,pVal,SLP))
	}
	sink()
	
