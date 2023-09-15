#!/share/apps/R-3.6.1/bin/Rscript

SLPFile="/cluster/project9/bipolargenomes/UKBB/UKBB.BMI.rec.20221221/UKBB.BMI.rec.20221221.summ.txt"
NumCompHetFile="/cluster/project9/bipolargenomes/UKBB/UKBB.BMI.rec.20221221/UKBB.BMI.rec.20221221.nCompHets.txt"
TestName="BMI.rec.20221221"
ResultsFolder="/cluster/project9/bipolargenomes/UKBB/UKBB.BMI.rec.20221221/results"
ResultsTemplate="UKBB.BMI.rec.20221221.%s.recsco"
TheseResultsFile="UKBB.BMI.rec.Over250.20230104.txt"
ScoreThreshold=250
quantitative=TRUE

is_try_error <- function(x) inherits(x, "try-error")

SLPs=data.frame(read.table(SLPFile,header=TRUE,sep="\t",stringsAsFactors=FALSE))
NumCompHets=data.frame(read.table(NumCompHetFile,header=TRUE,sep="\t",stringsAsFactors=FALSE))

print(sprintf("Number of genes with valid single variants = %d",nrow(SLPs)))
AllSLPs=merge(SLPs,NumCompHets,by="Gene",all=FALSE)
AllSLPs=AllSLPs[AllSLPs$numInformative>=200,]
print(sprintf("Number of genes with >=200 compound heterozygotes %d",nrow(AllSLPs)))

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

TheseResults=data.frame(matrix(ncol=4,nrow=nrow(AllSLPs),NA))
colnames(TheseResults)=c("Gene","All","Female","Male")
for (r in 1:nrow(AllSLPs)) {
	CarryOn=TRUE
	Gene=AllSLPs$Gene[r]
	print(Gene)
	ResFile=sprintf(ResultsTemplate,Gene)
	ResFile=sprintf("%s/%s",ResultsFolder,ResFile)
	Results=data.frame(read.table(ResFile,header=FALSE,sep="",stringsAsFactors=FALSE))
	colnames(Results)[1:3]=c("IID","Pheno","Score")
	Results$Carrier=(Results$Score>ScoreThreshold)*1
	AllData=merge(Covars,Results,by="IID",all=FALSE)
	Form=sprintf("%s Sex + Carrier",formulaString)
	TheseResults[r,1]=Gene
	m=NA
	if (quantitative) {
		m=try(glm(as.formula(Form), data = AllData))
	} else {
		m=glm(as.formula(Form), data = AllData, family = "binomial")
	}
	print(m)
	if (!is.na(m) & !is_try_error(m)) {
		coeffs=summary(m)$coefficients
		if (nrow(coeffs)>22) {
			TheseResults[r,2]=-log10(coeffs[23,4])
			if (coeffs[23,1]<0) {
				TheseResults[r,2]=TheseResults[r,2]*-1
			}
		}
	}
	Form=sprintf("%s Carrier",formulaString)
	Female=AllData[AllData$Sex==0,]
	m=NA
	if (quantitative) {
		m=try(glm(as.formula(Form), data = Female))
	} else {
		m=glm(as.formula(Form), data = Female, family = "binomial")
	}
	print(m)
	if (!is.na(m) & !is_try_error(m)) {
		coeffs=summary(m)$coefficients
		if (nrow(coeffs)>21) {
			TheseResults[r,3]=-log10(coeffs[22,4])
			if (coeffs[22,1]<0) {
				TheseResults[r,3]=TheseResults[r,3]*-1
			}
		}
	}
	Male=AllData[AllData$Sex==1,]
	m=NA
	if (quantitative) {
		m=try(glm(as.formula(Form), data = Male))
	} else {
		m=glm(as.formula(Form), data = Male, family = "binomial")
	}
	print(m)
	if (!is.na(m) & !is_try_error(m)) {
		coeffs=summary(m)$coefficients
		if (nrow(coeffs)>21) {
			TheseResults[r,4]=-log10(coeffs[22,4])
			if (coeffs[22,1]<0) {
				TheseResults[r,4]=TheseResults[r,4]*-1
			}
		}
	}
	write.table(TheseResults,TheseResultsFile,sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
}

write.table(TheseResults,TheseResultsFile,sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)


