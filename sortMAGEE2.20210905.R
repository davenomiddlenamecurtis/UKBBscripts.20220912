#!/share/apps/R-3.6.1/bin/Rscript

library(data.table)

genes="MAGEE2"
date="20210905"

model=sprintf("psychosocial.%s",date)

wd="C:/Users/dave/OneDrive/sharedseq/UKBB/MAGEE2"
wd="C:/Users/dave_000/OneDrive/sharedseq/UKBB/MAGEE2"
# setwd(wd)

PhenotypesFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb41465.exomes.20201103.txt"

SummarySchemeFile=sprintf("psychosocialFields.%s.txt",date)
OutputFile=sprintf("summary.%s.txt",date)
GenotypesFile="rs1343879.genotypes.raw"
NeededPhenotypesFile=sprintf("phenotypes.%s.txt",date)
PCsFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb23155.common.all.eigenvec"

SummaryScheme=data.frame(read.table(SummarySchemeFile,header=TRUE,sep="\t",stringsAsFactors=FALSE))

if (!file.exists(NeededPhenotypesFile)) {
	cmd=sprintf("head -n 1 %s > fields.txt",PhenotypesFile)
	system(cmd)
	fieldNames=data.frame(read.table("fields.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE))
	cmd="cut -f 1"
	for (r in 1:nrow(SummaryScheme)) {
		toMatch=sprintf("X%s",chartr("-",".",SummaryScheme$Field[r]))
		c=which(colnames(fieldNames)==toMatch)
		cmd=sprintf("%s,%d",cmd,c)
	}
	cmd=sprintf("%s %s > %s",cmd,PhenotypesFile,NeededPhenotypesFile)
	system(cmd)
}
phenotypes=data.frame(read.table(NeededPhenotypesFile,header=TRUE,sep="\t",stringsAsFactors=FALSE))

if (! file.exists(GenotypesFile)) {
	system("echo 23 75784694 75784694 rs1343879 > rs1343879.keep.txt")
	system("plink --bed /SAN/ugi/UGIbiobank/data/downloaded/ukb23155_cX_b0_v1.bed  --fam /SAN/ugi/UGIbiobank/data/downloaded/ukb23155_cX_b0_v1_s200632.20210202.fam --bim /SAN/ugi/UGIbiobank/data/downloaded/UKBexomeOQFE_chrX.bim --extract 'range' rs1343879.keep.txt --recode A --out rs1343879.genotypes")
}

genotypes=data.frame(na.omit(read.table(GenotypesFile,header=TRUE,sep=" ",stringsAsFactors=FALSE)))
colnames(genotypes)[7]="geno"
genotypes$LOF=0
genotypes$LOF[(genotypes$SEX==1) & (genotypes$geno!=0)]=1
genotypes$LOF[(genotypes$SEX==2) & (genotypes$geno==2)]=1

colnames(phenotypes)[1]="IID"
PCsTable=data.frame(read.table(PCsFile,header=FALSE,sep="\t"))
colnames(PCsTable)[1:2]=c("FID","IID")
PCs=("")
for (p in 1:20) {
  PCs=sprintf("%s PC%d +",PCs,p)
  colnames(PCsTable)[2+p]=sprintf("PC%d",p)
}

allData=merge(phenotypes,genotypes,by="IID",all=FALSE)
allData=merge(allData,PCsTable,by="IID",all=FALSE)

largestCoding=200000
codings= vector("list", largestCoding)

rowStart=vector()
rNames=vector()
rNames[1]="N"
totRow=1
for (r in 1:nrow(SummaryScheme)) {
	totRow=totRow+1
	rowStart[r]=totRow
	rNames[totRow]=SummaryScheme$Name[r]
	if (SummaryScheme$Coding[r]!=0) {
		fn=sprintf("coding%d.tsv",SummaryScheme$Coding[r])
		coding=data.frame(read.table(fn,header=TRUE,sep="\t",quote="",stringsAsFactors=FALSE))
		coding=coding[coding$coding>=SummaryScheme$Minimum[r] & coding$coding<=SummaryScheme$Maximum[r],]
	    codings[[SummaryScheme$Coding[r]]]=coding
		for (rr in 1:nrow(coding)) {
			totRow=totRow+1
			rNames[totRow]=coding$meaning[rr]
		}
	}
}

sexString=c("male","female")
ancString=c("All","SouthAsian","British","Chinese")
for (anc in 1:4) {
for (sx in 1:2) {
SummaryTab=data.frame(matrix(ncol=5,nrow=length(rNames),0))
colnames(SummaryTab)=c("Variable","Normal","LOF","FreqLOF","beta")
SummaryTab[,1]=rNames
SummaryTab$FreqLOF=""
SummaryTab$beta=""

sexToUse=allData[allData$SEX==sx,]
sexToUse=sexToUse[sexToUse[,2]==2-sx,] # sex should be 1 for male, 2 for female
if (anc!=1) {
	sexToUse=sexToUse[!is.na(sexToUse$X21000.0.0),]
}
if (anc==2) {
	sexToUse=sexToUse[sexToUse$X21000.0.0=="3001" | sexToUse$X21000.0.0=="3002",]
}
if (anc==3) {
	sexToUse=sexToUse[sexToUse$X21000.0.0=="1001",]
}
if (anc==4) {
	sexToUse=sexToUse[sexToUse$X21000.0.0=="5",]
}
for (c in 1:2) {
	SummaryTab[1,1+c]=sum(sexToUse$LOF==c-1)
}
SummaryTab[1,4]=sprintf("%.4f",SummaryTab[1,3]/(SummaryTab[1,2]+SummaryTab[1,3]))
for (r in 1:nrow(SummaryScheme)) {
	toMatch=sprintf("X%s",chartr("-",".",SummaryScheme$Field[r])) # because R converts column names
	if (SummaryScheme$Coding[r]!=0) {
		codes=codings[[SummaryScheme$Coding[r]]]$coding
	}
	st=which(colnames(sexToUse)==toMatch)
	toUse=sexToUse[!is.na(sexToUse[,st]),]
	toUse=toUse[toUse[,st]>=SummaryScheme$Minimum[r]&toUse[,st]<=SummaryScheme$Maximum[r],]
	if (SummaryScheme$Coding[r]==0) {
		for (c in 1:2) {
			m=mean(na.omit(toUse[toUse$LOF==c-1,st]))
			std=sd(na.omit(toUse[toUse$LOF==c-1,st]))
			SummaryTab[rowStart[r],1+c]=sprintf("%.2f (%.2f)",m,std)
		}
		if (r!=1) {# do not try to fit sex
			formulaString=sprintf("%s ~ %s LOF",toMatch,PCs)
			fittedModel=glm(as.formula(formulaString), data = toUse)
			coeff=summary(fittedModel)$coefficients
			SummaryTab$beta[rowStart[r]]=sprintf("%.3f (SE=%.3f), p=%f",coeff[22,1],coeff[22,2],coeff[22,4])
		}
	} else {
		for (rr in 1:length(codes)) {
			SummaryTab[rowStart[r],2:3]=""
			for (c in 1:2) {
				SummaryTab[rowStart[r]+rr,1+c]=sprintf("%.4f",sum(toUse[,st]==codes[rr] & toUse$LOF==c-1)/sum(toUse$LOF==c-1))
			}
			SummaryTab[rowStart[r]+rr,4]=sprintf("%.4f",sum(toUse[,st]==codes[rr] & toUse$LOF==1)/(sum(toUse[,st]==codes[rr] & toUse$LOF==0)+sum(toUse[,st]==codes[rr] & toUse$LOF==1)))
			toUse$outcome=0
			toUse$outcome[toUse[,st]==codes[rr]]=1
			formulaString=sprintf("outcome ~ %s LOF",PCs)
			fittedModel=glm(as.formula(formulaString), data = toUse, family="binomial")
			coeff=summary(fittedModel)$coefficients
			SummaryTab$beta[rowStart[r]+rr]=sprintf("%.3f (SE=%.3f), p=%f",coeff[22,1],coeff[22,2],coeff[22,4])
		}
	}
}

write.table(SummaryTab,sprintf("summary.%s.%s.%s.txt",sexString[sx],ancString[anc],model),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
}

}

eth="X21000.0.0" # row of Ethnicity
ethRow=2
qual="X6138.0.0" # field for Qualifications
degVal=1
codes=codings[[SummaryScheme$Coding[ethRow]]]$coding
for (sx in 1:2) {
degreeTab=data.frame(matrix(ncol=5,nrow=length(codes),0))
colnames(degreeTab)=c("Ethnicity","Normal","LOF","FreqLOF","beta")
degreeTab$Ethnicity=codings[[SummaryScheme$Coding[ethRow]]]$meaning
sexToUse=allData[allData$SEX==sx,]
sexToUse=sexToUse[sexToUse[,2]==2-sx,] # sex should be 1 for male, 2 for female
sexToUse=sexToUse[!is.na(sexToUse[[eth]]),]
sexToUse=sexToUse[!is.na(sexToUse[[qual]]),]
sexToUse=sexToUse[sexToUse[[eth]]>=SummaryScheme$Minimum[ethRow]&sexToUse[[eth]]<=SummaryScheme$Maximum[ethRow],]
for (rr in 1:length(codes)) {
	toUse=sexToUse[sexToUse[[eth]]==codes[rr],]
	degTab=toUse[toUse[[qual]]==degVal,]
	for (c in 1:2) {
		degreeTab[rr,1+c]=sprintf("%.4f",sum(degTab$LOF==c-1)/max(sum(toUse$LOF==c-1),1)) # avoid divide by zero
	}
	degreeTab[rr,4]=sprintf("%.4f",sum(degTab$LOF==1)/max((sum(degTab$LOF==0)+sum(degTab$LOF==1)),1))
	if (sum(degTab$LOF==0) | sum(degTab$LOF==1)) {
		degreeTab$beta[rr]=""
		next
	}
	toUse$outcome=0
	toUse$outcome[toUse[[qual]]==degVal]=1
	formulaString=sprintf("outcome ~ %s LOF",PCs)
	fittedModel=glm(as.formula(formulaString), data = toUse, family="binomial")
	coeff=summary(fittedModel)$coefficients
	degreeTab$beta[rr]=sprintf("%.3f (SE=%.3f), p=%f",coeff[22,1],coeff[22,2],coeff[22,4])
}
}
