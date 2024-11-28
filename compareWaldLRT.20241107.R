#!/share/apps/R-3.6.1/bin/Rscript
# script to get data files to analyse contribution of different types of variant to phenotype

wd=""/home/rejudcu/UKBB/WaldTest"
setwd(wd)

PhenoFile="/home/rejudcu/UKBB/psoriasis.22K.20241104/UKBB.psoriasis.22K.txt"

PCsFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb23158.common.all.20230806.eigenvec.txt"
SexFile="/home/rejudcu/UKBB/UKBB.sex.20230807.txt"
GWASFileTemplate="/home/rejudcu/UKBB/psoriasis.22K.20241104/UKBB.psoriasis.22K.txt.results.%d.PHENO2.glm.logistic"
SNPTemplate="/SAN/ugi/UGIbiobank/data/downloaded/ukb_cal_chr%d_v2"
SNPNamesTemplate="SNPsToExtract.%d.txt"
GenosFileTemplate="ExtractedGenos.%d"
ResultsFile="UKBB.psoriasis.22K.Wald_v_LR.txt"
MLPThreshold=4

for (chr in 1:22) {
	GWASFile=sprintf(GWASFileTemplate,chr)
	SNPNamesFile=sprintf(SNPNamesTemplate,chr)
	GenosFile=sprintf(GenosFileTemplate,chr)
	SNPs=sprintf(SNPTemplate,chr)
	GWAS=na.omit(data.frame(read.table(GWASFile,header=TRUE,stringsAsFactors=FALSE,comment.char="")))
	colnames(GWAS)[1]="CHROM"
	GWAS$ID=gsub("-",".",GWAS$ID)
	GWAS$MLP=-log10(GWAS$P)
	GWAS=GWAS[GWAS$MLP>MLPThreshold,]
	GWAS$SLP=GWAS$MLP*sign(GWAS$Z_STAT)
	SNPsNames=GWAS[,3,drop=FALSE]
	if (nrow(SNPsNames)==0) {
		next
	}
	write.table(SNPsNames,SNPNamesFile,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
	CommStr=sprintf("/share/apps/genomics/plink-2.0/bin/plink2 --bfile %s --extract %s --recode A --out %s",SNPs,SNPNamesFile,GenosFile)
	print(CommStr)
	system(CommStr)
	Genos=data.frame(read.table(sprintf("%s.raw",GenosFile),header=TRUE,stringsAsFactors=FALSE))
	Genos=Genos[,c(-1,-3,-4,-5,-6)]
	if (chr==1) {
		AllGenos=Genos
		AllGWAS=GWAS
	} else {
		AllGenos=merge(AllGenos,Genos,all=TRUE,by="IID")
		AllGWAS=rbind(AllGWAS,GWAS)
	}
}

Phenos=data.frame(read.table(PhenoFile,header=TRUE,stringsAsFactors=FALSE))
colnames(Phenos)=c("IID","Pheno")

PCs=data.frame(read.table(PCsFile,header=TRUE,stringsAsFactors=FALSE))
Sex=data.frame(read.table(SexFile,header=TRUE,stringsAsFactors=FALSE))
# not going to actually use sex pro tem because did not for GWAS

FormulaString="Pheno ~ PC1"
for (p in 2:20) {
	FormulaString=sprintf("%s + PC%d",FormulaString,p)
}

AllData=merge(Phenos,PCs,by="IID")

AllGWAS$ID=gsub("-",".",AllGWAS$ID)
for (s in 2:ncol(AllGenos)) {
	TestGeno=na.omit(AllGenos[,c(1,s)])
	RSName=colnames(TestGeno)[2]
	Pos=unlist(gregexpr('_', RSName))[1]
	colnames(TestGeno)[2]=substr(RSName,1,Pos-1)
	Test=merge(AllData,TestGeno,by="IID")
	M0=glm(as.formula(FormulaString),data=Test,family="binomial") # must redo each time because missing genotypes
	print(summary(M0))
	LL0=as.numeric(logLik(M0))
	M1=glm(as.formula(sprintf("%s + %s",FormulaString,colnames(TestGeno)[2])),data=Test,family="binomial")
	print(summary(M1))
	Wald_P=coef(summary(M1))[22,4]
	LL1=as.numeric(logLik(M1))
	ch2=2*(LL1-LL0)
	LR_P=pchisq(ch2,1,lower.tail=FALSE)
	print(LR_P)
	SNPNum=which(colnames(TestGeno)[2]==AllGWAS$ID)[[1]]
	AllGWAS$Wald_MLP[SNPNum]=-log10(Wald_P)
	AllGWAS$LR_MLP[SNPNum]=-log10(LR_P)
}

AllGWAS$LR_SLP=AllGWAS$LR_MLP*sign(AllGWAS$Z_STAT)
write.table(AllGWAS,ResultsFile,row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")

