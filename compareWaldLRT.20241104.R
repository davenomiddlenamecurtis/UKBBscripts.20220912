#!/share/apps/R-3.6.1/bin/Rscript
# script to get data files to analyse contribution of different types of variant to phenotype


# args = commandArgs(trailingOnly=TRUE)
args=c("UKBB.psoriasis.22K.txt","rs28367584.raw")
if (length(args)!=2) {
  print("Provide phenotype and SNP files as arguments")
  quit()
}

PhenoFile=args[1]
SNPFile=args[2]
PCsFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb23158.common.all.20230806.eigenvec.txt"
SexFile="/home/rejudcu/UKBB/UKBB.sex.20230807.txt"

Phenos=data.frame(read.table(PhenoFile,header=TRUE,stringsAsFactors=FALSE))
colnames(Phenos)=c("IID","Pheno")
SNP=na.omit(data.frame(read.table(SNPFile,header=TRUE,stringsAsFactors=FALSE)))
colnames(SNP)[7]="SNP"
PCs=data.frame(read.table(PCsFile,header=TRUE,stringsAsFactors=FALSE))
Sex=data.frame(read.table(SexFile,header=TRUE,stringsAsFactors=FALSE))
# not going to actually use sex pro tem because did not for GWAS

FormulaString="Pheno ~ PC1"
for (p in 2:20) {
	FormulaString=sprintf("%s + PC%d",FormulaString,p)
}

AllData=merge(Phenos,SNP,by="IID")
AllData=merge(AllData,PCs,by="IID")

M0=glm(as.formula(FormulaString),data=AllData,family="binomial")
print(summary(M0))
LL0=as.numeric(logLik(M0))

M1=glm(as.formula(sprintf("%s + SNP",FormulaString)),data=AllData,family="binomial")
print(summary(M1))
LL1=as.numeric(logLik(M1))

ch2=2*(LL1-LL0)
p=pchisq(ch2,1,lower.tail=FALSE)
print(p)
