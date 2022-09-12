#!/share/apps/R-3.6.1/bin/Rscript
library(plyr)

argFile="~/pars/iva.UKBB.all.20210201.arg"
varFile="/home/rejudcu/UKBB/alcProb/variants/ADH.variants.20210202.txt"
PCsFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb23155.common.all.eigenvec"
sexFile="/home/rejudcu/UKBB/UKBB.sex.20201111.txt"
gene="ADH1C"

nVarTypes=12 # 10000000000
types=c(
"IntronicEtc",
"FivePrime",
"Synonymous",
"SpliceRegion",
"ThreePrime",
"ProteinAltering",
"InDel",
"Disruptive",
"SpliceSite",
"SIFT",
"PossDam",
"ProbDam")

vars=data.frame(read.table(varFile,header=TRUE,stringsAsFactors=FALSE,sep='\t'))
intervals=""
for (r in 1:nrow(vars))
{
  intervals=sprintf("%s%d:%d-%d\n",intervals,vars$Chr[r],vars$Position[r],vars$Position[r])
}

for (m in c("alcHigh","alcProb"))
{
  wd=sprintf("~/UKBB/%s/variants",m)
  setwd(wd)
  resultsFile=sprintf("analyseVarTypes.withADH.results.%s.txt",gene)
  if (file.exists(resultsFile)) { 
    next 
  }
  writeLines(intervals,"intervals.txt")
  commLine=sprintf("intVarAssoc --arg-file %s --use-transposed-file 0 --interval-list-file intervals.txt --ID-and-phenotype-file ~/UKBB/%s/UKBB.%s.20210202.txt",argFile,m,m)
  system(commLine)  
  alls=data.frame(read.table("iva.dat",header=FALSE))
  locNames=strsplit(readLines("iva.comm.par")[1],"\\s+")[[1]]
  for (ll in 1:length(locNames))
  {
    words=strsplit(locNames[ll],":")[[1]]
	v=strsplit(words[2],"-")[[1]]
	locNames[ll]=vars$Variant[match(v[1],vars$Position)]
  }
genos=data.frame(matrix(ncol=length(locNames)+1,nrow=nrow(alls)))
colnames(genos)[1]="IID"
genos[,1]=alls[1]
for (r in 1:length(locNames))
{
  colnames(genos)[r+1]=locNames[r]
  genos[,r+1]=alls[,r*2+1]+alls[,r*2+2]-2
  # if alls were 0 0 then geno will be -2
}
genos=genos[,-3]
locNames=locNames[-2]
# drop the second allele for rs1229984
for (r in 2:(length(locNames)+1))
{
  genos[genos[,r]==-2,r]=mean(genos[genos[,r]!=-2,r])
}
# replace missing genotypes with mean score
genoFile=sprintf("ADH.genos.%s.txt",m)
write.table(genos,genoFile,row.names=FALSE, quote=FALSE)
genos=data.frame(read.table(genoFile,header=TRUE,stringsAsFactors=FALSE))

saoFile=sprintf("../genes/UKBB.%s.varCounts.20210202.%s.sao",m,gene)
varScoreFile=sprintf("../genes/UKBB.%s.varCounts.20210202.%s.sco",m,gene)
lines=readLines(saoFile)
varLoci=data.frame(matrix(ncol=2+nVarTypes,nrow=length(lines),NA))
colnames(varLoci)[1:2]=c("Locus","VarScore")
r=1
for (ll in 1:length(lines))
		{
		words=strsplit(lines[ll],"\\s+")[[1]] 
		if (length(words)<17 || !grepl(":",words[1])) next
		varLoci$Locus[r]=words[1]
		varLoci$VarScore[r]=as.numeric(words[16])
		r=r+1
		}
varLoci=varLoci[!is.na(varLoci[,1]),]

varScores=data.frame(read.table(varScoreFile,header=FALSE,sep=""))
PCsTable=data.frame(read.table(PCsFile,header=TRUE,sep="\t"))
sexTable=data.frame(read.table(sexFile,header=TRUE,sep="\t"))

colnames(PCsTable)[1:2]=c("FID","IID")
formulaString=("Pheno ~")
for (p in 1:20) {
  formulaString=sprintf("%s PC%d +",formulaString,p)
  colnames(PCsTable)[2+p]=sprintf("PC%d",p)
}

varTypes=data.frame(matrix(ncol=3+nVarTypes,nrow=nrow(varScores)))
varTypes[,1:3]=varScores
colnames(varTypes)[1:3]=c("IID","Pheno","VarScore")
typeNames=""
for (v in 1:nVarTypes) {
  colnames(varTypes)[3+v]=sprintf(types[v])
  colnames(varLoci)[2+v]=sprintf(types[v])
  typeNames=sprintf("%s + %s",typeNames,types[v])
  varTypes[,3+v]=varTypes[,3]%%10
  varTypes[,3]=floor(varTypes[,3]/10)
  varLoci[,2+v]=varLoci[,2]%%10
  varLoci[,2]=floor(varLoci[,2]/10)
}

nVar=length(locNames)
varNameString=""
for (l in 1:nVar) {
  varNameString=sprintf("%s + %s ",varNameString,locNames[l])
}

formulaString=sprintf("%s Sex %s %s",formulaString,varNameString,typeNames)

allData=merge(varTypes,PCsTable,by="IID",all=FALSE)
allData=merge(allData,sexTable,by="IID",all=FALSE)
allData=merge(allData,genos,by="IID",all=FALSE)

model=glm(as.formula(formulaString), data = allData, family = "binomial")
summary(model)

nr <- nVarTypes+nVar+22 - nrow(summary(model)$coefficients) 
nc <- 4
rnames <- names(which(summary(model)$aliased))
cnames <- colnames(summary(model)$coefficients)
mat_na <- matrix(data = NA,nrow = nr,ncol = nc,
           dimnames = list(rnames,cnames))
mat_coef <- rbind(summary(model)$coefficients,mat_na)

coeffs=data.frame(matrix(ncol=5,nrow=nrow(mat_coef)-22))
colnames(coeffs)[1]="Variant"
coeffs[,2:5]=mat_coef[23:nrow(mat_coef),]
coeffs[,1]=rownames(mat_coef)[23:nrow(mat_coef)]
coeffs$OR=exp(coeffs[,2])
coeffs$LCL=exp(coeffs[,2]-coeffs[,3]*2)
coeffs$HCL=exp(coeffs[,2]+coeffs[,3]*2)
coeffs
# Now do a join using dply to get rows in right order
orderTypes=data.frame(matrix(ncol=1,nrow=nVar+nVarTypes))
colnames(orderTypes)[1]="Variant"
orderTypes[1:nVar,1]=locNames
orderTypes[(nVar+1):(nVar+nVarTypes),1]=types
coeffs=join(orderTypes,coeffs,by="Variant")

resultColumns=c("Variant","NumLoci","NumCont","MeanCont","NumCase","MeanCase","OR","SLP")
results=data.frame(matrix(ncol=length(resultColumns),nrow=nVarTypes+nVar))
colnames(results)=resultColumns
results$Variant[1:nVar]=locNames
results$Variant[(nVar+1):(nVar+nVarTypes)]=types
for (v in 1:nVarTypes) {
  results$NumLoci[nVar+v]=sum(varLoci[,2+v])
  results$NumCont[nVar+v]=sum(varTypes[varTypes$Pheno==0,3+v])
  results$MeanCont[nVar+v]=sprintf("%.6f",results$NumCont[nVar+v]/sum(varTypes$Pheno==0))
  results$NumCase[nVar+v]=sum(varTypes[varTypes$Pheno==1,3+v])
  results$MeanCase[nVar+v]=sprintf("%.6f",results$NumCase[nVar+v]/sum(varTypes$Pheno==1))
}

results$OR=sprintf("%.2f (%.2f - %.2f)",coeffs$OR,coeffs$LCL,coeffs$HCL)
results$SLP=sprintf("%.2f",-log10(coeffs[,5])*sign(as.numeric(coeffs$OR)-1))

for (r in (nVar+1):nrow(results)) {
  nCont=as.numeric(results$NumCont[r])
  nCase=as.numeric(results$NumCase[r])
  if (nCont==0 || nCase==0) results$OR[r]=""
  if (nCont+nCase<50) {
    fVals=matrix(c(nCont,sum(varTypes$Pheno==0)-nCont,nCase,sum(varTypes$Pheno==1)-nCase),nrow=2,ncol=2)
	fTest=fisher.test(fVals)
  results$SLP[r]=sprintf("%.2f",log10(fTest$p.value)*sign(as.numeric(results$MeanCont[r])-as.numeric(results$MeanCase))[r])
  print(results$SLP[r])
  }
}
results

write.table(results,resultsFile,row.names=FALSE,quote=FALSE,sep="\t")

}


