#!/share/apps/R-3.6.1/bin/Rscript
# script to get data files to analyse contribution of different types of variant to phenotype
library(plyr)

pheno="HL.20201103"
gene="LDLR"

version="compAnnot.20210226"
catArgFile="gva.UKBB.varCounts.20210226.arg"
gerpArgFile="gva.UKBB.gerp.20210226.arg"
gerpLinArgFile="gva.UKBB.gerpLin.20210226.arg"

args = commandArgs(trailingOnly=TRUE)
if (length(args)==2) {
  pheno=args[1]
  gene=args[2]
}

PCsFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb23155.common.all.eigenvec"
sexFile="/home/rejudcu/UKBB/UKBB.sex.20201111.txt"

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

resultsFile=sprintf("%s.results.%s.%s.txt",version,pheno,gene)
if (file.exists(resultsFile)) { 
  quit() 
}

phenoFile=sprintf("UKBB.%s.txt",pheno)

testName=sprintf("%s.%s.cat",version,pheno)
catScoreFile=sprintf("%s.%s.sco",testName,gene)
catSaoFile=sprintf("%s.%s.sao",testName,gene)
commLine=sprintf("if [ ! -e %s ];then geneVarAssoc --arg-file %s --gene %s --test-name %s --ID-and-phenotype-file %s; fi",catScoreFile,catArgFile,gene,testName,phenoFile)
system(commLine)

testName=sprintf("%s.%s.gerp",version,pheno)
gerpScoreFile=sprintf("%s.%s.sco",testName,gene)
gerpSaoFile=sprintf("%s.%s.sao",testName,gene)
commLine=sprintf("if [ ! -e %s ];then geneVarAssoc --arg-file %s --gene %s --test-name %s --ID-and-phenotype-file %s; fi",gerpScoreFile,gerpArgFile,gene,testName,phenoFile)
system(commLine)

testName=sprintf("%s.%s.gerpLin",version,pheno)
gerpLinScoreFile=sprintf("%s.%s.sco",testName,gene)
gerpLinSaoFile=sprintf("%s.%s.sao",testName,gene)
commLine=sprintf("if [ ! -e %s ];then geneVarAssoc --arg-file %s --gene %s --test-name %s --ID-and-phenotype-file %s; fi",gerpLinScoreFile,gerpLinArgFile,gene,testName,phenoFile)
system(commLine)

lines=readLines(catSaoFile)
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

varScores=data.frame(read.table(catScoreFile,header=FALSE,sep=""))
gerpScores=data.frame(read.table(gerpScoreFile,header=FALSE,sep=""))
colnames(gerpScores)=c("IID","gerpPheno","gerp")
gerpLinScores=data.frame(read.table(gerpLinScoreFile,header=FALSE,sep=""))
colnames(gerpLinScores)=c("IID","gerpLinPheno","gerpLin")

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

allData=merge(varTypes,gerpScores,by="IID",all=FALSE)
allData=merge(allData,gerpLinScores,by="IID",all=FALSE)
typeNames=sprintf("%s + gerp + gerpLin",typeNames)
nVarTypes=nVarTypes+2
types=append(types,c("gerp","gerpLin"))

formulaString=sprintf("%s Sex %s",formulaString,typeNames)

allData=merge(allData,PCsTable,by="IID",all=FALSE)
allData=merge(allData,sexTable,by="IID",all=FALSE)

model=glm(as.formula(formulaString), data = allData, family = "binomial")
summary(model)

nr <- nVarTypes+22 - nrow(summary(model)$coefficients) 
nc <- 4
rnames <- names(which(summary(model)$aliased))
cnames <- colnames(summary(model)$coefficients)
mat_na <- matrix(data = NA,nrow = nr,ncol = nc,
           dimnames = list(rnames,cnames))
mat_coef <- rbind(summary(model)$coefficients,mat_na)

coeffs=data.frame(matrix(ncol=5,nrow=nrow(mat_coef)-22))
colnames(coeffs)[1]="Type"
coeffs[,2:5]=mat_coef[23:nrow(mat_coef),]
coeffs[,1]=rownames(mat_coef)[23:nrow(mat_coef)]
coeffs$OR=exp(coeffs[,2])
coeffs$LCL=exp(coeffs[,2]-coeffs[,3]*2)
coeffs$HCL=exp(coeffs[,2]+coeffs[,3]*2)
coeffs
# Now do a join using dply to get rows in right order
orderTypes=data.frame(matrix(ncol=1,nrow=nVarTypes))
colnames(orderTypes)[1]="Type"
orderTypes[,1]=types
coeffs=join(orderTypes,coeffs,by="Type")

resultColumns=c("Type","NumLoci","NumCont","MeanCont","NumCase","MeanCase","OR","SLP")
results=data.frame(matrix(ncol=length(resultColumns),nrow=nVarTypes))
colnames(results)=resultColumns
results$Type=types
for (v in 1:(nVarTypes-2)) {
  results$NumLoci[v]=sum(varLoci[,2+v])
  results$NumCont[v]=sum(varTypes[varTypes$Pheno==0,3+v])
  results$MeanCont[v]=sprintf("%.6f",results$NumCont[v]/sum(varTypes$Pheno==0))
  results$NumCase[v]=sum(varTypes[varTypes$Pheno==1,3+v])
  results$MeanCase[v]=sprintf("%.6f",results$NumCase[v]/sum(varTypes$Pheno==1))
}

results$OR=sprintf("%.2f (%.2f - %.2f)",coeffs$OR,coeffs$LCL,coeffs$HCL)
results$SLP=sprintf("%.2f",log10(coeffs[,5])*sign(1-coeffs$OR))

for (r in 1:(nrow(results)-2)) {
  nCont=as.numeric(results$NumCont[r])
  nCase=as.numeric(results$NumCase[r])
  if (nCont==0 || nCase==0) results$OR[r]=""
  if (nCont+nCase<50) {
    fVals=matrix(c(nCont,sum(varTypes$Pheno==0)-nCont,nCase,sum(varTypes$Pheno==1)-nCase),nrow=2,ncol=2)
	fTest=fisher.test(fVals)
  results$SLP[r]=sprintf("%.2f",log10(fTest$p.value)*sign(as.numeric(results$MeanCont[r])-as.numeric(results$MeanCase))[r])
  }
}
results

write.table(results,resultsFile,row.names=FALSE,quote=FALSE,sep="\t")


coVars=" ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + Sex"
tStats=data.frame(matrix(nrow=22,ncol=nVarTypes))
colnames(tStats)=types
for (v in 1:nVarTypes) {
  f=sprintf("%s %s",types[v],coVars)
  t=glm(as.formula(f), data = allData)
  tStats[,v]=summary(t)$coefficients[,3]
}
tStats

corrs=data.frame(matrix(nrow=21,ncol=nVarTypes))
MLPs=data.frame(matrix(nrow=21,ncol=nVarTypes))
colnames(corrs)=types
rownames(corrs)=colnames(allData)[(nVarTypes+5):(nVarTypes+25)]
colnames(MLPs)=types
rownames(MLPs)=colnames(allData)[(nVarTypes+5):(nVarTypes+25)]
for (v in 1:nVarTypes) {
  for (p in 1:21) {
    c=cor.test(allData[,v+3],allData[,nVarTypes+4+p])
    corrs[p,v]=c$estimate
    MLPs[p,v]=-log10(max(c$p.value,2.2e-16))
  }
}
corrs
MLPs

corrTableFile=sprintf("analyseVarTypes.corrs.%s.txt",gene)
MLPTableFile=sprintf("analyseVarTypes.MLPs.%s.txt",gene)
write.table(corrs,corrTableFile)
write.table(MLPs,MLPTableFile)


