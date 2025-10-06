#!/share/apps/R-3.6.1/bin/Rscript
# script to get data files to analyse contribution of different types of variant to phenotype

library(plyr)

test="UKBB.BMI.separateVarCounts.20240405"
wd="/home/rejudcu/UKBB/BMI.20240329/genes"
geneListFile="/home/rejudcu/UKBB/BMI.20240329/BMI.confirmed.genes.20240403.txt"


setwd(wd)

if (file.exists(geneListFile)) {
 genes=data.frame(read.table(geneListFile,header=FALSE,stringsAsFactors=FALSE))[,1]
}

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=0) {
  genes=args
}

PCsFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb23158.common.all.20230806.eigenvec.txt"
sexFile="/home/rejudcu/UKBB/UKBB.sex.20230807.txt"

nVarTypes=11
types=c(
"IntronicEtc",
"FivePrime",
"Synonymous",
"SpliceRegion",
"ThreePrime",
"ProteinAltering",
"InDel",
"LOF",
"SIFT",
"PossDam",
"ProbDam")
# as in gva.UKBB.separateVarCounts.20230825.arg

for (sex in c("male","female")) {

for (gene in genes) {
print(gene)
varScoreFile=sprintf("%s.%s.sco",test,gene)
saoFile=sprintf("%s.%s.sao",test,gene)
resultsFile=sprintf("analyseVarTypes.results.%s.%s.txt",gene,sex)
if (file.exists(resultsFile)) { 
  next 
}

print(varScoreFile)
if (!file.exists(varScoreFile)) {
	print(sprintf("Score file %s not found",varScoreFile))
	quit()
}

lines=readLines(saoFile)
numLoci=rep(0,nVarTypes)
for (ll in 1:length(lines))
		{
		words=strsplit(lines[ll],"\\s+")[[1]] 
		if (length(words)<15 || !grepl(":",words[1])) next
		numLoci=numLoci+as.numeric(words[14:(13+nVarTypes)]) 
		}

varScores=data.frame(read.table(varScoreFile,header=FALSE,sep=""))
PCsTable=data.frame(read.table(PCsFile,header=TRUE,sep="\t"))
sexTable=data.frame(read.table(sexFile,header=TRUE,sep="\t"))
colnames(sexTable)[2]="Sex"
if (sex=="male") {
	sexTable=sexTable[sexTable$Sex==0,]
} else {
	sexTable=sexTable[sexTable$Sex==1,]
}

colnames(PCsTable)[1:2]=c("FID","IID")
formulaString=("Pheno ~")
for (p in 1:20) {
  formulaString=sprintf("%s PC%d +",formulaString,p)
  colnames(PCsTable)[2+p]=sprintf("PC%d",p)
}

varTypes=data.frame(read.table(varScoreFile,header=FALSE,sep=""))
colnames(varTypes)=c("IID","Pheno",types)
varTypes$None=0
varTypes$None[rowSums(varTypes[,3:(nVarTypes+2)])==0]=1
typeNames=""
for (v in 1:nVarTypes) {
  typeNames=sprintf("%s + %s",typeNames,types[v])
}
formulaString=sprintf("%s %s",formulaString,typeNames)

allData=merge(varTypes,PCsTable,by="IID",all=FALSE)
allData=merge(allData,sexTable,by="IID",all=FALSE)

model=glm(as.formula(formulaString), data = allData)
summary(model)

nr <- nVarTypes+21 - nrow(summary(model)$coefficients) 
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
coeffs$beta=coeffs[,2]
coeffs$LCL=coeffs[,2]-coeffs[,3]*2
coeffs$HCL=coeffs[,2]+coeffs[,3]*2
coeffs
# Now do a join using dply to get rows in right order
orderTypes=data.frame(matrix(ncol=1,nrow=nVarTypes))
colnames(orderTypes)[1]="Type"
orderTypes[,1]=types
coeffs=join(orderTypes,coeffs,by="Type")

resultColumns=c("Type","NumLoci","TotLoad","AveLoad","CarrierMean","CarrierSD","SLP","Effect")
results=data.frame(matrix(ncol=length(resultColumns),nrow=nVarTypes+1))
colnames(results)=resultColumns
results$Type[1:nVarTypes]=types
results$Type[nVarTypes+1]="None"
results$NumLoci[1:nVarTypes]=numLoci
for (v in 1:(nVarTypes+1)) {
  results$TotLoad[v]=sum(varTypes[,2+v])
  results$AveLoad[v]=sprintf("%.6f",results$TotLoad[v]/nrow(varTypes))
  results$CarrierMean[v]=mean(varTypes$Pheno[varTypes[,2+v]!=0])
  results$CarrierSD[v]=sd(varTypes$Pheno[varTypes[,2+v]!=0])
}

results$Effect[1:nVarTypes]=sprintf("%.2f (%.2f - %.2f)",coeffs$beta,coeffs$LCL,coeffs$HCL)
results$SLP[1:nVarTypes]=sprintf("%.2f",-log10(coeffs[,5])*sign(coeffs$beta))

write.table(results,resultsFile,row.names=FALSE,quote=FALSE,sep="\t")

coVars=" ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20"

corrs=data.frame(matrix(nrow=21,ncol=nVarTypes+1))
MLPs=data.frame(matrix(nrow=21,ncol=nVarTypes+1))
colnames(corrs)=c("Covar",types)
rownames(corrs)=colnames(allData)[(nVarTypes+5):(nVarTypes+25)]
corrs[,1]=rownames(corrs)
colnames(MLPs)=c("Covar",types)
rownames(MLPs)=colnames(allData)[(nVarTypes+5):(nVarTypes+25)]
MLPs[,1]=rownames(MLPs)
for (v in 1:nVarTypes) {
  for (p in 1:21) {
    c=cor.test(allData[,v+2],allData[,nVarTypes+4+p])
    corrs[p,v+1]=c$estimate
    MLPs[p,v+1]=-log10(max(c$p.value,2.2e-16))
  }
}
corrs
MLPs

# corrTableFile=sprintf("analyseVarTypes.corrs.%s.txt",gene)
# MLPTableFile=sprintf("analyseVarTypes.MLPs.%s.txt",gene)
# write.table(corrs,corrTableFile,row.names=FALSE,quote=FALSE)
# write.table(MLPs,MLPTableFile,row.names=FALSE,quote=FALSE)
print(gene)

}

}

