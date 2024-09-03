#!/share/apps/R-3.6.1/bin/Rscript
# script to get data files to analyse contribution of different types of variant to phenotype
library(plyr)

gene="ATP1A2"
model="UKBB.migraine.separateVarCounts.20211015"
model="UKBB.psoriasis.separateVarCounts.20230329" 
gene="NPC1L1"
args=c("UKBB.T2D.separateVarCounts.20231006","TCF7L2")
args = commandArgs(trailingOnly=TRUE)
# args=c("UKBB.HT.separateVarCounts.20230825","DNMT3A")
if (length(args)!=2) {
  print("Provide model and gene as arguments")
  quit()
}
model=args[1]
gene=args[2]


PCsFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb23158.common.all.20230806.eigenvec.txt"
sexFile="/home/rejudcu/UKBB/UKBB.sex.20230807.txt"
# covariates for all 470K exomes

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

varScoreFile=sprintf("%s.%s.sco",model,gene)
saoFile=sprintf("%s.%s.sao",model,gene)
resultsFile=sprintf("analyseVarTypes.results.%s.txt",gene)
if (file.exists(resultsFile)) { 
  print(sprintf("Results file %s already exists",resultsFile))
  quit()
}

if (!file.exists(varScoreFile)) {
	print(sprintf("Cannot find score file %s",varScoreFile))
	quit()
}

lines=readLines(saoFile)
numLoci=rep(0,nVarTypes)
for (ll in 1:length(lines))
		{
		words=strsplit(lines[ll],"\\s+")[[1]] 
		if (length(words)<17 || !grepl(":",words[1])) next
		numLoci=numLoci+as.numeric(words[16:(15+nVarTypes)])
		}

PCsTable=data.frame(read.table(PCsFile,header=TRUE,sep="\t"))
sexTable=data.frame(read.table(sexFile,header=TRUE,sep="\t"))
colnames(sexTable)=c("IID","Sex")

colnames(PCsTable)[1:2]=c("FID","IID")
formulaString=("Pheno ~")
for (p in 1:20) {
  formulaString=sprintf("%s PC%d +",formulaString,p)
  colnames(PCsTable)[2+p]=sprintf("PC%d",p)
}

varTypes=data.frame(read.table(varScoreFile,header=FALSE,sep=""))
colnames(varTypes)=c("IID","Pheno",types)
typeNames=""
for (v in 1:nVarTypes) {
  typeNames=sprintf("%s + %s",typeNames,types[v])
}
formulaString=sprintf("%s Sex %s",formulaString,typeNames)

allData=merge(varTypes,PCsTable,by="IID",all=FALSE)
allData=merge(allData,sexTable,by="IID",all=FALSE)

m=glm(as.formula(formulaString), data = allData, family = "binomial")
summary(m)

nr <- nVarTypes+22 - nrow(summary(m)$coefficients) 
nc <- 4
rnames <- names(which(summary(m)$aliased))
cnames <- colnames(summary(m)$coefficients)
mat_na <- matrix(data = NA,nrow = nr,ncol = nc,
           dimnames = list(rnames,cnames))
mat_coef <- rbind(summary(m)$coefficients,mat_na)

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
results$NumLoci=numLoci
for (v in 1:nVarTypes) {
  results$NumCont[v]=sum(varTypes[varTypes$Pheno==0,2+v])
  results$MeanCont[v]=sprintf("%.6f",results$NumCont[v]/sum(varTypes$Pheno==0))
  results$NumCase[v]=sum(varTypes[varTypes$Pheno==1,2+v])
  results$MeanCase[v]=sprintf("%.6f",results$NumCase[v]/sum(varTypes$Pheno==1))
}

results$OR=sprintf("%.2f (%.2f - %.2f)",coeffs$OR,coeffs$LCL,coeffs$HCL)
results$SLP=sprintf("%.2f",log10(coeffs[,5])*sign(1-coeffs$OR))

for (r in 1:nrow(results)) {
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
rownames(corrs)=colnames(allData)[(nVarTypes+4):(nVarTypes+24)]
colnames(MLPs)=types
rownames(MLPs)=colnames(allData)[(nVarTypes+4):(nVarTypes+24)]
for (v in 1:nVarTypes) {
  for (p in 1:21) {
    c=cor.test(allData[,v+2],allData[,nVarTypes+2+p])
    corrs[p,v]=c$estimate
    MLPs[p,v]=-log10(max(c$p.value,2.2e-16))
  }
}
corrs
MLPs

corrTableFile=sprintf("analyseVarTypes.corrs.%s.txt",gene)
MLPTableFile=sprintf("analyseVarTypes.MLPs.%s.txt",gene)
write.table(corrs,corrTableFile,quote=FALSE,sep="\t")
write.table(MLPs,MLPTableFile,quote=FALSE,sep="\t")


