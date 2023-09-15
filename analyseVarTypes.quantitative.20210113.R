#!/share/apps/R-3.6.1/bin/Rscript

# script to get data files to analyse contribution of different types of variant to phenotype
library(plyr)
scoreFileTemplate="UKBB.HL.varCounts.20201211.%s.sco"
argFile="~/pars/gva.UKBB.HL.varCounts.20201211.arg"
wd="~/UKBB/lipids/HL.20201103/genes"
scoreFileTemplate="UKBB.alcHigh.varCounts.20201231.%s.sco"
argFile="~/pars/gva.UKBB.alcHigh.varCounts.20201231.arg"
wd="~/UKBB/alcohol/genes"

scoreFileTemplate="UKBB.T2D.varCounts.20210106.%s.sco"
argFile="~/pars/gva.UKBB.T2D.varCounts.20210106.arg"
wd="~/UKBB/T2D.20210104/genes"
geneListFile="genesToRedo.txt"

scoreFileTemplate="UKBB.alcProb.varCounts.20210111.%s.sco"
argFile="~/pars/gva.UKBB.alcProb.varCounts.20210111.arg"
wd="~/UKBB/alcoholProblems/genes"
geneListFile="allGenes.txt"

test="UKBB.BMI.varCounts.20210113"
scoreFileTemplate="UKBB.BMI.varCounts.20210113.%s.sco"
argFile="~/pars/gva.UKBB.BMI.varCounts.20210113.arg"
wd="/home/rejudcu/UKBB/BMI.20210111/genes"
# geneListFile="allGenes.BMI.20210113.txt"
geneListFile="genesFrom640K.txt"

test="UKBB.BMI.varCounts.forHTR2C.20230102"
scoreFileTemplate="UKBB.BMI.varCounts.forHTR2C.20230102.%s.sco"
argFile="gva.UKBB.BMI.varCounts.forHTR2C.20230102.arg"
wd="/home/rejudcu/UKBB/HTR2C"


genes=c(
"LDLR",
"PCSK9",
"G6PC",
"GCK",
"LIPG",
"ABCA1",
"CILP2",
"PPP1R3G",
"NPC1L1",
"ANGPTL3",
"ANGPTL4",
"DHCR24",
"LPL",
"APOC3",
"CETP",
"HMGCR",
"APOB",
"STAP1")

genes=c(
"LDLR",
"PCSK9",
"HMGCR",
"APOB",
"STAP1")

genes=c(
"HTR2C"
)


setwd(wd)

if (file.exists(geneListFile)) {
 genes=data.frame(read.table(geneListFile,header=FALSE,stringsAsFactors=FALSE))[,1]
}

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=0) {
  genes=args
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

for (gene in genes) {
print("on line 87")
print(gene)
varScoreFile=sprintf("%s.%s.sco",test,gene)
saoFile=sprintf("%s.%s.sao",test,gene)
resultsFile=sprintf("analyseVarTypes.results.%s.txt",gene)
if (file.exists(resultsFile)) { 
  next 
}

print(varScoreFile)
if (!file.exists(varScoreFile)) {
	commLine=sprintf("geneVarAssoc --arg-file %s --gene %s",argFile, gene)
	system(commLine)
}
print("on line 100")

lines=readLines(saoFile)
varLoci=data.frame(matrix(ncol=2+nVarTypes,nrow=length(lines),NA))
colnames(varLoci)[1:2]=c("Locus","VarScore")
r=1
for (ll in 1:length(lines))
		{
		words=strsplit(lines[ll],"\\s+")[[1]] 
		if (length(words)<15 || !grepl(":",words[1])) next
		varLoci$Locus[r]=words[1]
		varLoci$VarScore[r]=as.numeric(words[14])
		r=r+1
		}
varLoci=varLoci[!is.na(varLoci[,1]),]

varScores=data.frame(read.table(varScoreFile,header=FALSE,sep=""))
PCsTable=data.frame(read.table(PCsFile,header=TRUE,sep="\t"))
sexTable=data.frame(read.table(sexFile,header=TRUE,sep="\t"))

print("on line 120")

colnames(PCsTable)[1:2]=c("FID","IID")
formulaString=("Pheno ~")
for (p in 1:20) {
  formulaString=sprintf("%s PC%d +",formulaString,p)
  colnames(PCsTable)[2+p]=sprintf("PC%d",p)
}

varTypes=data.frame(matrix(ncol=3+nVarTypes,nrow=nrow(varScores)))
varTypes[,1:3]=varScores
colnames(varTypes)[1:3]=c("IID","Pheno","VarScore")
varTypes$None=0
varTypes$None[varTypes$VarScore==0]=1 # extra column
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
formulaString=sprintf("%s Sex %s",formulaString,typeNames)

allData=merge(varTypes,PCsTable,by="IID",all=FALSE)
allData=merge(allData,sexTable,by="IID",all=FALSE)

model=glm(as.formula(formulaString), data = allData)
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
for (v in 1:(nVarTypes+1)) {
  if (v<=nVarTypes){
    results$NumLoci[v]=sum(varLoci[,2+v])
  }
  results$TotLoad[v]=sum(varTypes[,3+v])
  results$AveLoad[v]=sprintf("%.6f",results$TotLoad[v]/nrow(varTypes))
  results$CarrierMean[v]=mean(varTypes$Pheno[varTypes[,3+v]!=0])
  results$CarrierSD[v]=sd(varTypes$Pheno[varTypes[,3+v]!=0])
}

results$Effect[1:nVarTypes]=sprintf("%.2f (%.2f - %.2f)",coeffs$beta,coeffs$LCL,coeffs$HCL)
results$SLP[1:nVarTypes]=sprintf("%.2f",-log10(coeffs[,5])*sign(coeffs$beta))

write.table(results,resultsFile,row.names=FALSE,quote=FALSE,sep="\t")

coVars=" ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + Sex"
print("on line 195")

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
print("on line 217")
print(gene)

}

