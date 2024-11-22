#!/share/apps/R-3.6.1/bin/Rscript
# perform multivariate analysis with small group of weights.
library(plyr)

geneModelsFn="/home/rejudcu/UKBB/annot.2023/geneModels.20231109.txt"
varCatsFn="/home/rejudcu/UKBB/annot.2023/varCats.20210401.txt"
catsFn="/home/rejudcu/UKBB/annot.2023/categories.txt"
extraWeightsFn="/home/rejudcu/UKBB/annot.2023/extraWeights.20210413.txt"
ver="forAnnot.20231112"

# PCsFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb23155.common.all.eigenvec"
PCsFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb23158.common.all.20230806.eigenvec.txt"
# sexFile="/home/rejudcu/UKBB/UKBB.sex.20201111.txt"
sexFile="/home/rejudcu/UKBB/UKBB.sex.20230807.txt"
vepAnnotFile="/home/rejudcu/UKBB/RAPfiles/annot/ukb23158.AM.annot.vcf.gz"

protectiveGenes=c("NPC1L1","PCSK9","ANGPTL3","APOC3","INPPL1","DBH")

vwToUse=c("ProteinAltering","LOF")
ewToUse=c("AlphaMissense_Score",
	"MutationTaster_converted_rankscore",
	"CADD_raw_rankscore_hg19")

resultsFile=sprintf("multivariateResults.%s.txt",ver)
wNames=c(vwToUse,ewToUse)

cats=data.frame(read.table(catsFn,header=FALSE,sep="",stringsAsFactors=FALSE))
vwNames=cats[cats[,1]!="Unused",1]
extraWeights=data.frame(read.table(extraWeightsFn,header=TRUE,sep="",stringsAsFactors=FALSE))
ewNames=c("AlphaMissense_Prediction","AlphaMissense_Score",extraWeights[,1])
allWeightNames=c(vwNames,ewNames)


geneModels=data.frame(read.table(geneModelsFn,header=FALSE,sep="",stringsAsFactors=FALSE))
genes=geneModels[,2]
phenos=geneModels[,1]

results=data.frame(matrix(ncol=1+length(genes),nrow=length(wNames)+5))
colnames(results)=c("Weight",genes)
results[1:length(wNames),1]=wNames
results[length(wNames)+1,1]="LL0"
results[length(wNames)+2,1]="LLV1"
results[length(wNames)+3,1]="LL1"
results[length(wNames)+4,1]="VEPonly"
results[length(wNames)+5,1]="Overall"

PCsTable=data.frame(read.table(PCsFile,header=TRUE,sep="\t"))
head(PCsTable)
# colnames(PCsTable)[1:2]=c("FID","IID")
covars=""
for (p in 1:20) {
  covars=sprintf("%sPC%d + ",covars,p)
#   colnames(PCsTable)[2+p]=sprintf("PC%d",p) # because first row of PCs file starts with hash sign FID
}

sexTable=data.frame(read.table(sexFile,header=TRUE,sep="\t"))
colnames(sexTable)=c("IID","Sex")
head(sexTable)
covars=sprintf("%s Sex ",covars)

# for (g in 1:1)
for (g in 1:length(genes))
{
gene=genes[g]
pheno=phenos[g]
scoresFn=sprintf("UKBB.%s.%s.%s.sco",pheno,ver,gene)
scores=data.frame(read.table(scoresFn,header=FALSE,sep=""))
head(scores)
colnames(scores)[1:3]=c("IID","Pheno","VEPweight")
colnames(scores)[4:ncol(scores)]=allWeightNames
head (scores)
head(PCsTable)
allData=merge(scores,PCsTable,by="IID",all=FALSE)
head(allData)
head(sexTable)
allData=merge(allData,sexTable,by="IID",all=FALSE)
head(allData)

formulaString=sprintf("Pheno ~ %s",covars)
model=glm(as.formula(formulaString), data = allData, family="binomial")
print(model)
LL0=logLik(model)

formulaString=sprintf("%s + VEPweight",formulaString)
model=glm(as.formula(formulaString), data = allData, family="binomial")
print(model)

LLV1=logLik(model)

formulaString=sprintf("Pheno ~ %s",covars)
wNames=c(vwToUse,ewToUse)
for (c in 1:length(wNames)) {
	formulaString=sprintf("%s + %s",formulaString,wNames[c])
}

if (pheno=="BMI") {
	model=glm(as.formula(formulaString), data = allData)
} else {
	model=glm(as.formula(formulaString), data = allData, family="binomial")
}

LL1=logLik(model)

nr <- length(wNames) + 22 - nrow(summary(model)$coefficients) 
nc <- 4
rnames <- names(which(summary(model)$aliased))
cnames <- colnames(summary(model)$coefficients)
mat_na <- matrix(data = NA,nrow = nr,ncol = nc,dimnames = list(rnames,cnames))
mat_coef <- rbind(summary(model)$coefficients,mat_na)

coeffs=data.frame(matrix(ncol=5,nrow=nrow(mat_coef)-22))
colnames(coeffs)[1]="Weight"
coeffs[,2:5]=mat_coef[23:nrow(mat_coef),]
coeffs[,1]=rownames(mat_coef)[23:nrow(mat_coef)]
coeffs
# Now do a join using dply to get rows in right order
orderWeights=data.frame(matrix(ncol=1,nrow=length(wNames)))
colnames(orderWeights)[1]="Weight"
orderWeights[,1]=wNames
coeffs=join(orderWeights,coeffs,by="Weight")
results[1:length(wNames),g+1]=sprintf("%.2f",-log10(coeffs[,5])*sign(coeffs[,2]))
if (gene %in% protectiveGenes) {
	results[1:length(wNames),g+1]=sprintf("%.2f",log10(coeffs[,5])*sign(coeffs[,2]))
	} else {
	results[1:length(wNames),g+1]=sprintf("%.2f",-log10(coeffs[,5])*sign(coeffs[,2]))
	}
results[length(wNames)+1,g+1]=sprintf("%.1f",LL0)
results[length(wNames)+2,g+1]=sprintf("%.1f",LLV1)
results[length(wNames)+3,g+1]=sprintf("%.1f",LL1)
results[length(wNames)+4,g+1]=sprintf("%.2f",-log10(pchisq(2*(LLV1-LL0),df=1,lower.tail=FALSE)))
results[length(wNames)+5,g+1]=sprintf("%.2f",-log10(pchisq(2*(LL1-LL0),df=length(wNames),lower.tail=FALSE)))

write.table(results,resultsFile,row.names=FALSE,quote=FALSE,sep="\t")
}

