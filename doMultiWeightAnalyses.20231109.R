#!/share/apps/R-3.6.1/bin/Rscript
# script to do analyses with different annotations contributing to weights
library(plyr)

geneModelsFn="/home/rejudcu/UKBB/annot.2023/geneModels.20231109.txt"
varCatsFn="/home/rejudcu/UKBB/annot.2023/varCats.20210401.txt"
catsFn="/home/rejudcu/UKBB/annot.2023/categories.txt"
extraWeightsFn="/home/rejudcu/UKBB/annot.2023/extraWeights.20210413.txt"
ver="forAnnot.20231109"

PCsFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb23155.common.all.eigenvec"
sexFile="/home/rejudcu/UKBB/UKBB.sex.20201111.txt"
vepAnnotFile="/home/rejudcu/UKBB/RAP/annot/ukb23158.AM.annot.vcf.gz"

pheno="T2D"
gene="GCK"

args = commandArgs(trailingOnly=TRUE)
if (length(args)==2) {
  pheno=args[1]
  gene=args[2]
}

scoresFn=sprintf("%s.%s.%s.sco",ver,pheno,gene)
resultsFile=sprintf("results.%s.%s.%s.txt",ver,pheno,gene)

cats=data.frame(read.table(catsFn,header=FALSE,sep="",stringsAsFactors=FALSE))
vwNames=cats[cats[,1]!="Unused",1]
extraWeights=data.frame(read.table(extraWeightsFn,header=TRUE,sep="",stringsAsFactors=FALSE))
ewNames=c("AM_prediction","AM_score",extraWeights[,1])
wNames=c(vwNames,ewNames)

if (!file.exists(scoresFn)) {
	print(sprintf("Error: file %s not found",scoresFn))
	quit()
}
	
scores=data.frame(read.table(scoresFn,header=FALSE,sep=""))
colnames(scores)[1:3]=c("IID","Pheno","VEPweight")
colnames(scores)[4:ncol(scores)]=wNames
PCsTable=data.frame(read.table(PCsFile,header=FALSE,sep="\t"))
colnames(PCsTable)[1:2]=c("FID","IID")
covars=""
for (p in 1:20) {
  covars=sprintf("%sPC%d + ",covars,p)
  colnames(PCsTable)[2+p]=sprintf("PC%d",p) # because first row of PCs file starts with hash sign FID
}

sexTable=data.frame(read.table(sexFile,header=TRUE,sep="\t"))
covars=sprintf("%s Sex ",covars)

allData=merge(scores,PCsTable,,by="IID",all=FALSE)
allData=merge(allData,sexTable,by="IID",all=FALSE)

formulaString=sprintf("Pheno ~ %s",covars)
for (c in 1:length(wNames)) {
	formulaString=sprintf("%s + %s",formulaString,wNames[c])
}

if (pheno=="BMI") {
	model=glm(as.formula(formulaString), data = allData)
} else {
	model=glm(as.formula(formulaString), data = allData, family="binomial")
}

nr <- length(wNames) + 22 - nrow(summary(model)$coefficients) 
nc <- 4
rnames <- names(which(summary(model)$aliased))
cnames <- colnames(summary(model)$coefficients)
mat_na <- matrix(data = NA,nrow = nr,ncol = nc,
           dimnames = list(rnames,cnames))
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

resultColumns=c("Weight","beta","multSLP","b","SLP")
results=data.frame(matrix(ncol=length(resultColumns),nrow=length(wNames)))
colnames(results)=resultColumns
results[,1:2]=coeffs[,1:2]
results$multSLP=sprintf("%.2f",-log10(coeffs[,5])*sign(results$beta))
for (r in 1:nrow(results)) {
	formulaString=sprintf("Pheno ~ %s + %s",covars,results[r,1])
	if (r>length(vwNames)) {
		formulaString=sprintf("%s + ProteinAltering",formulaString)
	}
	if (pheno=="BMI") {
		m=glm(as.formula(formulaString), data = allData)
	} else {
		m=glm(as.formula(formulaString), data = allData, family="binomial")
	}
	if ((nrow(summary(m)$coefficients)>22&&r<=length(vwNames)) || (nrow(summary(m)$coefficients)>23&&r>length(vwNames))) {
		results$b[r]=summary(m)$coefficients[23,1]
		results$SLP[r]=sprintf("%.2f",-log10(summary(m)$coefficients[23,4])*sign(results$b[r]))
	} else { # no valid data for this weight
		results$b[r]=0
		results$SLP[r]="0.00"
	}
}

write.table(results,resultsFile,row.names=FALSE,quote=FALSE,sep="\t")
}


