#!/share/apps/R-3.6.1/bin/Rscript
# script to do analyses with different annotations contributing to weights
library(plyr)

geneModelsFn="/home/rejudcu/UKBB/annot/geneModels.txt"
varCatsFn="/home/rejudcu/UKBB/annot/varCats.20210401.txt"
catsFn="/home/rejudcu/UKBB/annot/categories.txt"
extraWeightsFn="/home/rejudcu/UKBB/annot/extraWeights.20210413.txt"
ver="forAnnot.20210330"

PCsFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb23155.common.all.eigenvec"
sexFile="/home/rejudcu/UKBB/UKBB.sex.20201111.txt"
vepAnnotFile="/SAN/ugi/UGIbiobank/data/downloaded/UKBexomeOQFE.filtered_chr.annot.vcf.gz"

pheno="BMI"
gene="MC4R"

args = commandArgs(trailingOnly=TRUE)
if (length(args)==2) {
  pheno=args[1]
  gene=args[2]
}

argFn=sprintf("/home/rejudcu/UKBB/annot/gva.UKBB.%s.%s.arg",pheno,ver)
for (lw in seq(0,2,0.2)) {
wf=sprintf("%.2f",10^lw)
scoresFn=sprintf("%s.%s.%s.%s.sco",ver,pheno,gene,wf)
resultsFile=sprintf("results.%s.%s.%s.%s.txt",ver,pheno,gene,wf)

vepWeights=""
vepTest=""
cats=data.frame(read.table(catsFn,header=FALSE,sep="",stringsAsFactors=FALSE))
varCats=data.frame(read.table(varCatsFn,header=FALSE,sep="",stringsAsFactors=FALSE))
for (r in 1:nrow(cats)) {
	if (cats[r,1]=="Unused") next
	fn=sprintf("catWeight.%s.txt",cats[r,1])
	output="DEFAULT 0\n"
	for (rr in 1:nrow(varCats)) {
		if (varCats[rr,3]==cats[r,1]) output=sprintf("%s%s 1\n",output,varCats[rr,1])
	}
	writeLines(output,fn)
	vepWeights=sprintf("%s--weight-expression '( ( GETVEP (%sCSQ%s VCFLOOKUP %s%s%s)) GETWEIGHT %s%s%s )'\n",vepWeights,"\"","\"","\"",vepAnnotFile,"\"","\"",fn,"\"")
	vepWeights=sprintf("%s--weight-name %s\n",vepWeights,cats[r,1])
	vepTest=sprintf("%s%s 0.0 0 1\n",vepTest,cats[r,1])
}
writeLines(vepWeights,"gva.vepWeights.arg",sep="")
writeLines(vepTest,"vepTest.tst",sep="")
vwNames=cats[cats[,1]!="Unused",1]

extraWeights=data.frame(read.table(extraWeightsFn,header=TRUE,sep="",stringsAsFactors=FALSE))
ewTest=""
for (r in 1:nrow(extraWeights)) {
	ewTest=sprintf("%s%s 0.0 0 1\n",ewTest,extraWeights[r,1])
}
writeLines(ewTest,"ewTest.tst",sep="")
ewNames=extraWeights[,1]

commLine=sprintf("cat /home/rejudcu/UKBB/annot/coVars.tst vepTest.tst ewTest.tst > %s.tst",ver)
system(commLine)

commLine=sprintf("if [ ! -e %s -o ! -s %s ] ; then geneVarAssoc --arg-file %s --arg-file gva.vepWeights.arg --multi-weight-file %s --weight-factor %s --gene %s --keep-temp-files 1; fi",scoresFn,scoresFn,argFn,extraWeightsFn,wf,gene)
system(commLine)
commLine=sprintf("if [ ! -e %s -o ! -s %s ] ; then mv %s %s; mv %s %s; fi",
	scoresFn,scoresFn,
	sprintf("%s.%s.%s.sco",pheno,ver,gene),scoresFn,
	sprintf("%s.%s.%s.sao",pheno,ver,gene),sprintf("%s.%s.%s.%s.sao",ver,pheno,gene,wf))
system(commLine)
commLine=sprintf("rm %s.%s.%s.*",pheno,ver,gene)
system(commLine)


scores=data.frame(read.table(scoresFn,header=FALSE,sep=""))
colnames(scores)[1:3]=c("IID","Pheno","VEPweight")
wNames=c(vwNames,ewNames)
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


