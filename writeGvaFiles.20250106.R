#!/share/apps/R-3.6.1/bin/Rscript
# script to do analyses with different annotations contributing to weights
library(plyr)

varCatsFn="/home/rejudcu/UKBB/annot.2023/varCats.20210401.txt"
catsFn="/home/rejudcu/UKBB/annot.2023/categories.txt"
extraWeightsFn="/home/rejudcu/UKBB/annot.2023/extraWeights.20210413.txt"
ver="forAnnot.20250106"

PCsFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb23155.common.all.eigenvec"
sexFile="/home/rejudcu/UKBB/UKBB.sex.20201111.txt"
vepAnnotFile="/home/rejudcu/UKBB/RAPfiles/annot/ukb23158.AM.annot.vcf.gz"

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
	vepWeights=sprintf("%s--weight-expression '( ( GETVEP (\"CSQ\" VCFLOOKUP \"%s\")) GETWEIGHT \"%s\" )'\n",vepWeights,vepAnnotFile,fn)
	vepWeights=sprintf("%s--weight-name %s\n",vepWeights,cats[r,1])
	vepTest=sprintf("%s%s 0.0 0 1\n",vepTest,cats[r,1])
}
vepWeights=sprintf("%s--weight-expression '( ( GETVEP (\"CSQ\" VCFLOOKUP \"%s\")) GETWEIGHT \"weights.AM.20231026.txt\" )'\n",vepWeights,vepAnnotFile)
vepWeights=sprintf("%s--weight-name AM_prediction\n",vepWeights)
vepWeights=sprintf("%s--weight-expression ' ( ( \"CSQ\" VCFLOOKUP \"%s\") GETVEPFIELD 37 ) * 10 '\n",vepWeights,vepAnnotFile)
vepWeights=sprintf("%s--weight-name AM_score\n",vepWeights)
writeLines(vepWeights,"gva.vepWeights.arg",sep="")
writeLines(vepTest,"vepTest.tst",sep="")
vwNames=cats[cats[,1]!="Unused",1]

extraWeights=data.frame(read.table(extraWeightsFn,header=TRUE,sep="",stringsAsFactors=FALSE))
ewTest="AM_prediction 0.0 0 1\nAM_score 0.0 0 1\n"
for (r in 1:nrow(extraWeights)) {
	ewTest=sprintf("%s%s 0.0 0 1\n",ewTest,extraWeights[r,1])
}
ewTest=sprintf("%sMINUSGPN 0.0 0 1\n",ewTest)
writeLines(ewTest,"ewTest.tst",sep="")
