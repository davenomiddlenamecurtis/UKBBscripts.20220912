#!/share/apps/R-3.6.1/bin/Rscript

# script to get data files to analyse association with high alcohol consumption

# note that the column number to provide is one higher than that given in http://www.davecurtis.net/UKBB/ukb41465.html

exomesFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb41465.exomes.20201103.txt"
alcoholSourcesFile="AlcoholIntake.txt"
alcoholCodesFile="coding100006.txt"
sexCol=5
wd="/home/rejudcu/UKBB/alcohol"
# need to add 1 to these because first column is numbered 0

setwd(wd)

alcoholSources=data.frame(read.table(alcoholSourcesFile,header=TRUE,sep="\t"))

diagnosisFile="hyperlipidaemiaDiagnosis.txt" # ICD10 codes for various hyperlipidaemia diagnoses
ICD10SourcesFile="ICD10Sources.txt" # sources of ICD10 codes - hospital admissions, causes of death
cmd=sprintf("cut %s -f 1,%d",exomesFile,sexCol+1)
for (c in 2:ncol(alcoholSources)) {
  cmd=sprintf("%s,%d-%d",cmd,alcoholSources[1,c]+1,alcoholSources[2,c]+1)
}
cmd=sprintf("%s > alcoholData.txt",cmd)
system(cmd)

alcoholData=data.frame(read.table("alcoholData.txt",header=TRUE,sep="\t"))
codes=data.frame(read.table(alcoholCodesFile,header=TRUE,sep="\t"))
alcoholData=alcoholData[which(!is.na(alcoholData[,3])|!is.na(alcoholData[,4])|!is.na(alcoholData[,5])|!is.na(alcoholData[,6])|!is.na(alcoholData[,7])),]
for (r in 1:nrow(codes)) {
  alcoholData[alcoholData==codes[r,1]]=codes[r,2]
}
alcoholData[is.na(alcoholData)]=0

# just ignore the first yes/no question 
alcoholUnits=data.frame(matrix(ncol=9,nrow=nrow(alcoholData),0))
colnames(alcoholUnits)=c("IID","Sex","Inst1","Inst2","Inst3","Inst4","Inst5","MaxUnits","High")
alcoholUnits[,1:2]=alcoholData[,1:2]
for (i in 1:5) {
  for (d in 3:ncol(alcoholSources)) {
    alcoholUnits[,2+i]=alcoholUnits[,2+i]+alcoholData[7+(d-3)*5+i]*alcoholSources[3,c]
	# skip first seven columns
  }
}
alcoholUnits$MaxUnits=apply(alcoholUnits[,3:7],1,max)
head(alcoholUnits,n=100)
alcoholUnits$High[which(alcoholUnits$Sex==0 & alcoholUnits$MaxUnits>35/7)]=1
alcoholUnits$High[which(alcoholUnits$Sex==1 & alcoholUnits$MaxUnits>50/7)]=1
toWrite=alcoholUnits[c("IID","High")]
write.table(toWrite,"UKBB.alcHigh.20201229.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
toWrite=alcoholUnits[c("IID","MaxUnits")]
write.table(toWrite,"UKBB.alcUnits.20201229.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

