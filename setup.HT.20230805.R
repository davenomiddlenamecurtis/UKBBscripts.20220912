#!/share/apps/R-3.6.1/bin/Rscript

# script to get data files to analyse association with hypertension

# note that the column number to provide is one higher than that given in http://www.davecurtis.net/UKBB/ukb41465.html
exomesFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb43357.txt"
# the raw UK Biobank file has a leading tab, FFS
# so add 2 to all the column numbers below

diagnosisFile="hypertensionDiagnosis.txt" # ICD10 codes for various hypertension diagnoses
# http://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=19
ICD10SourcesFile="ICD10Sources.20230811.txt" # sources of ICD10 codes - hospital admissions, causes of death
wd="/home/rejudcu/UKBB/UKBBscripts.20220912"

# http://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=6
# Data Coding 6
# 1065	hypertension
# 1072	essential hypertension
selfReportedCodes=c(1065,1072)
selfReportFirst=5634	
selfReportLast=5769

femaleRxCols=c(5325:5340)
maleRxCols=c(5437:5448)

prescribedRxFirst=5770
prescribedRxLast=5961

setwd(wd)

antihypertensivesFile="antihypertensives.txt"

# first get subjects by self-reported diagnosis
cmd=sprintf("tail -n +2 %s | cut -f 2,%d-%d > ../HT.20230805/selfDiagnoses.txt",exomesFile,selfReportFirst+2,selfReportLast+2)
system(cmd)
selfReport=data.frame(read.table("../HT.20230805/selfDiagnoses.txt",header=FALSE,sep="\t"))
hypertension=data.frame(matrix(ncol=6,nrow=nrow(selfReport)))
colnames(hypertension)=c("IID","HT","HTself","HTICD10","RxSelf","RxIV")
hypertension$IID=selfReport[,1]
hypertension$HTself=0
for (r in 1:length(selfReportedCodes)) {
  hypertension$HTself[rowSums(selfReport==selfReportedCodes[r], na.rm = TRUE)>0]=1
}

ICD10Sources=data.frame(read.table(ICD10SourcesFile,header=TRUE,sep="\t"))
cmd=sprintf("tail -n +2 %s | cut -f 2",exomesFile)
for (r in 1:nrow(ICD10Sources)) {
  cmd=sprintf("%s,%d-%d",cmd,ICD10Sources$First[r]+2,ICD10Sources$Last[r]+2)
}
cmd=sprintf("%s > ../HT.20230805/ICD10Codes.txt",cmd)
system(cmd)

ICD10=data.frame(read.table("../HT.20230805/ICD10Codes.txt",header=FALSE,sep="\t"))
hypertension$HTICD10=0

diagnoses=data.frame(read.table(diagnosisFile,header=FALSE,sep="\t")) # ICD10 diagnoses
for (r in 1:nrow(diagnoses)) {
  hypertension$HTICD10[rowSums(ICD10==as.character(diagnoses[r,2]), na.rm = TRUE)>0]=1
}

cmd=sprintf("tail -n +2 %s | cut -f 2",exomesFile)
for (r in 1:length(femaleRxCols)) {
  cmd=sprintf("%s,%d",cmd,femaleRxCols[r]+2)
}
for (r in 1:length(maleRxCols)) {
  cmd=sprintf("%s,%d",cmd,maleRxCols[r]+2)
}
cmd=sprintf("%s > ../HT.20230805/selfReportRx.txt",cmd)
system(cmd)

selfReportRx=data.frame(read.table("../HT.20230805/selfReportRx.txt",header=FALSE,sep="\t"))
hypertension$RxSelf=0
hypertension$RxSelf[rowSums(selfReportRx==2, na.rm = TRUE)>0]=1
# 2 is the code for saying they are taking blood pressure medication
# http://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=100626

cmd=sprintf("tail -n +2 %s | cut -f 2,%d-%d > ../HT.20230805/prescribedRx.txt",exomesFile,prescribedRxFirst+2,prescribedRxLast+2)
system(cmd)
prescribedRx=data.frame(read.table("../HT.20230805/prescribedRx.txt",header=FALSE,sep="\t"))
antihypertensives=data.frame(read.table(antihypertensivesFile,header=FALSE,sep="\t")) # list of antihypertensive codes
hypertension$RxIV=0
for (r in 1:nrow(antihypertensives)) {
  hypertension$RxIV[rowSums(prescribedRx==antihypertensives[r,1], na.rm = TRUE)>0]=1
}

hypertension$HT=0
hypertension$HT[rowSums(hypertension[,3:6]==1, na.rm = TRUE)>0]=1
for (c in 2:ncol(hypertension)) {
  toWrite=hypertension[,c(1,c)]
  write.table(toWrite,sprintf("../HT.20230805/UKBB.%s.txt",colnames(hypertension)[c]),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
}
write.table(hypertension,"../HT.20230805/UKBB.HTall.20230805.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

phenos=as.matrix(hypertension[,2:6])
print(cor(phenos))