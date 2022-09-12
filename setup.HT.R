#!/share/apps/R-3.6.1/bin/Rscript

# script to get data files to analyse association with hypertension

# note that the column number to provide is two higher than that given in http://www.davecurtis.net/UKBB/ukb41465.html
exomesFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb41465.exomes.txt"

diagnosisFile="hypertensionDiagnosis.txt" # ICD10 codes for various hypertension diagnoses
# http://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=19
ICD10SourcesFile="ICD10Sources.txt" # sources of ICD10 codes - hospital admissions, causes of death

# http://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=6
# Data Coding 6
# 1065	hypertension
# 1072	essential hypertension
selfReportedCodes=c(1065,1072)
selfReportFirst=1730
selfReportLast=1763

femaleRxCols=c(1472:1475)
maleRxCols=c(1584:1586)

prescribedRxFirst=1866
prescribedRxLast=1913

antihypertensivesFile="antihypertensives.txt"

# first get subjects by self-reported diagnosis
cmd=sprintf("tail -n +2 %s | cut -f 1,%d-%d > ../HT/selfDiagnoses.txt",exomesFile,selfReportFirst+2,selfReportLast+2)
system(cmd)
selfReport=data.frame(read.table("../HT/selfDiagnoses.txt",header=FALSE,sep="\t"))
hypertension=data.frame(matrix(ncol=6,nrow=nrow(selfReport)))
colnames(hypertension)=c("IID","HT","HTself","HTICD10","RxSelf","RxIV")
hypertension$IID=selfReport[,1]
hypertension$HTself=0
for (r in 1:length(selfReportedCodes)) {
  hypertension$HTself[rowSums(selfReport==selfReportedCodes[r], na.rm = TRUE)>0]=1
}

ICD10Sources=data.frame(read.table(ICD10SourcesFile,header=TRUE,sep="\t"))
cmd=sprintf("tail -n +2 %s | cut -f 1",exomesFile)
for (r in 1:nrow(ICD10Sources)) {
  cmd=sprintf("%s,%d-%d",cmd,ICD10Sources$First[r]+2,ICD10Sources$Last[r]+2)
}
cmd=sprintf("%s > ../HT/ICD10Codes.txt",cmd)
system(cmd)

ICD10=data.frame(read.table("../HT/ICD10Codes.txt",header=FALSE,sep="\t"))
hypertension$HTICD10=0

diagnoses=data.frame(read.table(diagnosisFile,header=FALSE,sep="\t")) # ICD10 diagnoses
for (r in 1:nrow(diagnoses)) {
  hypertension$HTICD10[rowSums(ICD10==as.character(diagnoses[r,2]), na.rm = TRUE)>0]=1
}

cmd=sprintf("tail -n +2 %s | cut -f 1",exomesFile)
for (r in 1:length(femaleRxCols)) {
  cmd=sprintf("%s,%d",cmd,femaleRxCols[r]+2)
}
for (r in 1:length(maleRxCols)) {
  cmd=sprintf("%s,%d",cmd,maleRxCols[r]+2)
}
cmd=sprintf("%s > ../HT/selfReportRx.txt",cmd)
system(cmd)

selfReportRx=data.frame(read.table("../HT/selfReportRx.txt",header=FALSE,sep="\t"))
hypertension$RxSelf=0
hypertension$RxSelf[rowSums(selfReportRx==2, na.rm = TRUE)>0]=1
# 2 is the code for saying they are taking blood pressure medication
# http://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=100626

cmd=sprintf("tail -n +2 %s | cut -f 1,%d-%d > ../HT/prescribedRx.txt",exomesFile,prescribedRxFirst+2,prescribedRxLast+2)
system(cmd)
prescribedRx=data.frame(read.table("../HT/prescribedRx.txt",header=FALSE,sep="\t"))
antihypertensives=data.frame(read.table(antihypertensivesFile,header=FALSE,sep="\t")) # list of antihypertensive codes
hypertension$RxIV=0
for (r in 1:nrow(antihypertensives)) {
  hypertension$RxIV[rowSums(prescribedRx==antihypertensives[r,1], na.rm = TRUE)>0]=1
}

hypertension$HT=0
hypertension$HT[rowSums(hypertension[,3:6]==1, na.rm = TRUE)>0]=1
for (c in 2:ncol(hypertension)) {
  toWrite=hypertension[,c(1,c)]
  write.table(toWrite,sprintf("../HT/UKBB.%s.txt",colnames(hypertension)[c]),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
}
write.table(hypertension,"../HT/UKBB.HTall.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)


