#!/share/apps/R-3.6.1/bin/Rscript

# script to get data files to analyse association with hypercholesterolaemia

# note that the column number to provide is one higher than that given in http://www.davecurtis.net/UKBB/ukb41465.html

exomesFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb41465.470Kexomes.txt"

diagnosisFile="hyperlipidaemiaDiagnosis.txt" # ICD10 codes for various hyperlipidaemia diagnoses
ICD10SourcesFile="ICD10Sources.txt" # sources of ICD10 codes - hospital admissions, causes of death

selfReportedCodes=c(1473)
# Data Coding 6
# 28	1473	high cholesterol
selfReportFirst=1730
selfReportLast=1763

sexCol=5

femaleRxCols=c(1472:1475)
maleRxCols=c(1584:1586)

prescribedRxFirst=1866
prescribedRxLast=1913

statinsFile="statins.txt"

# first get subjects by self-reported diagnosis
cmd=sprintf("tail -n +2 %s | cut -f 1,%d-%d > ../HL.2023/selfDiagnoses.txt",exomesFile,selfReportFirst+1,selfReportLast+1)
system(cmd)
selfReport=data.frame(read.table("../HL.2023/selfDiagnoses.txt",header=FALSE,sep="\t"))
hyperlipidaemia=data.frame(matrix(ncol=6,nrow=nrow(selfReport)))
colnames(hyperlipidaemia)=c("IID","HL","HCself","HLICD10","RxSelf","RxIV")
hyperlipidaemia$IID=selfReport[,1]
hyperlipidaemia$HCself=0
for (r in 1:length(selfReportedCodes)) {
  hyperlipidaemia$HCself[rowSums(selfReport==selfReportedCodes[r], na.rm = TRUE)>0]=1
}

ICD10Sources=data.frame(read.table(ICD10SourcesFile,header=TRUE,sep="\t"))
cmd=sprintf("tail -n +2 %s | cut -f 1",exomesFile)
for (r in 1:nrow(ICD10Sources)) {
  cmd=sprintf("%s,%d-%d",cmd,ICD10Sources$First[r]+1,ICD10Sources$Last[r]+1)
}
cmd=sprintf("%s > ../HL.2023/ICD10Codes.txt",cmd)
system(cmd)

ICD10=data.frame(read.table("../HL.2023/ICD10Codes.txt",header=FALSE,sep="\t"))
hyperlipidaemia$HLICD10=0

diagnoses=data.frame(read.table(diagnosisFile,header=FALSE,sep="\t")) # ICD10 diagnoses
for (r in 1:nrow(diagnoses)) {
  hyperlipidaemia$HLICD10[rowSums(ICD10==as.character(diagnoses[r,2]), na.rm = TRUE)>0]=1
}


cmd=sprintf("tail -n +2 %s | cut -f 1",exomesFile)
for (r in 1:length(femaleRxCols)) {
  cmd=sprintf("%s,%d",cmd,femaleRxCols[r]+1)
}
for (r in 1:length(maleRxCols)) {
  cmd=sprintf("%s,%d",cmd,maleRxCols[r]+1)
}
cmd=sprintf("%s > ../HL.2023/selfReportRx.txt",cmd)
system(cmd)

selfReportRx=data.frame(read.table("../HL.2023/selfReportRx.txt",header=FALSE,sep="\t"))
hyperlipidaemia$RxSelf=0
hyperlipidaemia$RxSelf[rowSums(selfReportRx==1, na.rm = TRUE)>0]=1
# 1 is the code for saying they are taking cholesterol-lowering medication

cmd=sprintf("tail -n +2 %s | cut -f 1,%d-%d > ../HL.2023/prescribedRx.txt",exomesFile,prescribedRxFirst+1,prescribedRxLast+1)
system(cmd)
prescribedRx=data.frame(read.table("../HL.2023/prescribedRx.txt",header=FALSE,sep="\t"))
statins=data.frame(read.table(statinsFile,header=FALSE,sep="")) # list of statin codes
hyperlipidaemia$RxIV=0
for (r in 1:nrow(statins)) {
  hyperlipidaemia$RxIV[rowSums(prescribedRx==statins[r,2], na.rm = TRUE)>0]=1
}

hyperlipidaemia$HL=0
hyperlipidaemia$HL[rowSums(hyperlipidaemia[,3:6]==1, na.rm = TRUE)>0]=1
for (c in 2:ncol(hyperlipidaemia)) {
  toWrite=hyperlipidaemia[,c(1,c)]
  write.table(toWrite,sprintf("../HL.2023/UKBB.%s.20230915.txt",colnames(hyperlipidaemia)[c]),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
}

write.table(hyperlipidaemia,"../HL.2023/UKBB.HLall.20230915.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

phenos=as.matrix(hyperlipidaemia[,2:5])
print(cor(phenos))



