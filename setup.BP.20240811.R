#!/share/apps/R-3.6.1/bin/Rscript

# script to get data files to analyse association with bipolar disorder

# note that the column number to provide is one higher than that given in http://www.davecurtis.net/UKBB/ukb43357.html

exomesFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb43357.470Kexomes.txt"

diagnosisFile="bipolarDiagnosis.txt" # ICD10 codes for various bipolar diagnoses
# http://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=19
ICD10SourcesFile="ICD10Sources.20230811.txt" # sources of ICD10 codes - hospital admissions, causes of death
wd="/home/rejudcu/UKBB/UKBBscripts.20220912"

# http://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=6
# Data Coding 6
# 1291	mania/bipolar disorder/manic depression
selfReportedCodes=c(1291)
selfReportFirst=5634	
selfReportLast=5769

prescribedRxFirst=5770
prescribedRxLast=5961

setwd(wd)

bipolarRxFile="bipolarRx.txt"

# first get subjects by self-reported diagnosis
# cmd=sprintf("tail -n +2 %s | cut -f 1,%d-%d > ../HT.20230805/selfDiagnoses.txt",exomesFile,selfReportFirst+1,selfReportLast+1)
# system(cmd)
selfReport=data.frame(read.table("../HT.20230805/selfDiagnoses.txt",header=FALSE,sep="\t"))
bipolar=data.frame(matrix(ncol=5,nrow=nrow(selfReport)))
colnames(bipolar)=c("IID","BP.20240811","BPself","BPICD10","RxIV")
bipolar$IID=selfReport[,1]
bipolar$BPself=0
for (r in 1:length(selfReportedCodes)) {
  bipolar$BPself[rowSums(selfReport==selfReportedCodes[r], na.rm = TRUE)>0]=1
}

ICD10Sources=data.frame(read.table(ICD10SourcesFile,header=TRUE,sep="\t"))
cmd=sprintf("tail -n +2 %s | cut -f 1",exomesFile)
for (r in 1:nrow(ICD10Sources)) {
  cmd=sprintf("%s,%d-%d",cmd,ICD10Sources$First[r]+1,ICD10Sources$Last[r]+1)
}
cmd=sprintf("%s > ../HT.20230805/ICD10Codes.txt",cmd)
# system(cmd)

ICD10=data.frame(read.table("../HT.20230805/ICD10Codes.txt",header=FALSE,sep="\t"))
bipolar$BPICD10=0

diagnoses=data.frame(read.table(diagnosisFile,header=FALSE,sep="\t")) # ICD10 diagnoses
for (r in 1:nrow(diagnoses)) {
  bipolar$BPICD10[rowSums(ICD10==as.character(diagnoses[r,1]), na.rm = TRUE)>0]=1
}

cmd=sprintf("tail -n +2 %s | cut -f 1,%d-%d > ../HT.20230805/prescribedRx.txt",exomesFile,prescribedRxFirst+1,prescribedRxLast+1)
# system(cmd)
prescribedRx=data.frame(read.table("../HT.20230805/prescribedRx.txt",header=FALSE,sep="\t"))
bipolarRx=data.frame(read.table(bipolarRxFile,header=FALSE,sep="\t")) # list of lithium codes
bipolar$RxIV=0
for (r in 1:nrow(bipolarRx)) {
  bipolar$RxIV[rowSums(prescribedRx==bipolarRx[r,1], na.rm = TRUE)>0]=1
}

bipolar[,2]=0
bipolar[rowSums(bipolar[,3:5]==1, na.rm = TRUE)>0,2]=1
for (c in 2:ncol(bipolar)) {
  toWrite=bipolar[,c(1,c)]
  write.table(toWrite,sprintf("../BP.20240811/UKBB.%s.txt",colnames(bipolar)[c]),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
}
write.table(bipolar,"../BP.20240811/UKBB.BPall.20240811.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

phenos=as.matrix(bipolar[,2:5])
print(cor(phenos))