#!/share/apps/R-3.6.1/bin/Rscript

# script to get data files to analyse association with LOAD

# note that the column number to provide is one higher than that given in http://www.davecurtis.net/UKBB/ukb41465.html
exomesFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb41465.exomes.20201103.txt"
diagnosisFile="/home/rejudcu/UKBB/UKBBscripts/LOADDiagnosis.txt" # ICD10 codes for Alzheimer/dementia
ICD10SourcesFile="/home/rejudcu/UKBB/UKBBscripts/ICD10Sources.txt" # sources of ICD10 codes - hospital admissions, causes of death, etc
wd="/home/rejudcu/UKBB/LOAD"
# need to add 1 to these because first column is numbered 0
prescribedRxFirst=1866
prescribedRxLast=1913
dementiaRxFile="/home/rejudcu/UKBB/UKBBscripts/dementiaRx.txt" # codes for dementia treatments

setwd(wd)

ICD10Sources=data.frame(read.table(ICD10SourcesFile,header=TRUE,sep="\t"))
cmd=sprintf("cut %s -f 1",exomesFile)
for (r in 1:nrow(ICD10Sources)) {
  cmd=sprintf("%s,%d-%d",cmd,ICD10Sources$First[r]+1,ICD10Sources$Last[r]+1)
}
cmd=sprintf("%s > ICD10Codes.txt",cmd)
system(cmd)
ICD10=data.frame(read.table("ICD10Codes.txt",header=TRUE,sep="\t"))

cmd=sprintf("cut %s -f 1,%d-%d > selfReportCodes.txt",exomesFile,1730,1865)
system(cmd)
selfReport=data.frame(read.table("selfReportCodes.txt",header=TRUE,sep="\t"))
 # data coding 87 is for ICD9 - not many people have it
 
cmd=sprintf("tail -n +2 %s | cut -f 1,%d-%d > prescribedRx.txt",exomesFile,prescribedRxFirst+1,prescribedRxLast+1)
system(cmd)
prescribedRx=data.frame(read.table("prescribedRx.txt",header=FALSE,sep="\t"))
dementiaRx=data.frame(read.table(dementiaRxFile,header=FALSE,sep="\t")) # list of dementiaRx codes


LOADtable=data.frame(matrix(ncol=5,nrow=nrow(ICD10)))
colnames(LOADtable)=c("IID","LOAD","LOADICD10","SelfReport","RxIV")

diagnoses=data.frame(read.table(diagnosisFile,header=FALSE,sep="\t")) # ICD10 diagnoses
LOADtable$IID=ICD10$eid # I think they are all in the same order
LOADtable$LOADICD10=0
for (r in 1:nrow(diagnoses)) {
  LOADtable$LOADICD10[rowSums(ICD10==as.character(diagnoses[r,1]), na.rm = TRUE)>0]=1
}

LOADtable$LOADSelfReport=0
LOADtable$LOADSelfReport[rowSums(selfReport==as.character(1263), na.rm = TRUE)>0]=1
# 1263 is code for self-report of dementia/alzheimers/cognitive impairment

LOADtable$RxIV=0
for (r in 1:nrow(dementiaRx)) {
  LOADtable$RxIV[rowSums(prescribedRx==dementiaRx[r,1], na.rm = TRUE)>0]=1
}

LOADtable$LOAD=0
LOADtable$LOAD[rowSums(LOADtable[,3:5]==1, na.rm = TRUE)>0]=1

for (c in 2:ncol(LOADtable)) {
  toWrite=LOADtable[,c(1,c)]
  write.table(toWrite,sprintf("UKBB.%s.20210104.txt",colnames(LOADtable)[c]),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
}

write.table(LOADtable,"UKBB.allLOAD.20210104.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)


