#!/share/apps/R-3.6.1/bin/Rscript

# script to get phenotype for eczema / atopic dermatitis

# note that the column number to provide is one higher than that given in http://www.davecurtis.net/UKBB/ukb41465.html
exomesFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb41465.470Kexomes.txt"
diagnosisFile="/home/rejudcu/UKBB/asthma.2024/eczema.ICD10Diagnosis.txt" # ICD10 codes for dementia
diagnosisFileICD9="/home/rejudcu/UKBB/asthma.2024/eczema.ICD10Diagnosis.txt" # ICD10 codes for dementia
ICD10SourcesFile="/home/rejudcu/UKBB/UKBBscripts.20220912/ICD10Sources.txt" # sources of ICD10 codes - hospital admissions, causes of death, etc
ICD9SourcesFile="/home/rejudcu/UKBB/UKBBscripts.20220912/ICD9Sources.txt" # sources of ICD9 codes

wd="/home/rejudcu/UKBB/asthma.2024"

selfReportedCodes=c(1452)
selfReportFirst=1730
selfReportLast=1865

# need to add 1 to these because first column is numbered 0

setwd(wd)

ICD10Sources=data.frame(read.table(ICD10SourcesFile,header=TRUE,sep="\t"))
cmd=sprintf("cut %s -f 1",exomesFile)
for (r in 1:nrow(ICD10Sources)) {
  cmd=sprintf("%s,%d-%d",cmd,ICD10Sources$First[r]+1,ICD10Sources$Last[r]+1)
}
cmd=sprintf("%s > ICD10Codes.txt",cmd)
system(cmd)

ICD10=data.frame(read.table("ICD10Codes.txt",header=TRUE,sep="\t"))

ICD9Sources=data.frame(read.table(ICD9SourcesFile,header=TRUE,sep="\t"))
cmd=sprintf("cut %s -f 1",exomesFile)
for (r in 1:nrow(ICD9Sources)) {
  cmd=sprintf("%s,%d-%d",cmd,ICD9Sources$First[r]+1,ICD9Sources$Last[r]+1)
}
cmd=sprintf("%s > ICD9Codes.txt",cmd)
system(cmd)

ICD9=data.frame(read.table("ICD9Codes.txt",header=TRUE,sep="\t"))

eczema=data.frame(matrix(ncol=5,nrow=nrow(ICD10)))
colnames(eczema)=c("IID","eczema","eczemaICD10","eczemaICD9","eczemaSelf")

diagnoses=data.frame(read.table(diagnosisFile,header=TRUE,sep="\t",quote="")) # ICD10 diagnoses
eczema$IID=ICD10$eid # I think they are all in the same order
eczema$eczemaICD10=0
for (r in 1:nrow(diagnoses)) {
  eczema$eczemaICD10[rowSums(ICD10==as.character(diagnoses[r,1]), na.rm = TRUE)>0]=1
}

diagnosesICD9=data.frame(read.table(diagnosisFileICD9,header=TRUE,sep="\t",quote="")) # ICD9 diagnoses
eczema$eczemaICD9=0
for (r in 1:nrow(diagnosesICD9)) {
  eczema$eczemaICD9[rowSums(ICD9==as.character(diagnosesICD9[r,1]), na.rm = TRUE)>0]=1
}

# get subjects by self-reported diagnosis
cmd=sprintf("tail -n +2 %s | cut -f 1,%d-%d > selfDiagnoses.txt",exomesFile,selfReportFirst+1,selfReportLast+1)
system(cmd)
selfReport=data.frame(read.table("selfDiagnoses.txt",header=FALSE,sep="\t"))
eczema$eczemaSelf=0
for (r in 1:length(selfReportedCodes)) {
  eczema$eczemaSelf[rowSums(selfReport==selfReportedCodes[r], na.rm = TRUE)>0]=1
}

eczema[,2]=0
eczema[rowSums(eczema[,3:5]==1, na.rm = TRUE)>0,2]=1
for (c in 2:ncol(eczema)) {
  toWrite=eczema[,c(1,c)]
  write.table(toWrite,sprintf("UKBB.%s.txt",colnames(eczema)[c]),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
}

write.table(eczema,"UKBB.allEczema.20241128.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)



