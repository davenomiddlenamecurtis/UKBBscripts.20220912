#!/share/apps/R-3.6.1/bin/Rscript

# script to get data files to analyse association with FAMDEM

exomesFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb41465.470Kexomes.txt"
ParentalIllnessFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb43357.470Kexomes.txt"
diagnosisFile="/home/rejudcu/UKBB/dementia.2024/DEMDiagnosis.txt" # ICD10 codes for dementia
diagnosisFileICD9="/home/rejudcu/UKBB/dementia.2024/DEMDiagnosisICD9.txt" # ICD9 codes for dementia
ICD10SourcesFile="/home/rejudcu/UKBB/UKBBscripts.20220912/ICD10Sources.txt" # sources of ICD10 codes - hospital admissions, causes of death, etc
ICD9SourcesFile="/home/rejudcu/UKBB/UKBBscripts.20220912/ICD9Sources.txt" # sources of ICD9 codes
wd="/home/rejudcu/UKBB/dementia.2024"
FatherIllnessCols=c(7002,7041)
MotherIllnessCols=c(7042,7085)
SiblingIllnessCols=c(7086,7133)
# 10 is Alzheimer/dementia

selfReportedCodes=c(1263)
selfReportFirst=1730
selfReportLast=1865
# for dementia, this will include reports at follow-up visits

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

DEM=data.frame(matrix(ncol=9,nrow=nrow(ICD10)))
colnames(DEM)=c("IID","DEMICD10","DEMICD9","DEMFather","DEMMother","DEMSibling","DEMSelf","DEMScore","DEMMax2")

diagnoses=data.frame(read.table(diagnosisFile,header=TRUE,sep="\t",quote="")) # ICD10 diagnoses
DEM$IID=ICD10$eid # I think they are all in the same order
DEM$DEMICD10=0
for (r in 1:nrow(diagnoses)) {
  DEM$DEMICD10[rowSums(ICD10==as.character(diagnoses[r,1]), na.rm = TRUE)>0]=1
}

diagnosesICD9=data.frame(read.table(diagnosisFileICD9,header=TRUE,sep="\t",quote="")) # ICD9 diagnoses
DEM$DEMICD9=0
for (r in 1:nrow(diagnosesICD9)) {
  DEM$DEMICD9[rowSums(ICD9==as.character(diagnosesICD9[r,1]), na.rm = TRUE)>0]=1
}

DEM$DEMFather=0
cmd=sprintf("cut %s -f 1,%d-%d > %s",ParentalIllnessFile,FatherIllnessCols[1]+1,FatherIllnessCols[2]+1,"FatherIllness.txt")
system(cmd)
FatherIllness=data.frame(read.table("FatherIllness.txt",header=TRUE,sep="\t"))
DEM$DEMFather[rowSums(FatherIllness==10, na.rm = TRUE)>0]=1 # 10 is code for dementia

DEM$DEMMother=0
cmd=sprintf("cut %s -f 1,%d-%d > %s",ParentalIllnessFile,MotherIllnessCols[1]+1,MotherIllnessCols[2]+1,"MotherIllness.txt")
system(cmd)
MotherIllness=data.frame(read.table("MotherIllness.txt",header=TRUE,sep="\t"))
DEM$DEMMother[rowSums(MotherIllness==10, na.rm = TRUE)>0]=1 # 10 is code for dementia

DEM$DEMSibling=0
cmd=sprintf("cut %s -f 1,%d-%d > %s",ParentalIllnessFile,SiblingIllnessCols[1]+1,SiblingIllnessCols[2]+1,"SiblingIllness.txt")
system(cmd)
SiblingIllness=data.frame(read.table("SiblingIllness.txt",header=TRUE,sep="\t"))
DEM$DEMSibling[rowSums(SiblingIllness==10, na.rm = TRUE)>0]=1 # 10 is code for dementia

# get subjects by self-reported diagnosis
cmd=sprintf("tail -n +2 %s | cut -f 1,%d-%d > selfDiagnoses.txt",exomesFile,selfReportFirst+1,selfReportLast+1)
system(cmd)
selfReport=data.frame(read.table("selfDiagnoses.txt",header=FALSE,sep="\t"))
DEM$DEMSelf=0
for (r in 1:length(selfReportedCodes)) {
  DEM$DEMSelf[rowSums(selfReport==selfReportedCodes[r], na.rm = TRUE)>0]=1
}

DEM$DEMScore=((DEM$DEMICD10+DEM$DEMICD9+DEM$DEMSelf)>0)*2+DEM$DEMFather+DEM$DEMMother+DEM$DEMSibling
DEM$DEMMax2=DEM$DEMScore
DEM$DEMMax2[DEM$DEMMax2>2]=2
for (c in 2:ncol(DEM)) {
  toWrite=DEM[,c(1,c)]
  write.table(toWrite,sprintf("UKBB.%s.20240901.txt",colnames(DEM)[c]),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
}

write.table(DEM,"UKBB.allFAMDEM.20240829.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)


