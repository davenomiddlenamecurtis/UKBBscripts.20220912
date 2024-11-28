#!/share/apps/R-3.6.1/bin/Rscript

# script to get data files to analyse association with FAMDEM

# note that the column number to provide is one higher than that given in http://www.davecurtis.net/UKBB/ukb41465.html
exomesFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb41465.470Kexomes.txt"
ParentalIllnessFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb43357.470Kexomes.txt"
diagnosisFile="/home/lgibbons/FAMDEM/DEMDiagnosis.txt" # ICD10 codes for dementia
ICD10SourcesFile="/home/rejudcu/UKBB/UKBBscripts.20220912/ICD10Sources.txt" # sources of ICD10 codes - hospital admissions, causes of death, etc
wd="/home/lgibbons/FAMDEM"
FatherIllnessCols=c(7002,7041)
MotherIllnessCols=c(7042,7085)
# 10 is Alzheimer/dementia
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

DEM=data.frame(matrix(ncol=6,nrow=nrow(ICD10)))
colnames(DEM)=c("IID","DEMICD10","DEMFather","DEMMother","DEMScore","DEMMax2")

diagnoses=data.frame(read.table(diagnosisFile,header=FALSE,sep="\t")) # ICD10 diagnoses
DEM$IID=ICD10$eid # I think they are all in the same order
DEM$DEMICD10=0
for (r in 1:nrow(diagnoses)) {
  DEM$DEMICD10[rowSums(ICD10==as.character(diagnoses[r,1]), na.rm = TRUE)>0]=1
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

# remove rows with missing parental data, or use an average. or just use all, I think

DEM$DEMScore=DEM$DEMICD10*2+DEM$DEMFather+DEM$DEMMother
DEM$DEMMax2=DEM$DEMScore
DEM$DEMMax2[DEM$DEMMax2>2]=2
for (c in 2:ncol(DEM)) {
  toWrite=DEM[,c(1,c)]
  write.table(toWrite,sprintf("UKBB.%s.20240528.txt",colnames(DEM)[c]),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
}

write.table(DEM,"UKBB.allFAMDEM.20240528.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)


