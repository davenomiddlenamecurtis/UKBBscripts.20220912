#!/share/apps/R-3.6.1/bin/Rscript

# script to get data files to analyse association with FAMLOAD

# note that the column number to provide is one higher than that given in http://www.davecurtis.net/UKBB/ukb41465.html
exomesFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb41465.exomes.20201103.txt"
ParentalIllnessFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb43357.exomes.20201103.txt"
diagnosisFile="/home/rejudcu/UKBB/UKBBscripts/LOADDiagnosis.txt" # ICD10 codes for Alzheimer/dementia
ICD10SourcesFile="/home/rejudcu/UKBB/UKBBscripts/ICD10Sources.txt" # sources of ICD10 codes - hospital admissions, causes of death, etc
wd="/home/rejudcu/UKBB/LOAD"
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

LOAD=data.frame(matrix(ncol=6,nrow=nrow(ICD10)))
colnames(LOAD)=c("IID","LOADICD10","LOADFather","LOADMother","LOADScore","LOADMax2")

diagnoses=data.frame(read.table(diagnosisFile,header=FALSE,sep="\t")) # ICD10 diagnoses
LOAD$IID=ICD10$eid # I think they are all in the same order
LOAD$LOADICD10=0
for (r in 1:nrow(diagnoses)) {
  LOAD$LOADICD10[rowSums(ICD10==as.character(diagnoses[r,1]), na.rm = TRUE)>0]=1
}

LOAD$LOADFather=0
cmd=sprintf("cut %s -f 1,%d-%d > %s",ParentalIllnessFile,FatherIllnessCols[1]+1,FatherIllnessCols[2]+1,"FatherIllness.txt")
system(cmd)
FatherIllness=data.frame(read.table("FatherIllness.txt",header=TRUE,sep="\t"))
LOAD$LOADFather[rowSums(FatherIllness==10, na.rm = TRUE)>0]=1 # 10 is code for dementia

LOAD$LOADMother=0
cmd=sprintf("cut %s -f 1,%d-%d > %s",ParentalIllnessFile,MotherIllnessCols[1]+1,MotherIllnessCols[2]+1,"MotherIllness.txt")
system(cmd)
MotherIllness=data.frame(read.table("MotherIllness.txt",header=TRUE,sep="\t"))
LOAD$LOADMother[rowSums(MotherIllness==10, na.rm = TRUE)>0]=1 # 10 is code for dementia

# remove rows with missing parental data, or use an average. or just use all, I think

LOAD$LOADScore=LOAD$LOADICD10*2+LOAD$LOADFather+LOAD$LOADMother
LOAD$LOADMax2=LOAD$LOADScore
LOAD$LOADMax2[LOAD$LOADMax2>2]=2
for (c in 2:ncol(LOAD)) {
  toWrite=LOAD[,c(1,c)]
  write.table(toWrite,sprintf("UKBB.%s.20201208.txt",colnames(LOAD)[c]),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
}

write.table(LOAD,"UKBB.allFAMLOAD.20201208.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)


