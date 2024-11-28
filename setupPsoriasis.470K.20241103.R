#!/share/apps/R-3.6.1/bin/Rscript

# script to get data files to analyse association with Psoriasis

# note that the column number to provide is one higher than that given in http://www.davecurtis.net/UKBB/ukb41465.html
exomesFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb41465.470Kexomes.txt"
diagnosisFile="psoriasisDiagnosis.txt" # ICD10 codes for various Psoriasis diagnoses
PsoriasisCodes=c(1453) # code for self-reported psoriasis

# http://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=19
ICD10SourcesFile="ICD10Sources.txt" # sources of ICD10 codes - hospital admissions, causes of death
wd="/home/rejudcu/UKBB/UKBBscripts.20220912"

SelfReportFirst=1730
SelfReportLast=1763


setwd(wd)

# uses this: http://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=4&nl=1
# cmd=sprintf( "grep -Ff %s %s > %s",drugNamesFile,drugCodesFile,PsoriasisDrugsFile)
# system(cmd)

# first get subjects by Self-reported diagnosis
cmd=sprintf("tail -n +2 %s | cut -f 1,%d-%d > ../psoriasis.470K.20241103/SelfDiagnoses.txt",exomesFile,SelfReportFirst+1,SelfReportLast+1)
system(cmd)
SelfReport=data.frame(read.table("../psoriasis.470K.20241103/SelfDiagnoses.txt",header=FALSE,sep="\t"))
Psoriasis=data.frame(matrix(ncol=4,nrow=nrow(SelfReport)))
colnames(Psoriasis)=c("IID","Psoriasis","PsoriasisSelf","PsoriasisICD10")
Psoriasis$IID=SelfReport[,1]
Psoriasis$PsoriasisSelf=0
for (r in 1:length(PsoriasisCodes)) {
  Psoriasis$PsoriasisSelf[rowSums(SelfReport==PsoriasisCodes[r], na.rm = TRUE)>0]=1
}

ICD10Sources=data.frame(read.table(ICD10SourcesFile,header=TRUE,sep="\t"))
cmd=sprintf("tail -n +2 %s | cut -f 1",exomesFile)
for (r in 1:nrow(ICD10Sources)) {
  cmd=sprintf("%s,%d-%d",cmd,ICD10Sources$First[r]+1,ICD10Sources$Last[r]+1)
}
cmd=sprintf("%s >../psoriasis.470K.20241103/ICD10Codes.txt",cmd)
system(cmd)
ICD10=data.frame(read.table("../psoriasis.470K.20241103/ICD10Codes.txt",header=FALSE,sep="\t"))

Psoriasis$PsoriasisICD10=0

diagnoses=data.frame(read.table(diagnosisFile,header=TRUE,sep="\t")) # ICD10 diagnoses
for (r in 1:nrow(diagnoses)) {
  Psoriasis$PsoriasisICD10[rowSums(ICD10==as.character(diagnoses[r,1]), na.rm = TRUE)>0]=1
}
# for HT, diagnosis codes were in second column but here they are in first

Psoriasis$Psoriasis=0
Psoriasis$Psoriasis[rowSums(Psoriasis[,3:4]==1, na.rm = TRUE)>0]=1
for (c in 2:ncol(Psoriasis)) {
  toWrite=Psoriasis[,c(1,c)]
  write.table(toWrite,sprintf("../psoriasis.470K.20241103/UKBB.%s.txt",colnames(Psoriasis)[c]),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
}
write.table(Psoriasis,"../psoriasis.470K.20241103/UKBB.Psoriasisall.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)


