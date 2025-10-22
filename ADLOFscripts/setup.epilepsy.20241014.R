#!/share/apps/R-3.6.1/bin/Rscript

# script to get data files to analyse association with epilepsy

# note that the column number to provide is one higher than that given in http://www.davecurtis.net/UKBB/ukb41465.html
# exomesFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb43357.txt"
# the raw UK Biobank file has a leading tab, FFS
# so add 2 to all the column numbers below
# change to using text file with only data for exome-sequenced subjects, so only need to offset by 1
exomesFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb43357.470Kexomes.txt"

# I changed to use ukb43357 instead of ukb41465 so I had to change to the new column numbers

diagnosisFile="/home/mayuan/scripts/epilepsyDiagnosis.txt" # ICD10 codes for various epilepsy diagnoses
# http://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=19
ICD10SourcesFile="/home/mayuan/scripts/ICD10Sources.20230811.txt" # sources of ICD10 codes - hospital admissions, causes of death
wd="/home/mayuan/Epilepsy.20241014"

# http://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=6
# Data Coding 6
# 1264	epilepsy
selfReportedCodes=c(1264)
selfReportFirst=5634	
selfReportLast=5769

femaleRxCols=c(5325:5340)
maleRxCols=c(5437:5448)

prescribedRxFirst=5770
prescribedRxLast=5961

setwd(wd)


# first get subjects by self-reported diagnosis
cmd=sprintf("tail -n +2 %s | cut -f 1,%d-%d > /home/mayuan/Epilepsy.20241014/selfDiagnoses.txt",exomesFile,selfReportFirst+1,selfReportLast+1)
system(cmd)
selfReport=data.frame(read.table("/home/mayuan/Epilepsy.20241014/selfDiagnoses.txt",header=FALSE,sep="\t"))
epilepsy=data.frame(matrix(ncol=4,nrow=nrow(selfReport)))
colnames(epilepsy)=c("IID","Epilepsy","EpilepsySelf","EpilepsyICD10")
epilepsy$IID=selfReport[,1]
epilepsy$EpilepsySelf=0
for (r in 1:length(selfReportedCodes)) {
  epilepsy$EpilepsySelf[rowSums(selfReport==selfReportedCodes[r], na.rm = TRUE)>0]=1
}

ICD10Sources=data.frame(read.table(ICD10SourcesFile,header=TRUE,sep="\t"))
cmd=sprintf("tail -n +2 %s | cut -f 1",exomesFile)
for (r in 1:nrow(ICD10Sources)) {
  cmd=sprintf("%s,%d-%d",cmd,ICD10Sources$First[r]+1,ICD10Sources$Last[r]+1)
}
cmd=sprintf("%s > /home/mayuan/Epilepsy.20241014/ICD10Codes.txt",cmd)
system(cmd)

ICD10=data.frame(read.table("/home/mayuan/Epilepsy.20241014/ICD10Codes.txt",header=FALSE,sep="\t"))
epilepsy$EpilepsyICD10=0

diagnoses=data.frame(read.table(diagnosisFile,header=FALSE,sep="\t")) # ICD10 diagnoses
for (r in 1:nrow(diagnoses)) {
  epilepsy$EpilepsyICD10[rowSums(ICD10==as.character(diagnoses[r,1]), na.rm = TRUE)>0]=1
}

epilepsy$Epilepsy=0
epilepsy$Epilepsy[rowSums(epilepsy[,3:4]==1, na.rm = TRUE)>0]=1
for (c in 2:ncol(epilepsy)) {
  toWrite=epilepsy[,c(1,c)]
  write.table(toWrite,sprintf("/home/mayuan/Epilepsy.20241014/UKBB.%s.txt",colnames(epilepsy)[c]),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
}
write.table(epilepsy,"/home/mayuan/Epilepsy.20241014/UKBB.EpilepsyAll.20241014.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

phenos=as.matrix(epilepsy[,2:4])
print(cor(phenos))

