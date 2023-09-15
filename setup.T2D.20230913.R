#!/share/apps/R-3.6.1/bin/Rscript

# script to get data files to analyse association with T2D

# note that the column number to provide is one higher than that given in http://www.davecurtis.net/UKBB/ukb41465.html
exomesFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb41465.470Kexomes.txt"
diagnosisFile="T2DDiagnosis.txt" # ICD10 codes for various T2D diagnoses
# http://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=19
ICD10SourcesFile="ICD10Sources.txt" # sources of ICD10 codes - hospital admissions, causes of death
drugNamesFile="T2DDrugNames.txt" # list of hypoglycaemics
drugCodesFile="UKBB.coding4.tsv"
wd="/home/rejudcu/UKBB/UKBBscripts.20220912"
# http://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=6
# Data Coding 6
# 1220	diabetes	Yes	1245	1075
# 1221	gestational diabetes	Yes	1246	1245
# 1222	type 1 diabetes	Yes	1247	1245
# 1223	type 2 diabetes	Yes	1248	1245
selfReportedCodes=c(1220,1223)
selfReportFirst=1730
selfReportLast=1763

femaleRxCols=c(1472:1475)
maleRxCols=c(1584:1586)

prescribedRxFirst=1866
prescribedRxLast=1913

setwd(wd)

hypoglycaemicsFile="hypoglycaemics.txt"
# uses this: http://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=4&nl=1
cmd=sprintf( "grep -Ff %s %s > %s",drugNamesFile,drugCodesFile,hypoglycaemicsFile)
system(cmd)

# first get subjects by self-reported diagnosis
cmd=sprintf("tail -n +2 %s | cut -f 1,%d-%d > ../T2D.20210104/selfDiagnoses.txt",exomesFile,selfReportFirst+1,selfReportLast+1)
system(cmd)
selfReport=data.frame(read.table("../T2D.20210104/selfDiagnoses.txt",header=FALSE,sep="\t"))
T2D=data.frame(matrix(ncol=5,nrow=nrow(selfReport)))
colnames(T2D)=c("IID","T2D","T2Dself","T2DICD10","RxIV")
T2D$IID=selfReport[,1]
T2D$T2Dself=0
for (r in 1:length(selfReportedCodes)) {
  T2D$T2Dself[rowSums(selfReport==selfReportedCodes[r], na.rm = TRUE)>0]=1
}

ICD10Sources=data.frame(read.table(ICD10SourcesFile,header=TRUE,sep="\t"))
cmd=sprintf("tail -n +2 %s | cut -f 1",exomesFile)
for (r in 1:nrow(ICD10Sources)) {
  cmd=sprintf("%s,%d-%d",cmd,ICD10Sources$First[r]+1,ICD10Sources$Last[r]+1)
}
cmd=sprintf("%s > ../T2D.20210104/ICD10Codes.txt",cmd)
system(cmd)

ICD10=data.frame(read.table("../T2D.20210104/ICD10Codes.txt",header=FALSE,sep="\t"))
T2D$T2DICD10=0

diagnoses=data.frame(read.table(diagnosisFile,header=FALSE,sep="\t")) # ICD10 diagnoses
for (r in 1:nrow(diagnoses)) {
  T2D$T2DICD10[rowSums(ICD10==as.character(diagnoses[r,1]), na.rm = TRUE)>0]=1
}
# for HT, diagnosis codes were in second column but here they are in first

cmd=sprintf("tail -n +2 %s | cut -f 1,%d-%d > ../T2D.20210104/prescribedRx.txt",exomesFile,prescribedRxFirst+1,prescribedRxLast+1)
system(cmd)
prescribedRx=data.frame(read.table("../T2D.20210104/prescribedRx.txt",header=FALSE,sep="\t"))
hypoglycaemics=data.frame(read.table(hypoglycaemicsFile,header=FALSE,sep="\t")) # list of hypoglycaemic codes
T2D$RxIV=0
for (r in 1:nrow(hypoglycaemics)) {
  T2D$RxIV[rowSums(prescribedRx==hypoglycaemics[r,1], na.rm = TRUE)>0]=1
}

T2D$T2D=0
T2D$T2D[rowSums(T2D[,3:5]==1, na.rm = TRUE)>0]=1
for (c in 2:ncol(T2D)) {
  toWrite=T2D[,c(1,c)]
  write.table(toWrite,sprintf("../T2D.2023/UKBB.%s.txt",colnames(T2D)[c]),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
}
write.table(T2D,"../T2D.2023/UKBB.T2Dall.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

phenos=as.matrix(T2D[,2:5])
print(cor(phenos))

