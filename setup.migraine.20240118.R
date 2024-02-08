#!/share/apps/R-3.6.1/bin/Rscript

# script to get data files to analyse association with migraine, using 470K exomes

# note that the column number to provide is one higher than that given in http://www.davecurtis.net/UKBB/ukb41465.html
exomesFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb41465.470Kexomes.txt"
diagnosisFile="migraineICDcodes.tsv" # ICD10 codes for various migraine diagnoses
OutputDir="/home/rejudcu/UKBB/migraine.20240118"

# http://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=19
ICD10SourcesFile="ICD10Sources.txt" # sources of ICD10 codes - hospital admissions, causes of death
drugNamesFile="migraineDrugs.tsv" # list of migraineDrugs
drugCodesFile="UKBB.coding4.tsv"
wd="/home/kmarkel/UKBB/UKBBscripts"
# http://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=6
# Data Coding 6
# 1220	diabetes	Yes	1245	1075
# 1221	gestational diabetes	Yes	1246	1245
# 1222	type 1 diabetes	Yes	1247	1245
# 1223	type 2 diabetes	Yes	1248	1245
# 1265	migraine	Yes	1291	1265
selfReportedCodes=c(1265)
selfReportFirst=1730
selfReportLast=1763

femaleRxCols=c(1472:1475)
maleRxCols=c(1584:1586)

prescribedRxFirst=1866
prescribedRxLast=1913

setwd(wd)

migraineDrugsFile=drugNamesFile

# uses this: http://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=4&nl=1
# cmd=sprintf( "grep -Ff %s %s > %s",drugNamesFile,drugCodesFile,migraineDrugsFile)
# system(cmd)

# first get subjects by self-reported diagnosis
cmd=sprintf("tail -n +2 %s | cut -f 1,%d-%d > %s/selfDiagnoses.txt",exomesFile,selfReportFirst+1,selfReportLast+1,OutputDir)
system(cmd)
selfReport=data.frame(read.table(sprintf("%s/selfDiagnoses.txt",OutputDir),header=FALSE,sep="\t"))
migraine=data.frame(matrix(ncol=5,nrow=nrow(selfReport)))
colnames(migraine)=c("IID","migraine","migraineself","migraineICD10","RxIV")
migraine$IID=selfReport[,1]
migraine$migraineself=0
for (r in 1:length(selfReportedCodes)) {
  migraine$migraineself[rowSums(selfReport==selfReportedCodes[r], na.rm = TRUE)>0]=1
}

ICD10Sources=data.frame(read.table(ICD10SourcesFile,header=TRUE,sep="\t"))
cmd=sprintf("tail -n +2 %s | cut -f 1",exomesFile)
for (r in 1:nrow(ICD10Sources)) {
  cmd=sprintf("%s,%d-%d",cmd,ICD10Sources$First[r]+1,ICD10Sources$Last[r]+1)
}
cmd=sprintf("%s > %s/ICD10Codes.txt",cmd,OutputDir)
system(cmd)

ICD10=data.frame(read.table(sprintf("%s/ICD10Codes.txt",OutputDir),header=FALSE,sep="\t"))
migraine$migraineICD10=0

diagnoses=data.frame(read.table(diagnosisFile,header=TRUE,sep="\t")) # ICD10 diagnoses
for (r in 1:nrow(diagnoses)) {
  migraine$migraineICD10[rowSums(ICD10==as.character(diagnoses[r,1]), na.rm = TRUE)>0]=1
}
# for HT, diagnosis codes were in second column but here they are in first

cmd=sprintf("tail -n +2 %s | cut -f 1,%d-%d > %s/prescribedRx.txt",exomesFile,prescribedRxFirst+1,prescribedRxLast+1,OutputDir)
system(cmd)
prescribedRx=data.frame(read.table(sprintf("%s/prescribedRx.txt",OutputDir),header=FALSE,sep="\t"))
migraineDrugs=data.frame(read.table(drugNamesFile,header=TRUE,sep="\t")) # list of migraine drug codes
migraine$RxIV=0
for (r in 1:nrow(migraineDrugs)) {
  migraine$RxIV[rowSums(prescribedRx==migraineDrugs[r,1], na.rm = TRUE)>0]=1
}

migraine$migraine=0
migraine$migraine[rowSums(migraine[,3:5]==1, na.rm = TRUE)>0]=1
for (c in 2:ncol(migraine)) {
  toWrite=migraine[,c(1,c)]
  write.table(toWrite,sprintf("%s/UKBB.%s.txt",OutputDir,colnames(migraine)[c]),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
}
write.table(migraine,sprintf("%s/UKBB.migraineall.txt",OutputDir),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)


