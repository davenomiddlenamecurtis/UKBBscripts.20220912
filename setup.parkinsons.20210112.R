#!/share/apps/R-3.6.1/bin/Rscript

# script to get data files to analyse association with Parkinson's

# note that the column number to provide is one higher than that given in http://www.davecurtis.net/UKBB/ukb41465.html
exomesFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb41465.exomes.20210202.txt"
diagnosisFile="../parkinsons.2023/parkinsons.ICD10.codes.tsv" # ICD10 codes for various parkinsons diagnoses
# http://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=19
ICD10SourcesFile="ICD10Sources.txt" # sources of ICD10 codes - hospital admissions, causes of death
# http://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=6
# Data Coding 6
# 1220	diabetes	Yes	1245	1075
# 1221	gestational diabetes	Yes	1246	1245
# 1222	type 1 diabetes	Yes	1247	1245
# 1223	type 2 diabetes	Yes	1248	1245
# 1265	migraine	Yes	1291	1265
selfReportedCodes=c(1262) # self-report code(s) for Parkinsons
selfReportFirst=1730
selfReportLast=1763

femaleRxCols=c(1472:1475)
maleRxCols=c(1584:1586)

prescribedRxFirst=1866
prescribedRxLast=1913

# setwd(wd)


# first get subjects by self-reported diagnosis
if (!file.exists("../parkinsons.2023/selfDiagnoses.txt")) {
cmd=sprintf("tail -n +2 %s | cut -f 1,%d-%d > ../parkinsons.2023/selfDiagnoses.txt",exomesFile,selfReportFirst+1,selfReportLast+1)
system(cmd)
}
selfReport=data.frame(read.table("../parkinsons.2023/selfDiagnoses.txt",header=FALSE,sep="\t",stringsAsFactors=FALSE))
parkinsons=data.frame(matrix(ncol=4,nrow=nrow(selfReport)))
colnames(parkinsons)=c("IID","parkinsons","parkinsonsself","parkinsonsICD10")
parkinsons$IID=selfReport[,1]
parkinsons$parkinsonsself=0
for (r in 1:length(selfReportedCodes)) {
  parkinsons$parkinsonsself[rowSums(selfReport==selfReportedCodes[r], na.rm = TRUE)>0]=1
}

if (!file.exists("../parkinsons.2023/ICD10Codes.txt")) {
ICD10Sources=data.frame(read.table(ICD10SourcesFile,header=TRUE,sep="\t"))
cmd=sprintf("tail -n +2 %s | cut -f 1",exomesFile)
for (r in 1:nrow(ICD10Sources)) {
  cmd=sprintf("%s,%d-%d",cmd,ICD10Sources$First[r]+1,ICD10Sources$Last[r]+1)
}
cmd=sprintf("%s > ../parkinsons.2023/ICD10Codes.txt",cmd)
system(cmd)
}

ICD10=data.frame(read.table("../parkinsons.2023/ICD10Codes.txt",header=FALSE,sep="\t",stringsAsFactors=FALSE))
parkinsons$parkinsonsICD10=0

diagnoses=data.frame(read.table(diagnosisFile,header=TRUE,sep="\t",stringsAsFactors=FALSE)) # ICD10 diagnoses
for (r in 1:nrow(diagnoses)) {
  parkinsons$parkinsonsICD10[rowSums(ICD10==as.character(diagnoses[r,1]), na.rm = TRUE)>0]=1
}


parkinsons$parkinsons=0
parkinsons$parkinsons[rowSums(parkinsons[,3:4]==1, na.rm = TRUE)>0]=1
for (c in 2:ncol(parkinsons)) {
  toWrite=parkinsons[,c(1,c)]
  write.table(toWrite,sprintf("../parkinsons.2023/UKBB.%s.txt",colnames(parkinsons)[c]),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
}
write.table(parkinsons,"../parkinsons.2023/UKBB.parkinsonsall.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)


