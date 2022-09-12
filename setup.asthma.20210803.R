#!/share/apps/R-3.6.1/bin/Rscript

# script to get data files to analyse association with asthma

# note that the column number to provide is one higher than that given in http://www.davecurtis.net/UKBB/ukb41465.html
exomesFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb41465.exomes.20210202.txt"

# we will just use two fields as sources
# only use the first instance, which is initial assessment
# this is expected to miss some late onset cases
ageAtopy=908 # field 3761-0.0
ageAsthma=916 # field 3786-0.0

diagnosisFile="asthmaDiagnosis.txt" # ICD10 codes for various asthma diagnoses
# http://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=19
ICD10SourcesFile="ICD10Sources.txt" # sources of ICD10 codes - hospital admissions, causes of death
wd="/home/rejudcu/UKBB/UKBBscripts"

setwd(wd)

# first get ages of diagnosis
cmd=sprintf("tail -n +2 %s | cut -f 1,%d,%d > ../asthma/ageDiagnoses.txt",exomesFile,ageAtopy+1,ageAsthma+1)
system(cmd)
ageDiag=data.frame(read.table("../asthma/ageDiagnoses.txt",header=FALSE,sep="\t"))
ageDiag[is.na(ageDiag)]=-2
asthma=data.frame(matrix(ncol=6,nrow=nrow(ageDiag),0))
colnames(asthma)=c("IID","earlyAsthma","lateAsthma","atopicAsthma","nonatopicAsthma")

asthma$IID=ageDiag[,1]
asthma$earlyAsthma=1
asthma$earlyAsthma[ageDiag[,2]<0]=0
asthma$earlyAsthma[ageDiag[,2]>=16]=0

for (c in 2:ncol(asthma)) {
  toWrite=asthma[,c(1,c)]
  write.table(toWrite,sprintf("../asthma/UKBB.%s.txt",colnames(asthma)[c]),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
}
write.table(asthma,"../asthma/UKBB.asthmaAll.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)


