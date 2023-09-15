#!/share/apps/R-3.6.1/bin/Rscript

# script to get data files to analyse association with asthma, atopic dermatitis, hay fever / allergic rhinitis

# note that the column number to provide is one higher than that given in http://www.davecurtis.net/UKBB/ukb41465.html
exomesFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb41465.exomes.20210202.txt"
asthmaDiagnosisFile="asthmaDiagnosis.txt" # ICD10 codes for various asthma diagnoses
dermatitisDiagnosisFile="atopicDermatitisDiagnosis.txt" # ICD10 codes for various atopic dermatitis diagnoses
# http://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=19
AllergicRhinitisCodes=c(1387) # code for self-reported hay fever / allergic rhinitis

ICD10SourcesFile="ICD10Sources.txt" # sources of ICD10 codes - hospital admissions, causes of death
wd="/home/rejudcu/UKBB/UKBBscripts.20220912"

selfReportFirst=1730
selfReportLast=1763

setwd(wd)

# first get subjects by self-reported diagnosis
cmd=sprintf("tail -n +2 %s | cut -f 1,%d-%d > ../asthma.20230313/selfDiagnoses.txt",exomesFile,selfReportFirst+1,selfReportLast+1)
system(cmd)
selfReport=data.frame(read.table("../asthma.20230313/selfDiagnoses.txt",header=FALSE,sep="\t"))
AllergicRhinitis=data.frame(matrix(ncol=2,nrow=nrow(selfReport)))
colnames(AllergicRhinitis)=c("IID","AllergicRhinitis")
AllergicRhinitis$IID=selfReport[,1]
AllergicRhinitis$AllergicRhinitis=0
for (r in 1:length(AllergicRhinitisCodes)) {
# in fact, there is only one code
  AllergicRhinitis$AllergicRhinitis[rowSums(selfReport==AllergicRhinitisCodes[r], na.rm = TRUE)>0]=1
}

ICD10Sources=data.frame(read.table(ICD10SourcesFile,header=TRUE,sep="\t"))
cmd=sprintf("tail -n +2 %s | cut -f 1",exomesFile)
for (r in 1:nrow(ICD10Sources)) {
  cmd=sprintf("%s,%d-%d",cmd,ICD10Sources$First[r]+1,ICD10Sources$Last[r]+1)
}
cmd=sprintf("%s > ../asthma.20230313/ICD10Codes.txt",cmd)
system(cmd)

ICD10=data.frame(read.table("../asthma.20230313/ICD10Codes.txt",header=FALSE,sep="\t"))

Asthma=data.frame(matrix(ncol=2,nrow=nrow(ICD10)))
colnames(Asthma)=c("IID","Asthma")
diagnoses=data.frame(read.table(asthmaDiagnosisFile,header=TRUE,sep="\t",quote="")) # ICD10 diagnoses
Asthma$IID=ICD10[,1]
Asthma$Asthma=0
for (r in 1:nrow(diagnoses)) {
  Asthma$Asthma[rowSums(ICD10==as.character(diagnoses[r,1]), na.rm = TRUE)>0]=1
}

Dermatitis=data.frame(matrix(ncol=2,nrow=nrow(ICD10)))
colnames(Asthma)=c("IID","Dermatitis")
diagnoses=data.frame(read.table(dermatitisDiagnosisFile,header=TRUE,sep="\t",quote="")) # ICD10 diagnoses
Dermatitis$IID=ICD10[,1]
Dermatitis$Dermatitis=0
for (r in 1:nrow(diagnoses)) {
  Dermatitis$Dermatitis[rowSums(ICD10==as.character(diagnoses[r,1]), na.rm = TRUE)>0]=1
}

write.table(AllergicRhinitis,"../asthma.20230313/UKBB.AllergicRhinitis.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
write.table(Asthma,"../asthma.20230313/UKBB.Asthma.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
write.table(Dermatitis,"../asthma.20230313/UKBB.Dermatitis.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

