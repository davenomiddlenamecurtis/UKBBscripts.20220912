#!/share/apps/R-3.6.1/bin/Rscript

# script to get data files to analyse association with cancer

# note that the column number to provide is one higher than that given in http://www.davecurtis.net/UKBB/ukb41465.html
exomesFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb41465.exomes.20210202.txt"
CancerCodes=c(1453) # code for self-reported psoriasis

# http://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=19
ICD10SourcesFile="ICD10Sources.txt" # sources of ICD10 codes - hospital admissions, causes of death
wd="/home/rejudcu/UKBB/UKBBscripts.20220912"

SelfReportCancerFirst=1706
SelfReportCancerLast=1729
SelfReportNumCancersFirst=371
SelfReportNumCancersLast=374
YoBCol=6
DoACol=28
DoDFirst=4499 # from death register
DoDLast=4500
AgeAtDeathFirst=4565
AgeAtDeathLast=4566


setwd(wd)

# uses this: http://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=4&nl=1
# cmd=sprintf( "grep -Ff %s %s > %s",drugNamesFile,drugCodesFile,CancerDrugsFile)
# system(cmd)

# first get subjects by Self-reported diagnosis
cmd=sprintf("tail -n +2 %s | cut -f 1,%d-%d > ../cancer.20260928/SelfDiagnosesCancer.txt",exomesFile,SelfReportCancerFirst+1,SelfReportCancerLast+1)
system(cmd)
SelfReport=data.frame(read.table("../cancer.20260928/SelfDiagnosesCancer.txt",header=FALSE,sep="\t"))
Cancer=data.frame(matrix(ncol=5,nrow=nrow(SelfReport)))
colnames(Cancer)=c("IID","Cancer","CancerSelf","CancerNumber","CancerICD10")
Cancer$IID=SelfReport[,1]
Cancer$CancerSelf=0
Cancer$CancerSelf[rowSums(SelfReport[,-1], na.rm = TRUE)!=0]=1 # actually anything not NA}

cmd=sprintf("tail -n +2 %s | cut -f 1,%d-%d > ../cancer.20260928/SelfDiagnosesNumCancers.txt",exomesFile,SelfReportNumCancersFirst+1,SelfReportNumCancersLast+1)
system(cmd)
SelfReportNum=data.frame(read.table("../cancer.20260928/SelfDiagnosesNumCancers.txt",header=FALSE,sep="\t"))
Cancer$CancerNumber=0
Cancer$CancerNumber[rowSums(SelfReportNum[,-1], na.rm = TRUE)>0]=1
# corr(Cancer$CancerNumber,Cancer$CancerSelf) == 1


# ICD10 codes: C* but not Ch*, D0*
ICD10Sources=data.frame(read.table(ICD10SourcesFile,header=TRUE,sep="\t"))
cmd=sprintf("tail -n +2 %s | cut -f 1",exomesFile)
for (r in 1:nrow(ICD10Sources)) {
  cmd=sprintf("%s,%d-%d",cmd,ICD10Sources$First[r]+1,ICD10Sources$Last[r]+1)
}
cmd=sprintf("%s >../cancer.20260928/ICD10Codes.txt",cmd)
system(cmd)
ICD10=data.frame(read.table("../cancer.20260928/ICD10Codes.txt",header=FALSE,sep="\t"))

ICD10First=ICD10
ICD10First[,2:ncol(ICD10First)]=as.data.frame(
  lapply(
    X=ICD10First[,2:ncol(ICD10First)],
	FUN=function(x) substr(x,1,1)
  )
)

ICD10Second=ICD10
ICD10Second[,2:ncol(ICD10Second)]=as.data.frame(
  lapply(
    X=ICD10Second[,2:ncol(ICD10Second)],
	FUN=function(x) substr(x,2,2)
  )
)


ICD10Cancer=(ICD10First=="C" & ICD10Second!="h") | (ICD10First=="D" & ICD10Second=="0")

Cancer$CancerICD10=0
Cancer$CancerICD10[rowSums(ICD10Cancer, na.rm = TRUE)>0]=1

Cancer$Cancer=0
Cancer$Cancer[rowSums(Cancer[,3:5]==1, na.rm = TRUE)>0]=1
for (c in 2:ncol(Cancer)) {
  toWrite=Cancer[,c(1,c)]
  write.table(toWrite,sprintf("../cancer.20260928/UKBB.%s.txt",colnames(Cancer)[c]),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
}
write.table(Cancer,"../cancer.20260928/UKBB.Cancerall.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

# get useful dates
cmd=sprintf("tail -n +2 %s | cut -f 1,%d,%d,%d-%d,%d-%d > ../cancer.20260928/Dates.txt",exomesFile,YoBCol+1,DoACol+1,DoDFirst+1,DoDLast+1,AgeAtDeathFirst+1,AgeAtDeathLast+1)
system(cmd)
DatesRaw=data.frame(read.table("../cancer.20260928/Dates.txt",header=FALSE,sep="\t",stringsAsFactors=FALSE))
# there are no rows where V4 is blank and V5 is not (same for V6 and V7)
DatesRaw[,3:4]=as.data.frame(
  lapply(
    X=DatesRaw[,3:4],
	FUN=function(x) substr(x,1,4)
  ),stringsAsFactors=FALSE
)

Dates=data.frame(matrix(ncol=8,nrow=nrow(DatesRaw)))

colnames(Dates)=c("IID","YoB","YoA","AgeAtAssessment","YoD","AgeAtDeath","TimeToEvent","Status")
Dates[,1:3]=DatesRaw[,1:3]
Dates$AgeAtAssessment=as.numeric(Dates$YoA)-as.numeric(Dates$YoB)
Dates$YoD=DatesRaw$V4
Dates$AgeAteDeath=DatesRaw$V6
MaxDeathDate=max(as.numeric(DatesRaw$V4[DatesRaw$V4!=""]))
MaxDeathDate
Dates$Status=0
Dates$Status[DatesRaw$V4!=""]=1
DatesRaw$TimeToEvent=MaxDeathDate-as.numeric(DatesRaw$V3)
Dates$TimeToEvent[DatesRaw$V4==""]=DatesRaw$TimeToEvent[DatesRaw$V4==""]
DatesRaw$TimeToEvent=as.numeric(DatesRaw$V4)-as.numeric(DatesRaw$V3)
Dates$TimeToEvent[DatesRaw$V4!=""]=DatesRaw$TimeToEvent[DatesRaw$V4!=""]
Dates$TimeToEvent[Dates$TimeToEvent<0]=0

write.table(Dates,"../cancer.20260928/UKBB.Dates.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)




