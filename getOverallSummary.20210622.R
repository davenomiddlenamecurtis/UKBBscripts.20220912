#!/share/apps/R-3.6.1/bin/Rscript

# Rscript getClinicalSummary.R UKBB.gene.txt UKBB.gene.summ.txt SummaryTable.txt
# Summary table has fields Field First Last DataScheme Description
library(data.table)

genes=c(
"SETD1A",
"CUL1",
"XPO7",
"TRIO",
"CACNA1G",
"SP4",
"GRIA3",
"GRIN2A",
"HERC1",
"RB1CC1"
)

# genes="TRIO"
date="20210615"

model=sprintf("LOF.%s",date)

wd="C:/Users/dave/OneDrive/sharedseq/UKBB/LOF"
wd="C:/Users/dave_000/OneDrive/sharedseq/UKBB/LOF"
# setwd(wd)

SummarySchemeFile="overallSummaryScheme.20210622.txt"
OutputFile="overallSummary.20210622.txt"
SummaryScheme=data.frame(read.table(SummarySchemeFile,header=TRUE,sep="\t",stringsAsFactors=FALSE))

nSCHEMAcase=24248 
nSCHEMAcontrol=97322
nUKBB=200632 

largestCoding=200000
codings= vector("list", largestCoding)
SCHEMATable=data.frame(read.table("SCHEMA.counts.txt",header=FALSE,sep="\t",stringsAsFactors=FALSE))
colnames(SCHEMATable)=c("Gene","CaseCount","ControlCount")

rowStart=vector()
rNames=vector()
rNames[1]="N"
rNames[2]="UKBBPrev"
rNames[3]="SCHEMACasePrev"
rNames[4]="SCHEMAControlPrev"
totRow=4
for (r in 1:nrow(SummaryScheme)) {
	totRow=totRow+1
	rowStart[r]=totRow
	rNames[totRow]=SummaryScheme$Name[r]
	if (SummaryScheme$Coding[r]!=0) {
		fn=sprintf("coding%d.tsv",SummaryScheme$Coding[r])
		coding=data.frame(read.table(fn,header=TRUE,sep="\t",quote="",stringsAsFactors=FALSE))
	    codings[[SummaryScheme$Coding[r]]]=coding
		for (rr in 1:nrow(coding)) {
			totRow=totRow+1
			rNames[totRow]=coding$meaning[rr]
		}
	}
}

SummaryTab=data.frame(matrix(ncol=length(genes)+1,nrow=length(rNames),0))
colnames(SummaryTab)=c("Variable",genes)
SummaryTab[,1]=rNames

c=2 # start in second column
for (gene in genes) {
goodUKBBdataFile=sprintf("%s.%s.phenos.good.txt",model,gene)
UKBBdata=data.frame(fread(goodUKBBdataFile))
N=nrow(UKBBdata)
SummaryTab[1,c]=N
SummaryTab[2,c]=N*100000/nUKBB
SummaryTab[3,c]=SCHEMATable[which(SCHEMATable$Gene==gene),2]*100000/nSCHEMAcase
SummaryTab[4,c]=SCHEMATable[which(SCHEMATable$Gene==gene),3]*100000/nSCHEMAcontrol
for (r in 1:nrow(SummaryScheme)){
  toMatch=sprintf("X%s",chartr("-",".",SummaryScheme$First[r])) # because R converts column names
  if (SummaryScheme$Coding[r]!=0) {
	codes=codings[[SummaryScheme$Coding[r]]]$coding
  } else {
    n=0
  }
  st=which(colnames(UKBBdata)==toMatch)
  for (s in 1:nrow(UKBBdata)) {
	val=UKBBdata[s,st]
	if (is.na(val)) next
	if (SummaryScheme$Coding[r]==0) {
		if (val<0) next # code for missing values, none of the ones I am using are genuinely negative
		# NB this would not work for e.g. Townsend index
		SummaryTab[rowStart[r],c]=SummaryTab[rowStart[r],c]+val
		n=n+1
	} else{
		SummaryTab[rowStart[r]+which(codes==val),c]=SummaryTab[rowStart[r]+which(codes==val),c]+1
	}
  }
  if (SummaryScheme$Coding[r]==0) {
    if (n>0) {
      SummaryTab[rowStart[r],c]=SummaryTab[rowStart[r],c]/n
	} else {
	  SummaryTab[rowStart[r],c]=0
	}
  }
}
for (r in 1:nrow(SummaryScheme)){
  if (SummaryScheme$Coding[r]!=0) {
	SummaryTab[rowStart[r],c]=""
  }
}
c=c+1
}

write.table(SummaryTab,OutputFile,sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
