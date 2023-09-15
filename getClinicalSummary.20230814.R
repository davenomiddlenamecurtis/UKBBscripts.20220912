#!/share/apps/R-3.6.1/bin/Rscript

# Rscript getClinicalSummary.R IDlist.txt SummaryTable.txt OutFileRoot
# Summary table has fields Field First Last DataScheme Description

# does no QC
# does not match ID to varfiant
# see earlier version for these functions

library(data.table)
exomesFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb41465.470Kexomes.txt"

date="20230906"

args = commandArgs(trailingOnly=TRUE)
if (length(args)< 3) {
	IDlistFile="HERC1.LOF.carriers.new.20230906.txt"
	SummarySchemeFile="/home/rejudcu/UKBB/LOF/SummaryScheme.20210618.txt"
	OutFileRoot="HERC1.LOF.new"
} else {
	IDlistFile=args[1]
	SummarySchemeFile=args[2]
	OutFileRoot=args[3]
}

model=sprintf("%s.%s",OutFileRoot,date)

wd="C:/Users/dave/OneDrive/sharedseq/UKBB/LOF"
# setwd(wd)
SummarySchemeFile="/home/rejudcu/UKBB/LOF/SummaryScheme.20210618.txt"
SummaryScheme=data.frame(read.table(SummarySchemeFile,header=TRUE,sep="\t",stringsAsFactors=FALSE))

largestCoding=200000
codings= vector("list", largestCoding)

UKBBdataFile=sprintf("%s.phenos.txt",model)
OutputFile=sprintf("SummaryFor.%s.txt",model)

if (!file.exists(IDlistFile)) {
	print(sprintf("Could not find file %s",IDlistFile))
	quit()
}

if (!file.exists(UKBBdataFile)) {
CommLine="echo eid$'\t'extra > extraField.txt"
print(CommLine)
system(CommLine)
CommLine=sprintf("head -n 1 %s | join -t  $'\t' - extraField.txt > %s",exomesFile,UKBBdataFile)
# this is because the exomes file has got a tab appended to each line after the header
# and otherwise fread will not work
print(CommLine)
system(CommLine)
CommLine=sprintf("tail -n +2 %s | join -t  $'\t' %s - >> %s",exomesFile,IDlistFile,UKBBdataFile) 
print(CommLine)
system(CommLine)
}

UKBBdata=data.frame(fread(UKBBdataFile,stringsAsFactors=FALSE))

SummaryTab=data.frame(matrix(ncol=nrow(UKBBdata),nrow=nrow(SummaryScheme)))
rownames(SummaryTab)=SummaryScheme$Name
SummaryTab[2,]="BlankAnnot" 
SummaryTab[3,]="BlankAnnot"
for (r in 1:nrow(SummaryScheme)){
  if (r==2|| r==3) { next } # for variant
  st=NA
  en=NA

  if (r==1) {
    toMatch=SummaryScheme$First[r]
  } else {
    toMatch=sprintf("X%s",chartr("-",".",SummaryScheme$First[r])) # because R converts column names
  }
  for (s in 1:ncol(UKBBdata)) {
    if (colnames(UKBBdata)[s]==toMatch) {
      st=s
      break
    }
  }
  if (r==1) {
    toMatch=SummaryScheme$Last[r]
  } else {
    toMatch=sprintf("X%s",chartr("-",".",SummaryScheme$Last[r]))
  }
  for (s in st:ncol(UKBBdata)) {
    if (colnames(UKBBdata)[s]==toMatch) {
      en=s
      break
    }
  }
  for (s in 1:nrow(UKBBdata)) {
      if (s==1 & SummaryScheme$Coding[r] !=0) {
	    if (is.null(codings[[SummaryScheme$Coding[r]]])) {
	      fn=sprintf("coding%d.tsv",SummaryScheme$Coding[r])
	      codings[[SummaryScheme$Coding[r]]]=data.frame(read.table(fn,header=TRUE,sep="\t",quote="",stringsAsFactors=FALSE))
		}
		codes=codings[[SummaryScheme$Coding[r]]]
	  }
      first=TRUE
	  output="Missing"
	  for (c in st:en){
	    val=UKBBdata[s,c]
	    if (is.na(val)) {
	      next
	    }
	    if (val=="") {
	      next
	    }
	    if (SummaryScheme$Coding[r]!=0) {
	      code=codes[codes$coding==val,]
	      val=code[1,2]
		}
	    if (first==TRUE) {
		  output=val
		  first=FALSE
		} else {
		  output=sprintf("%s; %s",output,val)
		}
	  }
	  SummaryTab[r,s]=output
  }
}

write.table(SummaryTab,OutputFile,sep="\t",col.names=FALSE,row.names=TRUE,quote=FALSE)

