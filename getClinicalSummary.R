#!/share/apps/R-3.6.1/bin/Rscript

# Rscript getClinicalSummary.R UKBB.gene.txt UKBB.gene.summ.txt SummaryTable.txt
# Summary table has fields Field First Last DataScheme Description


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

args = commandArgs(trailingOnly=TRUE)
if (length(args)>0) {
  genes=args
}


wd="C:/Users/dave/OneDrive/sharedseq/UKBB/LOF"
# setwd(wd)
SummarySchemeFile="SummaryScheme.txt"
SummaryScheme=data.frame(read.table(SummarySchemeFile,header=TRUE,sep="\t"))

largestCoding=200000
codings= vector("list", largestCoding)
# genes="LDLR"
# genes="HIRA"

for (gene in genes) {

UKBBdataFile=sprintf("LOFs.%s.phenos.txt",gene)
annotsFile=sprintf("IDsAndAnnots.sorted.LOFs.%s.txt",gene)
UKBBdataFile=sprintf("LOF.20201227.%s.phenos.txt",gene)
annotsFile=sprintf("IDsAndAnnots.sorted.LOF.20201227.%s.txt",gene)
OutputFile=sprintf("SummaryFor.%s.txt",gene)

if (!file.exists(UKBBdataFile)) {
	next
}

UKBBdata=data.frame(read.table(UKBBdataFile,header=TRUE,sep="\t",fill=TRUE,quote="",stringsAsFactors=FALSE))
annots=data.frame(read.table(annotsFile,header=FALSE,sep=" ",stringsAsFactors=FALSE))
Summary=data.frame(matrix(ncol=nrow(UKBBdata),nrow=nrow(SummaryScheme)))
rownames(Summary)=SummaryScheme$Name
Summary[2,]=annots[,2] # just assume in same order
for (r in 1:nrow(SummaryScheme)){
  if (r==2) { next } # for variant
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
if (r==55) { print (codes) }
      first=TRUE
	  output="Missing"
	  for (c in st:en){
	    val=UKBBdata[s,c]
	    if (is.na(val)) {
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
if (r==55) { print (output) }
	  Summary[r,s]=output
  }
}

write.table(Summary,OutputFile,sep="\t",col.names=FALSE,row.names=TRUE,quote=FALSE)
}
