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
# setwd(wd)
SummarySchemeFile="SummaryScheme.20210618.txt"
SummaryScheme=data.frame(read.table(SummarySchemeFile,header=TRUE,sep="\t"))
afterIGVFile="afterIGV.20210702.txt"
afterIGV=data.frame(read.table(afterIGVFile,header=TRUE,sep="\t",stringsAsFactors=FALSE))
badIGVs=afterIGV[afterIGV$Keep=="N",1]

largestCoding=200000
codings= vector("list", largestCoding)

for (gene in genes) {
annotsFile=sprintf("IDsAndCounts.%s.%s.txt",model,gene)
UKBBdataFile=sprintf("%s.%s.phenos.txt",model,gene)
goodUKBBdataFile=sprintf("%s.%s.phenos.good.txt",model,gene)
OutputFile=sprintf("SummaryFor.%s.txt",gene)

if (!file.exists(UKBBdataFile)) {
	next
}

UKBBdata=data.frame(fread(UKBBdataFile))
annots=data.frame(read.table(annotsFile,header=FALSE,sep="\t",stringsAsFactors=FALSE))

# now do some QC
minDepth=20
minAB=0.25 # AB must exceed this
minFisherP=0.01
nBad=0
badAnnots=annots[0,]
goodOnes=vector()
for (r in 1:nrow(annots)) {
	bad=FALSE
	deepestReads=c(0,0)
	nextDeepestReads=c(0,0)
	totalDepth=0
	counts=strsplit(annots[r,3],"\\s+")[[1]]
	for (c in 1:length(counts)) {
		reads=as.numeric(strsplit(strsplit(counts[c],":")[[1]][2],"/")[[1]])
		totalDepth=totalDepth+reads[1]+reads[2]
		if (reads[1]+reads[2] > deepestReads[1]+deepestReads[2]) {
			nextDeepestReads=deepestReads
			deepestReads=reads
		} else if (reads[1]+reads[2] > nextDeepestReads[1]+nextDeepestReads[2]) {
			nextDeepestReads=reads
		}
	}
	if (totalDepth<minDepth){
		bad=TRUE
	}
	if ((deepestReads[1]+deepestReads[2])/totalDepth<=minAB) {
		bad=TRUE
	}
	if ((nextDeepestReads[1]+nextDeepestReads[2])/totalDepth<=minAB) {
		bad=TRUE
	}
	dat=data.frame(deepestReads,nextDeepestReads)
	if (fisher.test(dat)$p.value<minFisherP) {
		bad=TRUE
	}
	if (bad) {
		nBad=nBad+1
		badAnnots[nBad,]=annots[r,]
		goodOnes[r]=FALSE
	} else {
		goodOnes[r]=TRUE
	}
	secondColon=unlist(gregexpr(pattern=":",annots[r,4]))[2]
	variant=substr(annots[r,4],1,secondColon-1)
	transcripts=strsplit(substr(annots[r,4],secondColon+1,nchar(annots[r,4])),",")[[1]]
	transcripts=grep(gene,transcripts,value=TRUE)
	annots[r,4]=paste(variant,paste(transcripts,collapse=" "),sep=" ")
}

annots=annots[goodOnes,]
write.table(annots,sprintf("IDsAndCounts.good.%s.%s.txt",model,gene),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
write.table(badAnnots,sprintf("IDsAndCounts.bad.%s.%s.txt",model,gene),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)

UKBBdata=UKBBdata[UKBBdata[,1]  %in% annots[,1],]

UKBBdata=UKBBdata[!(UKBBdata[,1]  %in% badIGVs),] # remove ones with unconvincing IGV reads
write.table(UKBBdata,goodUKBBdataFile,sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
annots=annots[!(annots[,1]  %in% badIGVs),] 

# adjust the full variant annotation to only include those for this gene
# maybe try to find the field which identifies the canonical transcript
# no, better to rebuild the whole annotation just without those for other genes

SummaryTab=data.frame(matrix(ncol=nrow(UKBBdata),nrow=nrow(SummaryScheme)))
rownames(SummaryTab)=SummaryScheme$Name
SummaryTab[2,]=annots[,4] 
SummaryTab[3,]=annots[,3] 
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
}
