#!/share/apps/R-3.6.1/bin/Rscript

# accept file name as an argument and append a column indicating ancestry
# assume IDs are in first column and file is white space separated

args = commandArgs(trailingOnly=TRUE)
Data=data.frame(read.table(args[1],header=FALSE,sep="",stringsAsFactors=FALSE))
colnames(Data)[1]="IID"
ethnicityCodes=data.frame(read.table("/home/rejudcu/UKBB/UKBB.ethnicityCodes.tsv",header=TRUE,sep="\t"))
ethnicityTable=data.frame(read.table("/home/rejudcu/UKBB/UKBB.ethnicity.txt",header=TRUE,sep="\t"))
colnames(ethnicityTable)=c("IID","ethnicityCode")
ethnicityCodes$Meaning=as.character(ethnicityCodes$Meaning)
for (r in 1:nrow(ethnicityCodes)) {
  rows=ethnicityTable$ethnicityCode==ethnicityCodes$Coding[r]
  ethnicityTable$Ethnicity[rows]=ethnicityCodes$Meaning[r]
  }

Data=merge(Data,ethnicityTable,by="IID")
write.table(Data,sprintf("%s.anc.txt",args[1]),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)



