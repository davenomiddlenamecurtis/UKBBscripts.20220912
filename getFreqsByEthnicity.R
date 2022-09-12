ethnicityTable=data.frame(na.omit(read.table("UKBB.ethnicity.txt",header=TRUE,sep="\t")))

genotypeTable=data.frame(na.omit(read.table("genos.positive.rs71325088.raw",header=TRUE,sep=" ")))
table=merge(ethnicityTable,genotypeTable,by="IID")
table[,9]=1
a=aggregate(table[,8],list(table$ethnicity),mean)
summary=a
a=aggregate(table[,9],list(table$ethnicity),sum)
summary[,2]=summary[,2]/2
summary[,3]=a[,2]
colnames(summary)=c("Coding","MAF","N")
ethnicityCodes=data.frame(read.table("/home/rejudcu/UKBB/UKBB.ethnicityCodes.tsv",header=TRUE,sep="\t"))
for (r in 1:nrow(summary)) {
  summary[r,4]=ethnicityCodes[ethnicityCodes$Coding==summary[r,1],2]
}
print(summary)

genotypeTable=data.frame(na.omit(read.table("genos.negative.rs71325088.raw",header=TRUE,sep=" ")))
table=merge(ethnicityTable,genotypeTable,by="IID")
table[,9]=1
a=aggregate(table[,8],list(table$ethnicity),mean)
summary=a
a=aggregate(table[,9],list(table$ethnicity),sum)
summary[,2]=summary[,2]/2
summary[,3]=a[,2]
colnames(summary)=c("Coding","MAF","N")
ethnicityCodes=data.frame(read.table("/home/rejudcu/UKBB/UKBB.ethnicityCodes.tsv",header=TRUE,sep="\t"))
for (r in 1:nrow(summary)) {
  summary[r,4]=ethnicityCodes[ethnicityCodes$Coding==summary[r,1],2]
}
print(summary)
