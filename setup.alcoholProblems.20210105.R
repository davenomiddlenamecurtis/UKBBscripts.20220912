#!/share/apps/R-3.6.1/bin/Rscript

# script to get data files to analyse association with high alcohol consumption

# note that the column number to provide is one higher than that given in http://www.davecurtis.net/UKBB/ukb41465.html

exomesFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb41465.exomes.20201103.txt"
wd="/home/rejudcu/UKBB/UKBBscripts"
# need to add 1 to these because first column is numbered 0
dateCol=3854 # people who answered mental health questionnaire
yesNoCols=c(3857,3858,3859) # >0
monthlyCols=c(3860,3861,3862,3865,3866) # >=3

setwd(wd)

cmd=sprintf("cut %s -f 1,%d",exomesFile,dateCol+1)
for (c in 1:length(yesNoCols)) {
  cmd=sprintf("%s,%d",cmd,yesNoCols[c]+1)
}
for (c in 1:length(monthlyCols)) {
  cmd=sprintf("%s,%d",cmd,monthlyCols[c]+1)
}
cmd=sprintf("%s > ../alcoholProblems/alcoholProblemsData.txt",cmd)
system(cmd)

alcoholProblems=data.frame(read.table("../alcoholProblems/alcoholProblemsData.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE))
alcoholProblems=alcoholProblems[alcoholProblems[,2]!="",] # only use those with a date
toWrite=data.frame(matrix(ncol=2,nrow=nrow(alcoholProblems)))
colnames(toWrite)=c("IID","alcoholProblems")
toWrite$IID=alcoholProblems[,1]
toWrite$alcoholProblems=0
for (c in 3:(length(yesNoCols)+2)) {
  toWrite$alcoholProblems[alcoholProblems[,c]>=1]=1
}
for (c in (3+length(yesNoCols)):(length(yesNoCols)+2+length(monthlyCols))) {
  toWrite$alcoholProblems[alcoholProblems[,c]>=3]=1
}

write.table(toWrite,"../alcoholProblems/alcoholProblems.20210105.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
