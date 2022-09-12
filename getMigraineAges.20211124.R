#!/share/apps/R-3.6.1/bin/Rscript

# script to get data files to analyse association with Migraine

# note that the column number to provide is one higher than that given in http://www.davecurtis.net/UKBB/ukb41465.html
targetDir="/home/rejudcu/UKBB/migraine.20211015"
AgeCol=4047
MigraineFile="/home/rejudcu/UKBB/migraine.20211015/UKBB.migraine.txt"
AgeFile="UKBB.Age.20201103.txt"
SexFile="/home/rejudcu/UKBB/UKBB.sex.20201111.txt"
setwd(targetDir)
if (!file.exists(AgeFile)) {
	cmd=sprintf("bash /home/rejudcu/UKBB/UKBBscripts/extract.UKBB.var.exomes.20201207.sh Age %d",AgeCol+1)
	system(cmd)
}

Migraine=na.omit(data.frame(read.table(MigraineFile,header=TRUE,stringsAsFactors=FALSE,fill=TRUE)))
Age=na.omit(data.frame(read.table(AgeFile,header=TRUE,stringsAsFactors=FALSE,fill=TRUE)))
Sex=na.omit(data.frame(read.table(SexFile,header=TRUE,stringsAsFactors=FALSE,fill=TRUE)))
allData=merge(Migraine,Age,by="IID",,all=FALSE)
allData=merge(allData,Sex,by="IID",,all=FALSE)
# allData=!is.na(allData)
ccString=c("controls","cases")
output=""
for (s in 1:2) {
	toUse=allData[allData$migraine==s-1,]
	output=sprintf("%s%s of whom %.1f%% were female with mean age %.1f (SD=%.1f) ",
		output,ccString[s],sum(toUse$Sex==0)/nrow(toUse)*100,mean(toUse$Age),sd(toUse$Age))
}
print(output)
