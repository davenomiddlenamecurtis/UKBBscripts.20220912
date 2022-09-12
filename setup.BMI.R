#!/share/apps/R-3.6.1/bin/Rscript

# script to get data files to analyse association with BMI

# note that the column number to provide is one higher than that given in http://www.davecurtis.net/UKBB/ukb41465.html
targetDir="/home/rejudcu/UKBB/BMI.20210111"
BMICol=4039
AgeCol=4047
BMIFile="UKBB.BMI.20201103.txt"
AgeFile="UKBB.Age.20201103.txt"
SexFile="/home/rejudcu/UKBB/UKBB.sex.20201111.txt"
setwd(targetDir)
if (!file.exists(BMIFile)) {
	cmd=sprintf("bash /home/rejudcu/UKBB/UKBBscripts/extract.UKBB.var.exomes.20201207.sh BMI %d",BMICol+1)
	system(cmd)
}
if (!file.exists(AgeFile)) {
	cmd=sprintf("bash /home/rejudcu/UKBB/UKBBscripts/extract.UKBB.var.exomes.20201207.sh Age %d",AgeCol+1)
	system(cmd)
}

BMI=na.omit(data.frame(read.table(BMIFile,header=TRUE,stringsAsFactors=FALSE,fill=TRUE)))
Age=na.omit(data.frame(read.table(AgeFile,header=TRUE,stringsAsFactors=FALSE,fill=TRUE)))
Sex=na.omit(data.frame(read.table(SexFile,header=TRUE,stringsAsFactors=FALSE,fill=TRUE)))
allData=merge(BMI,Age,by="IID",,all=FALSE)
allData=merge(allData,Sex,by="IID",,all=FALSE)
# allData=!is.na(allData)
sexString=c("male","female")
output=""
for (s in 1:2) {
	toUse=allData[allData$Sex==s-1,]
	output=sprintf("%sThere were %d %s subjects with mean age %.1f (SD=%.1f) and mean BMI %.1f (SD=%.1f). ",output,nrow(toUse),sexString[s],mean(toUse$Age),sd(toUse$Age),mean(toUse$BMI),sd(toUse$BMI))
}
print(output)
