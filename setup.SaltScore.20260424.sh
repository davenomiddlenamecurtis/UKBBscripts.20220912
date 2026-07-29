#!/share/apps/R-3.6.1/bin/Rscript

# script to get data files to analyse association with BMI

# note that the column number to provide is one higher than that given in http://www.davecurtis.net/UKBB/ukb41465.html
targetDir="/home/rejudcu/UKBB/SaltScore.2026"
RawSaltCol=507
RawSaltFile="UKBB.RawSaltScore.txt"
SaltFile="UKBB.SaltScore.txt"
setwd(targetDir)
if (!file.exists(RawSaltFile)) {
	cmd=sprintf("bash /home/rejudcu/UKBB/UKBBscripts.20220912/extract.UKBB.var.41465.20260424.sh RawSaltScore %d",RawSaltCol+1)
	system(cmd)
}
if (!file.exists(SaltFile)) {
	SaltScore=na.omit(data.frame(read.table(RawSaltFile,header=TRUE,stringsAsFactors=FALSE,fill=TRUE)))
	SaltScore=SaltScore[SaltScore$RawSaltScore>=1 & SaltScore$RawSaltScore<=4,]
	colnames(SaltScore)[2]="SaltScore"
	write.table(SaltScore,SaltFile,sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
}
