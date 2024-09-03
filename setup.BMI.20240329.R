#!/share/apps/R-3.6.1/bin/Rscript

# script to get data files to analyse association with BMI

# note that the column number to provide is one higher than that given in http://www.davecurtis.net/UKBB/ukb41465.html
targetDir="/home/rejudcu/UKBB/BMI.20240329"
BMICol=4039
AgeCol=4047
BMIFile="UKBB.BMI.txt"
AgeFile="UKBB.Age.txt"
setwd(targetDir)
if (!file.exists(BMIFile)) {
	cmd=sprintf("bash /home/rejudcu/UKBB/UKBBscripts.20220912/extract.UKBB.var.41465.20230913.sh BMI %d",BMICol+1)
	system(cmd)
}
if (!file.exists(AgeFile)) {
	cmd=sprintf("bash /home/rejudcu/UKBB/UKBBscripts.20220912/extract.UKBB.var.41465.20230913.sh Age %d",AgeCol+1)
	system(cmd)
}
