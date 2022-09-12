#!/share/apps/R-3.6.1/bin/Rscript

# script to get data files to analyse association of sildenafil with LOAD

UKBBFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb41465.first100rows.tab"
UKBBFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb41465.tab"
drugCodesFile="UKBB.coding4.tsv" # but use this anyway for consistency? not at the moment
sildefanilCodesFile="/home/rejudcu/UKBB/UKBB.sildenafil/sildenafilCodes.txt"
ICD10SourcesFile="ICD10Sources.txt" # sources of ICD10 codes - hospital admissions, causes of death
wd="/home/rejudcu/UKBB/UKBBscripts"

selfReportedCodes=c(1265)
selfReportFirst=1730
selfReportLast=1763

femaleRxCols=c(1472:1475)
maleRxCols=c(1584:1586)

prescribedRxFirst=1866
prescribedRxLast=1913
medicationAtIntake=c("20003-0.0","20003-0.47")

setwd(wd)

medicationFile="/home/rejudcu/UKBB/UKBB.sildenafil/medication.txt"
sildefanilRxFile="/home/rejudcu/UKBB/UKBB.sildenafil/sildefanilRx.txt"

extractCols <- function(destFile,srcFile,first,last) {
	cmd=sprintf("head -n 1 %s > colNames.txt",srcFile)
	system(cmd)
	fields=data.frame(read.table("colNames.txt",header=TRUE))
	f=match(sprintf("f.%s",gsub("-",".",first)),colnames(fields))
	print(f)
	l=match(sprintf("f.%s",gsub("-",".",last)),colnames(fields))
	if (f==l) {
		cmd=sprintf("cut -f 1,%d %s > %s",f,srcFile,destFile)
	} else {
		cmd=sprintf("cut -f 1,%d-%d %s > %s",f,l,srcFile,destFile)
	}
	print(cmd)
	system(cmd)
}

extractCols(medicationFile,UKBBFile,medicationAtIntake[1],medicationAtIntake[2])

sildefanilCodes=data.frame(read.table(sildefanilCodesFile,header=FALSE,sep="\t")) # may be spaces in drug names

prescribedRx=data.frame(read.table(medicationFile,header=TRUE))
sildefanilRx=data.frame(matrix(nrow=nrow(prescribedRx),ncol=2,0))
colnames(sildefanilRx)=c("IID","sildenafil")
sildefanilRx$IID=prescribedRx[,1]
sildefanilRx$sildenafil=0
for (r in 1:nrow(sildefanilCodes)) {
	sildefanilRx$sildenafil[rowSums(prescribedRx==sildefanilCodes[r,1], na.rm = TRUE)>0]=1
}
write.table(sildefanilRx,sildefanilRxFile,quote=FALSE,sep="\t")
