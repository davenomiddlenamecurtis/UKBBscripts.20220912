#!/share/apps/R-3.6.1/bin/Rscript

# plink is on the PATH

SNPtemplate="/SAN/ugi/UGIbiobank/data/downloaded/ukb_cal_chr%d_v2"
PCsFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb23158.common.all.20230806.eigenvec.txt"

args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
	args=c("UKBB.HT.txt","22")
}
PhenoFile=args[1]
chr=as.numeric(args[2])
Phenos=data.frame(read.table(PhenoFile,header=FALSE,stringsAsFactors=FALSE,fill=TRUE)) # will work whether or not there is a header
Phenos=Phenos[Phenos[,2]==0 | Phenos[,2]==1 ,]
Phenos[,2]=as.numeric(Phenos[,2])+1
FamFile=sprintf(SNPtemplate,chr)
FamFile=sprintf("%s.fam",FamFile)
Fam=data.frame(read.table(FamFile,header=FALSE,stringsAsFactors=FALSE))
Phenos=Phenos[Phenos[,1] %in% Fam[,1],]
ToWrite=Phenos
ToWrite[,2]=Phenos[,1]
ToWrite[,3]=Phenos[,2] # too lazy to find out how to do this properly
PlinkPhenoFile=sprintf("%s.phenos.%d.txt",PhenoFile,chr) # separate one for each chromosome so parallel jobs do not overwrite each other
write.table(ToWrite,PlinkPhenoFile,col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")

# we will keep SNPs not in HWE

PlinkRoot=sprintf(SNPtemplate,chr)
cmd=sprintf("plink --bfile %s --keep %s --pheno %s --maf 0.05 --covar %s --hide-covar --logistic --out %s.results.%d",
		PlinkRoot,
		PlinkPhenoFile,
		PlinkPhenoFile,
		PCsFile,
		PhenoFile,
		chr)
system(cmd)

# for chr in {1..22}; do subComm.sh /share/apps/R-3.6.1/bin/Rscript ~/UKBB/UKBBscripts.20220912/runGWAS.20240222.R UKBB.HT.txt $chr; done
