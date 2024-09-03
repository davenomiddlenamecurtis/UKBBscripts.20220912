#!/share/apps/R-3.6.1/bin/Rscript

# script to make a QQ plot from GWAS results
# just provide the results file name as an argument

SNPtemplate="/SAN/ugi/UGIbiobank/data/downloaded/ukb_cal_chr%d_v2"

args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
	print("Run with GWAS results file as first argument and phenotype file as second argument, with optional threshold as third argument")
	quit()
}

if (length(args)>2) {
	threshold=as.numeric(args[3])
} else {
	threshold=3
}

ResFile=args[1]
results=data.frame(read.table(ResFile,header=TRUE,,stringsAsFactors=FALSE))
results=results[results$P<10^-threshold,]
results=results[,2,drop=FALSE] # SNP
SNPsFile=sprintf("%s.bestSNPs.txt",ResFile)
write.table(results,SNPsFile,col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")

system("rm mergelist.txt")
for (chr in 1:22) {
	PlinkRoot=sprintf(SNPtemplate,chr)
	PlinkOutRoot=sprintf("%s.bestSNPs.%d",ResFile,chr)
	cmd=sprintf("plink --bfile %s --extract %s --indep-pairwise 50 5 0.5 --out toprune",
		PlinkRoot,SNPsFile)
	system(cmd)
	cmd=sprintf("plink --bfile %s --extract toprune.prune.in --recode A --out %s --make-bed ",
		PlinkRoot,PlinkOutRoot)
	system(cmd)
	if (chr!=1) {
		cmd=sprintf("echo %s >> mergelist.txt",PlinkOutRoot)
		system(cmd)
	}
}

PlinkOutRoot=sprintf("%s.bestSNPs.%d",ResFile,1)

PhenoFile=args[2]
Phenos=data.frame(read.table(PhenoFile,header=FALSE,stringsAsFactors=FALSE,fill=TRUE)) # will work whether or not there is a header
Phenos=Phenos[Phenos[,2]==0 | Phenos[,2]==1 ,]
Phenos[,2]=as.numeric(Phenos[,2])+1
FamFile=sprintf("%s.fam",PlinkOutRoot)
Fam=data.frame(read.table(FamFile,header=FALSE,stringsAsFactors=FALSE))
Phenos=Phenos[Phenos[,1] %in% Fam[,1],]
ToWrite=Phenos
ToWrite[,2]=Phenos[,1]
ToWrite[,3]=Phenos[,2] # too lazy to find out how to do this properly
PlinkPhenoFile=sprintf("%s.phenos.txt",PhenoFile,chr) # separate one for each chromosome so parallel jobs do not overwrite each other
write.table(ToWrite,PlinkPhenoFile,col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")


PlinkOutAllRoot=sprintf("%s.bestSNPs.all",ResFile)
cmd=sprintf("plink --bfile %s --merge-list mergelist.txt --pheno %s --keep %s --out %s --make-bed",PlinkOutRoot,PlinkPhenoFile,PlinkPhenoFile,PlinkOutAllRoot)
system(cmd)

