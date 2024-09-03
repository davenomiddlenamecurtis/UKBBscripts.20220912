#!/share/apps/R-3.6.1/bin/Rscript
# analyse.TCF7L2.20240808.R

# plink is on the PATH

SNPtemplate="/SAN/ugi/UGIbiobank/data/downloaded/ukb_cal_chr%d_v2"
PCsFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb23158.common.all.20230806.eigenvec.txt"

PhenoFile="UKBB.T2D.txt"
margin=100000
MAF=0.05
Gene="TCF7L2"
start=114710006-margin # https://genome.ucsc.edu/cgi-bin/hgGene?hgg_gene=ENST00000355995.9_6&hgg_chrom=chr10&hgg_start=114710005&hgg_end=114927437&hgg_type=knownGene&db=hg19
end=114927437+margin
# GRCh37
chr=10

RawFile=sprintf("%s.results.%d.raw",PhenoFile,chr)
EthnicityFile="UKBB.ethnicity.txt"

if (!file.exists(RawFile)) {
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

Ranges=data.frame(matrix(ncol=4,nrow=1))
Ranges[1,1]=10
Ranges[1,2]=start
Ranges[1,3]=end
Ranges[1,4]=Gene
RangeFile="range.txt"
write.table(Ranges,RangeFile,sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)

PlinkRoot=sprintf(SNPtemplate,chr)
cmd=sprintf("plink --bfile %s --keep %s --pheno %s --maf %f --extract range range.txt --recode A --out %s.results.%d",
		PlinkRoot,
		PlinkPhenoFile,
		PlinkPhenoFile,
		MAF,
		PhenoFile,
		chr)
system(cmd)
} # create RawFile

Ethnicity=data.frame(read.table(EthnicityFile,header=TRUE,fill=TRUE,stringsAsFactors=FALSE))
Ethnicity=Ethnicity[(Ethnicity[,2]==1001 | Ethnicity[,2]==1),]

AllData=data.frame(read.table(RawFile,header=TRUE,stringsAsFactors=FALSE))
AllData=AllData[AllData[,1] %in% Ethnicity[,1],]
PCs=data.frame(read.table(PCsFile,header=TRUE,stringsAsFactors=FALSE))
PCs=PCs[,-1] # avoid duplicating FID
nSNPs=ncol(AllData)-6
AllData=merge(AllData,PCs,by="IID",all=FALSE)
Covars=" + SEX"
for (c in 1:20) {
	Covars=sprintf("%s + PC%d",Covars,c)
}

AllData$PHENOTYPE=AllData$PHENOTYPE-1
ResultsCols=c("SNP","SLP","OR")
Results=list()
for (r in 1:2) {
	Results[[r]]=data.frame(matrix(ncol=length(ResultsCols),nrow=nSNPs))
	colnames(Results[[r]])=ResultsCols
	Results[[r]][,1]=colnames(AllData)[7:(6+nSNPs)]
}
for (v in 1:nSNPs) {
	for (r in 1:2) {
		if (r==1) {
			Formula=sprintf("PHENOTYPE ~ %s",colnames(AllData)[6+v])
		} else {
			Formula=sprintf("PHENOTYPE ~ %s%s",colnames(AllData)[6+v],Covars)
		}
		m=glm(as.formula(Formula), data = AllData, family = "binomial")
		SLP=log10(summary(m)$coefficients[2,4])
		beta=summary(m)$coefficients[2,1]
		if (beta>0) {
			SLP=-SLP
		}
		betaSE=summary(m)$coefficients[2,2]
		OR=sprintf("%.2f (%.2f - %.2f}",exp(beta),exp(beta-2*betaSE),exp(beta+2*betaSE))
		Results[[r]][v,2]=SLP
		Results[[r]][v,3]=OR
		print(Results[[r]][v,])
	}
}

RR=c("uncorrected","corrected")
for (r in 1:2) {
	ResultsFile=sprintf("%s.results.%s,txt",Gene,RR[r])
	write.table(Results[[r]],ResultsFile,col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
}
