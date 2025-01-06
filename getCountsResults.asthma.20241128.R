#!/share/apps/R-3.6.1/bin/Rscript

# script to obtain counts for different variant types and get effect sizes

CountsModel="UKBB.asthma.counts.20241128"
CountsArgFile=sprintf("gva.%s.arg",CountsModel)
bashScript="run.asthma.counts.20241128.sh"
ResultsFile="Asthma.Results.20241118.txt"
DefaultTestFile="~/UKBB/RAPfiles/pars/test.counts.20241128.tst"

DxCommand="dx cd / ; dx run swiss-army-knife -y  --ignore-reuse --instance-type mem3_hdd2_v2_x4 -imount_inputs=FALSE -iin=/scripts/runOnePlinkGeneWithFileList.20241101.sh "

wd="/home/rejudcu/UKBB/asthma.2024"
setwd(wd)

MainWeights=c("VEPWeight","FivePrime","InDelEtc","IntronicEtc","LOF","ProteinAltering","SpliceRegion","Synonymous","ThreePrime")
ExtraWeightsFile="~/UKBB/RAPfiles/annot/extraWeights.20231106.txt"
CategoryNames=c("Five prime UTR","InDel etc","Intronic etc","LOF","Protein altering","Splice region","Synonymous","Three prime UTR")

ExtraWeights=data.frame(read.table(ExtraWeightsFile,header=TRUE,stringsAsFactors=FALSE,fill=TRUE))[,1]

Results=na.omit(data.frame(read.table(ResultsFile,header=TRUE,stringsAsFactors=FALSE,fill=TRUE)))
wd="genes"
setwd(wd)


# need to run  ~/UKBB/RAPfiles/RAPscripts/makePlinkGeneWESVCF.20241025.sh first on all neeeded genes so there is a gene.vcf.gz file

AllDone=TRUE
for (r in 1:nrow(Results)) {
	gene=Results[r,1]
	test=Results[r,2]
	test=gsub(".tst","",test)
	test=gsub("test.annot.","",test)
	SaoFile=sprintf("%s.%s.sao",CountsModel,gene)
	if (!file.exists(SaoFile)){
		commStr=sprintf("dx cd /results/%s.%s/geneResults;dx download %s",CountsModel,test,SaoFile)
		print(commStr)
		system(commStr)
	}	
	if (!file.exists(SaoFile)){
		AllDone=FALSE
		FileList=sprintf("for.%s.%s.lst",CountsModel,test)
		CommStr=sprintf("cp ~/UKBB/RAPfiles/fileLists/for.%s.lst %s",CountsModel,FileList)
		system(CommStr)
		TestFile=sprintf("test.counts.for.%s.tst",gene)
		CommStr=sprintf("cp %s %s",DefaultTestFile,TestFile)
		system(CommStr)
		write(sprintf("%s 0.0 0 1",test),TestFile,append=TRUE)
		CommStr=sprintf("dx cd /pars; dxupload %s",TestFile)
		system(CommStr)
		write(sprintf("/pars/%s",TestFile),FileList,append=TRUE)
		ArgFile=sprintf("gva.%s.%s.arg",CountsModel,test)
		CommStr=sprintf("cp ~/UKBB/RAPfiles/pars/%s %s",CountsArgFile,ArgFile)
		system(CommStr)
		write(sprintf("--testfile %s",TestFile),ArgFile,append=TRUE)
		CommStr=sprintf("dx cd /pars; dxupload %s",ArgFile)
		system(CommStr)
		CommStr=sprintf("dx cd /fileLists; dxupload %s",FileList)
		system(CommStr)
		CommStr=sprintf("dx mkdir /results/%s.%s",CountsModel,test)
		cmd=sprintf("bash runOnePlinkGeneWithFileList.20241101.sh %s.%s %s for.%s.%s.lst",
			CountsModel,test,gene,CountsModel,test)
		CommStr=sprintf("%s -icmd=\"%s\"",DxCommand,cmd)
		print(CommStr)
		system(CommStr)
	}
}

if (AllDone==FALSE) {
	quit()
}

for (r in 1:nrow(Results)) {
	gene=Results[r,1]
	Results[r,2]=gsub("test.annot.","",Results[r,2])
	Results[r,2]=gsub(".tst","",Results[r,2])
	test=Results[r,2]
	SaoFile=sprintf("%s.%s.sao",CountsModel,gene)
	SaoTable=na.omit(data.frame(read.table(SaoFile,header=TRUE,stringsAsFactors=FALSE,fill=TRUE)))
	colnames(SaoTable)[16:24]=MainWeights
	colnames(SaoTable)[25]="AM_prediction"
	colnames(SaoTable)[26]="AM_score"
	colnames(SaoTable)[27:(27+length(ExtraWeights)-1)]=ExtraWeights
	GlmTable=data.frame(read.table(SaoFile,header=TRUE,stringsAsFactors=FALSE,fill=TRUE))
	First=which("L1"==GlmTable[,1])[[1]]
	GlmTable=GlmTable[-First:-1,1:4]
	colnames(GlmTable)=GlmTable[1,]
	Cats=c(MainWeights[-1],test)
	GeneResults=data.frame(matrix(ncol=8,nrow=length(Cats),0))
	colnames(GeneResults)=c("Category","NumVars","TotCountControls","MeanCountControls","TotCountCases","MeanCountCases","OR","SLP")
	for (c in 1:length(Cats)) {
		Cat=Cats[c]
		GeneResults[c,1]=Cat
		Tab=SaoTable[as.numeric(SaoTable[,Cat])!=0,]
		for (cc in c(2,4,6,8,10,12)) {
			Tab[,cc]=as.numeric(Tab[,cc])
		}
		GeneResults$NumVars[c]=nrow(Tab)
		if (nrow(Tab)>0) {
			GeneResults$TotCountControls[c]=sum(Tab$contAB)+2*sum(Tab$contBB)
			GeneResults$MeanCountControls[c]=GeneResults$TotCountControls[c]/((sum(Tab$contAA)+sum(Tab$contAB)+sum(Tab$contBB))/nrow(Tab))
			GeneResults$TotCountCases[c]=sum(Tab$caseAB)+2*sum(Tab$caseBB)
			GeneResults$MeanCountCases[c]=GeneResults$TotCountCases[c]/((sum(Tab$caseAA)+sum(Tab$caseAB)+sum(Tab$caseBB))/nrow(Tab))
			gr=which(Cat==GlmTable$beta)[[1]]
			b=as.numeric(GlmTable$value[gr])
			SE=as.numeric(GlmTable$SE[gr])
			GeneResults$OR[c]=sprintf("%.2f (%.2f-%.2f)",exp(b),exp(b-2*SE),exp(b+2*SE))
			GeneResults$SLP[c]=log10(2*pnorm(abs(as.numeric(GlmTable$z[gr])),lower.tail=FALSE))
			if (b>0) {
				GeneResults$SLP[c]=GeneResults$SLP[c]* -1
			}
		}
	}
	print(GeneResults)
	GeneResults[1:length(CategoryNames),1]=CategoryNames
	write.table(GeneResults,sprintf("%s.counts.txt",gene),col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
}

