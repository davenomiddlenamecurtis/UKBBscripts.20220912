#!/share/apps/R-3.6.1/bin/Rscript

# script to obtain counts for different variant types and get effect sizes

CountsModel="UKBB.DEMScore.counts.20241031"
CountsArgFile=sprintf("gva.%s.arg",CountsModel)
bashScript="run.DEMScore.PlusExWAS.counts.20241109.sh"
ResultsFile="DEMScore.PlusExWAS.Results.txt"
DefaultTestFile="~/UKBB/RAPfiles/pars/test.counts.withAPOE.20241101.tst"

DxCommand="dx cd / ; dx run swiss-army-knife -y  --ignore-reuse --instance-type mem3_hdd2_v2_x4 -imount_inputs=FALSE -iin=/scripts/runOnePlinkGeneWithFileList.20241101.sh "

wd="/home/rejudcu/UKBB/dementia.2024"
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
	test=gsub(".withAPOE.tst","",test)
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
		TestFile=sprintf("test.counts.withAPOE.for.%s.tst",gene)
		CommStr=sprintf("cp %s %s",DefaultTestFile,TestFile)
		system(CommStr)
		write(sprintf("%s 0.0 0 1",test),TestFile,append=TRUE)
		CommStr=sprintf("dx cd /pars; dxupload %s",TestFile)
		system(CommStr)
		write(sprintf("/pars/%s",TestFile),FileList,append=TRUE)
		ArgFile=sprintf("gva.%s.%s.arg",CountsModel,test)
		CommStr=sprintf("cp ~/UKBB/RAPfiles/pars/%s %s",CountsArgFile,ArgFile)
		system(CommStr)
		write(sprintf("--lintestfile %s",TestFile),ArgFile,append=TRUE)
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
	Results[r,2]=gsub(".withAPOE.tst","",Results[r,2])
	test=Results[r,2]
	SaoFile=sprintf("%s.%s.sao",CountsModel,gene)
	SaoTable=na.omit(data.frame(read.table(SaoFile,header=TRUE,stringsAsFactors=FALSE,fill=TRUE)))
	colnames(SaoTable)[14:22]=MainWeights
	colnames(SaoTable)[23]="AM_prediction"
	colnames(SaoTable)[24]="AM_score"
	colnames(SaoTable)[25:(25+length(ExtraWeights)-1)]=ExtraWeights
	GlmTable=data.frame(read.table(SaoFile,header=TRUE,stringsAsFactors=FALSE,fill=TRUE))
	First=which("L1"==GlmTable[,1])[[1]]
	GlmTable=GlmTable[-First:-1,1:4]
	colnames(GlmTable)=GlmTable[1,]
	Cats=c(MainWeights[-1],test)
	GeneResults=data.frame(matrix(ncol=6,nrow=length(Cats)))
	colnames(GeneResults)=c("Category","NumVars","Carrier count","Carrier mean (SE)","Beta (SE)","Beta/SE")
	for (c in 1:length(Cats)) {
		Cat=Cats[c]
		GeneResults[c,1]=Cat
		Tab=SaoTable[as.numeric(SaoTable[,Cat])!=0,]
		Tab$contAB=as.numeric(Tab$contAB)
		Tab$meanAB=as.numeric(Tab$meanAB)
		GeneResults[c,2]=nrow(Tab)
		GeneResults[c,3]=sum(Tab$contAB)
		Tab$SigmaX=Tab$contAB*Tab$meanAB
		Tab$SigmaX2=Tab$contAB*Tab$meanAB*Tab$meanAB
		Mean=sum(Tab$SigmaX)/sum(Tab$contAB)
		Var=(sum(Tab$SigmaX2)-sum(Tab$contAB)*Mean*Mean)/(sum(Tab$contAB)-1)
		GeneResults[c,4]=sprintf("%.3f (%.3f)",Mean,sqrt(Var)/sqrt(sum(Tab$contAB)))
		gr=which(Cat==GlmTable$beta)[[1]]
		GeneResults[c,5]=sprintf("%.3f (%.3f)",as.numeric(GlmTable$value[gr]),as.numeric(GlmTable$SE[gr]))
		GeneResults[c,6]=sprintf("%.2f",as.numeric(GlmTable$z[gr]))
	}
	print(GeneResults)
	GeneResults[1:length(CategoryNames),1]=CategoryNames
	write.table(GeneResults,sprintf("%s.counts.txt",gene),col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
}

