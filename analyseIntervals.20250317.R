#!/share/apps/R-3.6.1/bin/Rscript

# script to run intVarAssoc on a number of intervals

# Tasks
# Create a VCF for each interval
# Produce annotations with GPN for that interval and upload them
# Run intVarAssoc on each interval 

wd="/home/rejudcu/UKBB/migraine.2025"
Model="UKBB.migraine.nonCoding.20250317"
IntervalsFile="migraine.regRegions.20250317.txt"

wd="/home/rejudcu/UKBB/LDLRpromoter.20250814"
Model="UKBB.HL.nonCoding.20250814"
IntervalsFile="LDLRpromoter.20250814.txt"


setwd(wd)

args = commandArgs(trailingOnly=TRUE)
if (length(args)>=2) {
	Model=args[1]
	IntervalsFile=args[2]
}

GPNScoreFile="/cluster/ref1/UGISharedData/GPN-MSA/scores.tsv.bgz"

Intervals=data.frame(read.table(IntervalsFile,header=FALSE,stringsAsFactors=FALSE))

AllDone=TRUE
for (r in 1:nrow(Intervals)) {
	Name=Intervals[r,1]
	VCF=sprintf("%s.vcf.gz",Name)
	system(sprintf("rm %s.exists.txt",VCF))
	system(sprintf("dx cd /WGS/VCFs;f=`dx ls %s`; if [ .$f == .%s ]; then echo got it > %s.exists.txt; fi",VCF,VCF,VCF))
	if (!file.exists(sprintf("%s.exists.txt",VCF))) {
		AllDone=FALSE
		system(sprintf("bash /home/rejudcu/UKBB/RAPfiles/RAPscripts/makeExtractedVCF.20240926.sh %s chr%s %d %d",Intervals[r,1],Intervals[r,2],Intervals[r,3],Intervals[r,4]))
		ScriptName=sprintf("extract.%s.sh",Name)
		system(sprintf("dx cd /scripts ; dxupload %s",ScriptName))
		CommStr=sprintf("dx cd / ; dx run swiss-army-knife -y --ignore-reuse --instance-type mem3_hdd2_v2_x4  -imount_inputs=FALSE -iin=/scripts/%s -icmd=\"bash %s\"",ScriptName,ScriptName)
		print(CommStr)
		system(CommStr)
	}
}

if (AllDone!=TRUE) {
	quit()
}

for (r in 1:nrow(Intervals)) {
	Name=Intervals[r,1]
	AnnotVCF=sprintf("%s.GPN.annot.vcf.gz",Name)
	system(sprintf("rm %s.exists.txt",AnnotVCF))
	system(sprintf("dx cd /annot;f=`dx ls %s`; if [ .$f == .%s ]; then echo got it > %s.exists.txt; fi",AnnotVCF,AnnotVCF,AnnotVCF))
	if (!file.exists(sprintf("%s.exists.txt",AnnotVCF))) {
		AllDone=FALSE
		system(sprintf("rm %s*",AnnotVCF))
		system(sprintf("tabix %s %s:%d-%d > %s.GPN.annot.vcf",GPNScoreFile,Intervals[r,2],Intervals[r,3],Intervals[r,4],Name))
		system(sprintf("bgzip %s.GPN.annot.vcf",Name))
		system(sprintf("tabix -p vcf %s",AnnotVCF))
		system(sprintf("dx cd /annot; dxupload %s",AnnotVCF))
		system(sprintf("dx cd /annot; dxupload %s.tbi",AnnotVCF))
	}
}

for (r in 1:nrow(Intervals)) {
	Name=Intervals[r,1]
	IntFile=sprintf("%s.int",Name)
	cat(sprintf("%s:%s-%s\n",Intervals[r,2],Intervals[r,3],Intervals[r,4]),file=IntFile,append=FALSE)
	system(sprintf("dx cd /intervals; dxupload %s",IntFile))
}

AllDone=TRUE
for (r in 1:nrow(Intervals)) {
	Name=Intervals[r,1]
	SaoFile=sprintf("%s.%s.sao",Model,Name)
	if (!file.exists(sprintf("/cluster/project9/bipolargenomes/UKBB/%s/results/%s",Model,SaoFile))) {
		CommStr=sprintf("dx cd /results/%s/intervalResults/;pushd /cluster/project9/bipolargenomes/UKBB/%s/results;dx download %s;popd",Model,Model,SaoFile)
		print(CommStr)
		system(CommStr)
		if (!file.exists(sprintf("/cluster/project9/bipolargenomes/UKBB/%s/results/%s",Model,SaoFile))) {
			AllDone=FALSE
			CommStr=sprintf("cp /home/rejudcu/UKBB/RAPfiles/fileLists/for.%s.lst for.%s.lst; echo /reference38/chr%s.fa>> for.%s.lst;dx cd /fileLists; dxupload for.%s.lst",Model,Name,Intervals[r,2],Name,Name)
			print(CommStr)
			system(CommStr)
			CommStr=sprintf("dx cd /;dx run swiss-army-knife -y --ignore-reuse --instance-type mem3_ssd2_v2_x4 -imount_inputs=FALSE -iin=\"/scripts/runOneIntervalWithFileList.20250418.sh\" -icmd=\"bash runOneIntervalWithFileList.20250418.sh %s %s for.%s.lst\" ",Model,Name,Name)
			print(CommStr)
			system(CommStr)
		}
	}
}

if (AllDone!=TRUE) {
	quit()
}




Summary$MaxMLP=apply(Summary[,2:ncol(Summary)], 1, max)# create new column
Summary=Summary[order(Summary$MaxMLP,decreasing=TRUE),]
Top=Summary[1:NumBest,]
TopGenes=Top[,1]
write.table(TopGenes,BestGenesFile,col.names=FALSE,row.names=FALSE,quote=FALSE)

AllDone=TRUE
for (gene in Top[,1]) {
	SaoFile=sprintf("/cluster/project9/bipolargenomes/UKBB/%s/results/%s.%s.sao",BestModel,BestModel,gene)
	if (!file.exists(SaoFile)) {
		AllDone=FALSE
		break
	}
}

if (!AllDone) {
	cmd=sprintf("bash /home/rejudcu/UKBB/RAPfiles/scripts/getAllSaos.20240815.sh %s %s for.%s.lst",BestModel,BestGenesFile,BestModel)
	print("Command to run analyses on RAP:")
	print(cmd)
	print("Should be repeated until all analyses have run and results downloaded")
	system(cmd)
	quit()
}

# only get to here if all new analyses ran and got downloaded OK, so no else statement

# there are only 100 rows so we will just do a loop so we know what we are doing
# will need to set up a single scoreassoc analysis for each gene

Tests=data.frame(read.table(TestsFile,header=FALSE,stringsAsFactors=FALSE))

TopTests=data.frame(matrix("",nrow=nrow(Top),ncol=3),stringsAsFactors=FALSE)
colnames(TopTests)=c("Gene","TestFile","MaxMLP")
TopTests[,1]=Top[,1]
for (r in 1:nrow(Top)) {
	for (t in 1:nrow(Tests)) {
		if (Top[r,t+1]==Top$MaxMLP[r]) {
			TopTests$TestFile[r]=Tests[t,1]
			TopTests$MaxMLP[r]=Top$MaxMLP[r]
			break
		}
	}
}

write.table(TopTests,TopTestsFile,col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

TopTests=data.frame(read.table(TopTestsFile,header=TRUE,sep="\t",stringsAsFactors=FALSE))

AllDone=TRUE
commStr=sprintf("dx cd /results/%s",NewModel)
print(commStr)
system(commStr)

for (gene in Top[,1]) {
	SummFile=sprintf("/home/rejudcu/UKBB/dementia.2024/genes/summ.UKBB.dementia.best.annot.20240620.%s.txt",gene)
	if (!file.exists(SummFile)) {
		commStr=sprintf("pushd /home/rejudcu/UKBB/dementia.2024/genes; dx download summ.UKBB.dementia.best.annot.20240620.%s.txt; popd",gene)
		print(commStr)
		system(commStr)
		if (!file.exists(SummFile)) {
			AllDone=FALSE
		}
	}
}

if (!AllDone) {
GenesToDo=TopTests
sink(bashScript)
for (r in 1:nrow(TopTests)) {
	gene=TopTests[r,1]
	GenesToDo$DoThis[r]=FALSE
	SummFile=sprintf("/home/rejudcu/UKBB/dementia.2024/genes/summ.UKBB.dementia.best.annot.20240620.%s.txt",gene)
	if (file.exists(SummFile)) {
		next
	}
	GenesToDo$DoThis[r]=TRUE
	cat(sprintf("Rscript scoreassoc.R --arg-file rsco.%s.rarg --gene %s --lintestfile %s --summaryoutputfile summ.UKBB.dementia.best.annot.20240620.%s.txt\n",
		NewModel,gene,TopTests[r,2],gene))
}
sink()
GenesToDo=GenesToDo[GenesToDo$DoThis==TRUE,1]
write.table(GenesToDo,ToDoGenesFile,col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
commStr=sprintf("dx cd /geneLists; dxupload %s",ToDoGenesFile)
print(commStr)
system(commStr)
commStr=sprintf("dx cd /scripts; dxupload %s",bashScript)
print(commStr)
system(commStr)
cmd=sprintf("bash runMultipleScoreassoc.20240822.sh %s for.%s.lst %s %s",
		NewModel,NewModel,bashScript,ToDoGenesFile)
commStr=sprintf("%s -icmd=\"%s\"",DxCommand,cmd)
print(commStr)
system(commStr)
quit()
}

SummaryTable=data.frame(read.table(TopTestsFile,header=TRUE,sep="\t",stringsAsFactors=FALSE))
SummaryTable$TestFile=gsub("/home/rejudcu/pars/","",SummaryTable$TestFile)
SummaryTable$MLP270K=0
for (r in 1:nrow(TopTests)) {
	gene=TopTests[r,1]
	SummFile=sprintf("/home/rejudcu/UKBB/dementia.2024/genes/summ.UKBB.dementia.best.annot.20240620.%s.txt",gene)
	Summ=read.table(SummFile,header=TRUE,sep="\t",stringsAsFactors=FALSE)
	SummaryTable$MLP270K[r]=Summ[1,2]
}
write.table(SummaryTable,ResultsFile,col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

SummaryTable=data.frame(read.table(ResultsFile,header=TRUE,sep="\t",stringsAsFactors=FALSE))
GenesToDo=SummaryTable
GenesToDo$DoThis=SummaryTable$MLP270K>MLPThreshold
GenesToDo=GenesToDo[GenesToDo$DoThis,]

AllDone=TRUE
for (r in 1:nrow(GenesToDo)) {
	gene=GenesToDo[r,1]
	SummFile=sprintf("/home/rejudcu/UKBB/dementia.2024/genes/summ.UKBB.DEMScore.470K.annot.20240904.%s.txt",gene)
	if (!file.exists(SummFile)) {
		commStr=sprintf("pushd /home/rejudcu/UKBB/dementia.2024/genes; dx download summ.UKBB.DEMScore.470K.annot.20240904.%s.txt; popd",gene)
		print(commStr)
		system(commStr)
		if (!file.exists(SummFile)) {
			AllDone=FALSE
		}
	}
}

if (!AllDone) {
sink(bashScript470K)
for (r in 1:nrow(TopTests)) {
	gene=GenesToDo[r,1]
	GenesToDo$DoThis[r]=FALSE
	SummFile=sprintf("/home/rejudcu/UKBB/dementia.2024/genes/summ.UKBB.DEMScore.470K.annot.20240904.%s.txt",gene)
	if (file.exists(SummFile)) {
		next
	}
	GenesToDo$DoThis[r]=TRUE
	cat(sprintf("Rscript scoreassoc.R --arg-file rsco.%s.rarg --gene %s --lintestfile %s --summaryoutputfile summ.UKBB.dementia.best.annot.20240620.%s.txt\n",
		New470KModel,gene,GenesToDo[r,2],gene))
}
sink()
GenesToDo=GenesToDo[GenesToDo$DoThis==TRUE,1]
write.table(GenesToDo,ToDoGenesFile,col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
commStr=sprintf("dx cd /geneLists; dxupload %s",ToDoGenesFile)
print(commStr)
system(commStr)
commStr=sprintf("dx cd /scripts; dxupload %s",bashScript470K)
print(commStr)
system(commStr)
cmd=sprintf("bash runMultipleScoreassoc.20240822.sh %s for.%s.lst %s %s",
		New470KModel,New470KModel,bashScript470K,ToDoGenesFile)
commStr=sprintf("%s -icmd=\"%s\"",DxCommand,cmd)
print(commStr)
system(commStr)
quit()
}

SummaryTable=data.frame(read.table(ResultsFile,header=TRUE,sep="\t",stringsAsFactors=FALSE))
SummaryTable$MLP470K=""
for (r in 1:nrow(SummaryTable)) {
	if (SummaryTable$MLP270K[r]>MLPThreshold) {
		gene=SummaryTable$Gene[r]
		SummFile=sprintf("/home/rejudcu/UKBB/dementia.2024/genes/summ.UKBB.DEMScore.470K.annot.20240904.%s.txt",gene)
		Summ=read.table(SummFile,header=TRUE,sep="\t",stringsAsFactors=FALSE)
		SummaryTable$MLP470K[r]=Summ[1,2]
	}
}
write.table(SummaryTable,ResultsFile,col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")






