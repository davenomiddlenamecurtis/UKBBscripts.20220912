#!/share/apps/R-3.6.1/bin/Rscript

# script to take summary file, select most significant gene/annotation results
# then run scoreassoc on second dataset with those pairs

NumBest=100
SummFile="/cluster/project9/bipolargenomes/UKBB/UKBB.DEMScore.annot.20240903/UKBB.DEMScore.annot.20240903.summ.txt"
BestModel="UKBB.DEMScore.forAnnot.20240904"
BestGenesFile="DEMScore.PlusExWAS.20241109.bestGenes.txt"
TestsFile="/home/rejudcu/pars/UKBB.annot.allTests.withAPOE.20240620.txt"
TopTestsFile="DemScore.topTests.PlusExWAS.20241109.txt"
NewModel="UKBB.DEMScore.best.annot.20240904"
NewArgFile=sprintf("/home/rejudcu/pars/rsco.%s.rarg",NewModel)
New470KModel="UKBB.DEMScore.470K.annot.20240904"
New470KArgFile=sprintf("/home/rejudcu/pars/rsco.%s.rarg",New470KModel)
bashScript="run.DEMScore.PlusExWAS.best.annot.20241109.sh"
bashScript470K="run.DEMScore.PlusExWAS.470K.annot.20241109.sh"
ToDoGenesFile="DEMScore.PlusExWAS.270K.genes.txt"
ResultsFile="DEMScore.PlusExWAS.Results.txt"
ExWASGenes=c("SORL1", "GRN", "PSEN1", "ABCA7", "GBA", "ADAM10", "FRMD8", "DDX1", "DNMT3L", "MORC1", "TGM2")

DxCommand="dx cd / ; dx run swiss-army-knife -y  --ignore-reuse --instance-type mem3_hdd2_v2_x4 -imount_inputs=FALSE -iin=/scripts/runMultipleScoreassoc.20240822.sh "

wd="/home/rejudcu/UKBB/dementia.2024"
setwd(wd)

MLPThreshold=-log10(0.05/NumBest)

Summary=data.frame(read.table(SummFile,header=TRUE,stringsAsFactors=FALSE))
Summary$MaxMLP=apply(Summary[,2:ncol(Summary)], 1, max)# create new column
Summary=Summary[order(Summary$MaxMLP,decreasing=TRUE),]
Top=Summary[1:NumBest,]
ExtraGenes=ExWASGenes[! (ExWASGenes %in% Top[,1])]
Extra=Summary[Summary[,1] %in% ExtraGenes, ]
Top=rbind(Top,Extra)
TopGenes=Top[,1,drop=FALSE]
write.table(TopGenes,BestGenesFile,col.names=FALSE,row.names=FALSE,quote=FALSE)

AllDone=TRUE
LeftToDo=c("")
for (gene in Top[,1]) {
	SaoFile=sprintf("/cluster/project9/bipolargenomes/UKBB/%s/results/%s.%s.sao",BestModel,BestModel,gene)
	if (!file.exists(SaoFile)) {
		if (LeftToDo[1]=="") {
			LeftToDo=c(gene)
		} else {
			LeftToDo=c(LeftToDo,gene)
		}
		AllDone=FALSE
	}
}

if (!AllDone) {
	write(LeftToDo[1],"GenesLeftToDo.txt")
	if (length(LeftToDo)>1) {
		for (g in 2:length(LeftToDo)) {
			write(LeftToDo[g],"GenesLeftToDo.txt",append=TRUE)
		}
	}
	cmd=sprintf("bash /home/rejudcu/UKBB/RAPfiles/RAPscripts/getAllSaos.20240815.sh %s %s for.%s.lst",BestModel,"GenesLeftToDo.txt",BestModel)
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
GenesToDo=SummaryTable[SummaryTable$MLP270K>MLPThreshold,]

ExtraGenes=ExWASGenes[! (ExWASGenes %in% GenesToDo[,1])]
Extra=SummaryTable[SummaryTable[,1] %in% ExtraGenes, ]
GenesToDo=rbind(GenesToDo,Extra)

AllDone=TRUE
commStr=("dx cd /results/UKBB.DEMScore.470K.annot.20240904")
print(commStr)
system(commStr)
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
for (r in 1:nrow(GenesToDo)) {
	gene=GenesToDo[r,1]
	GenesToDo$DoThis[r]=FALSE
	SummFile=sprintf("/home/rejudcu/UKBB/dementia.2024/genes/summ.UKBB.DEMScore.470K.annot.20240904.%s.txt",gene)
	if (file.exists(SummFile)) {
		next
	}
	GenesToDo$DoThis[r]=TRUE
	cat(sprintf("Rscript scoreassoc.R --arg-file rsco.%s.rarg --gene %s --lintestfile %s --summaryoutputfile summ.UKBB.DEMScore.470K.annot.20240904.%s.txt\n",
		New470KModel,gene,GenesToDo[r,2],gene))
}
sink()
GenesToDo=GenesToDo[GenesToDo$DoThis==TRUE,1,drop=FALSE]
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
	if (SummaryTable$Gene[r] %in% GenesToDo$Gene) {
		gene=SummaryTable$Gene[r]
		SummFile=sprintf("/home/rejudcu/UKBB/dementia.2024/genes/summ.UKBB.DEMScore.470K.annot.20240904.%s.txt",gene)
		Summ=read.table(SummFile,header=TRUE,sep="\t",stringsAsFactors=FALSE)
		SummaryTable$MLP470K[r]=Summ[1,2]
	}
}
write.table(SummaryTable,ResultsFile,col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")






