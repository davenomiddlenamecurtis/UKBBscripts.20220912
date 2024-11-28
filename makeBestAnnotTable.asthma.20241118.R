#!/share/apps/R-3.6.1/bin/Rscript

# script to take summary file, select most significant gene/annotation results
# then run scoreassoc on second dataset with those pairs

NumBest=100
SummFile="/cluster/project9/bipolargenomes/UKBB/UKBB.asthma.annot.20240301/UKBB.asthma.annot.20240301.summ.sorted.txt"
BestModel="UKBB.asthma.forAnnot.20240610"
BestGenesFile="asthma.20241118.bestGenes.txt"
TestsFile="/home/rejudcu/pars/UKBB.annot.allTests.20231126.txt"
TopTestsFile="asthma.topTests.20241118.txt"
NewModel="UKBB.asthma.best.annot.20240610"
NewArgFile=sprintf("/home/rejudcu/pars/rsco.%s.rarg",NewModel)
New470KModel="UKBB.asthma.470K.annot.20241125"
New470KArgFile=sprintf("/home/rejudcu/pars/rsco.%s.rarg",New470KModel)
ResultsFile="Asthma.Results.20241118.txt"
ScoreAssocCommand="/home/rejudcu/scripts/subComm.sh /share/apps/R-3.6.1/bin/Rscript /home/rejudcu/scoreassoc/scoreassoc.R"

wd="/home/rejudcu/UKBB/asthma.2024"
setwd(wd)

MLPThreshold=-log10(0.05/NumBest)

Summary=data.frame(read.table(SummFile,header=TRUE,stringsAsFactors=FALSE))
Summary$MaxMLP=apply(Summary[,2:ncol(Summary)], 1, max)# create new column
Summary=Summary[order(Summary$MaxMLP,decreasing=TRUE),]
Top=Summary[1:NumBest,]
TopGenes=Top[,1]
write.table(TopGenes,BestGenesFile,col.names=FALSE,row.names=FALSE,quote=FALSE)

AllDone=TRUE
for (gene in Top[,1]) {
	ScoreFile=sprintf("/cluster/project9/bipolargenomes/UKBB/%s/results/%s.%s.sco",BestModel,BestModel,gene)
	if (!file.exists(ScoreFile)) {
		AllDone=FALSE
		break
	}
}

if (!AllDone) {
	cmd=sprintf("bash /home/rejudcu/UKBB/RAPfiles/scripts/getAllScores.20240708.sh %s %s for.%s.lst",BestModel,BestGenesFile,BestModel)
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

AllDone=TRUE
TopTests=data.frame(read.table(TopTestsFile,header=TRUE,sep="\t",stringsAsFactors=FALSE))
for (r in 1:nrow(TopTests)) {
	gene=TopTests[r,1]
	SummFile=sprintf("summ.%s.%s.txt",NewModel,gene)
	if (file.exists(SummFile)) {
		next
	}
	AllDone=FALSE
	commStr=sprintf("%s --arg-file %s --gene %s --summaryoutputfile %s --testfile %s",
		ScoreAssocCommand,NewArgFile,gene,SummFile,TopTests[r,2])
	print(commStr)
	system(commStr)
}

if (AllDone==FALSE) {
	quit()
}

# only get to here if all analyses have been run for 270K

SummaryTable=data.frame(read.table(TopTestsFile,header=TRUE,sep="\t",stringsAsFactors=FALSE))
SummaryTable$MLP270K=0
for (r in 1:nrow(TopTests)) {
	gene=TopTests[r,1]
	SummFile=sprintf("summ.%s.%s.txt",NewModel,gene)
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
	SummFile=sprintf("summ.%s.%s.txt",New470KModel,gene)
	if (!file.exists(SummFile)) {
		AllDone=FALSE
	}
}

if (!AllDone) {
for (r in 1:nrow(GenesToDo)) {
	gene=GenesToDo[r,1]
	GenesToDo$DoThis[r]=FALSE
	SummFile=sprintf("summ.%s.%s.txt",New470KModel,gene)
	if (file.exists(SummFile)) {
		next
	}
	GenesToDo$DoThis[r]=TRUE
	commStr=sprintf("%s --arg-file %s --gene %s --summaryoutputfile %s --testfile %s",
		ScoreAssocCommand,New470KArgFile,gene,SummFile,GenesToDo[r,2])
	print(commStr)
	system(commStr)
}
quit() # leave time for analyses to run
}

SummaryTable=data.frame(read.table(ResultsFile,header=TRUE,sep="\t",stringsAsFactors=FALSE))
SummaryTable$MLP470K=""
for (r in 1:nrow(SummaryTable)) {
	if (SummaryTable$MLP270K[r]>MLPThreshold) {
		gene=SummaryTable$Gene[r]
		SummFile=sprintf("summ.%s.%s.txt",New470KModel,gene)
		Summ=read.table(SummFile,header=TRUE,sep="\t",stringsAsFactors=FALSE)
		SummaryTable$MLP470K[r]=Summ[1,2]
	}
}
SummaryTable$TestFile=gsub("/home/rejudcu/pars/","",SummaryTable$TestFile)
write.table(SummaryTable,ResultsFile,col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")




