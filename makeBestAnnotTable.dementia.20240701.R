#!/share/apps/R-3.6.1/bin/Rscript

# script to take summary file, select most significant gene/annotation results
# then run scoreassoc on second dataset with those pairs

NumBest=100
SummFile="/cluster/project9/bipolargenomes/UKBB/UKBB.dementia.annot.20240620/UKBB.dementia.annot.20240620.summ.txt"
BestModel="UKBB.dementia.forAnnot.20240620"
BestGenesFile="dementia.20240620.bestGenes.txt"
TestsFile="/home/lgibbons/pars/UKBB.annot.allTests.withAPOE.20240620.txt"
TopTestsFile="dementia.topTests.20240620.txt"
NewModel="UKBB.dementia.best.annot.20240620"
NewArgFile=sprintf("/home/lgibbons/pars/rsco.%s.rarg",NewModel)

ScoreAssocCommand="bash ~/scripts/subComm.sh /share/apps/R-3.6.1/bin/Rscript /home/rejudcu/scoreassoc/scoreassoc.R"


wd="/home/lgibbons/dementia.20240701"
setwd(wd)

Summary=data.frame(read.table(SummFile,header=TRUE,sep="\t",stringsAsFactors=FALSE))
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
	cmd=sprintf("bash /home/rejudcu/UKBB/RAPfiles/scripts/getAllScores.20240118.sh %s %s for.%s.lst",BestModel,BestGenesFile,BestModel)
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
for (r in 1:nrow(TopTests)) {
	gene=TopTests[r,1]
	commStr=sprintf("%s --arg-file %s --gene %s --summaryoutputfile summ.%s.%s.txt --testfile %s",
		ScoreAssocCommand,NewArgFile,gene,NewModel,gene,TopTests[r,2])
	print(commStr)
	system(commStr)
}

