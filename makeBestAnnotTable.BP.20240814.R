#!/share/apps/R-3.6.1/bin/Rscript

# script to take summary file, select most significant gene/annotation results
# then run scoreassoc on second dataset with those pairs

NumBest=3
SummFile="/home/rejudcu/UKBB/BP.20240814/genes/UKBB.BP.annot.20240814.summ.txt"
BestModel="UKBB.BP.forAnnot.20240814"
BestGenesFile="BP.20240814.bestGenes.txt"
TestsFile="/home/rejudcu/UKBB/RAPfiles/pars/UKBB.annot.allTests.20240815.txt"
TopTestsFile="BP.topTests.20240814.txt"
NewModel="UKBB.BP.best.annot.20240815"
NewArgFile=sprintf("/home/rejudcu/pars/rsco.%s.rarg",NewModel)

DxCommand="dx cd / ; dx run swiss-army-knife -y --instance-type mem3_hdd2_v2_x4 -imount_inputs=FALSE -iin=/scripts/runScoreassoc.20240815.sh "


wd="/home/rejudcu/UKBB/BP.20240814"
setwd(wd)

Summary=data.frame(read.table(SummFile,header=TRUE,sep="\t",stringsAsFactors=FALSE))
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
for (r in 1:nrow(TopTests)) {
	gene=TopTests[r,1]
	icmd=sprintf("bash runScoreassoc.20240815.sh %s for.%s.lst --gene %s --summaryoutputfile %s.%s.summ.txt --testfile %s",
		NewModel,NewModel,gene,NewModel,gene,TopTests[r,2])
	commStr=sprintf("%s -icmd=\"%s\"",DxCommand,icmd)
	print(commStr)
	system(commStr)
}

