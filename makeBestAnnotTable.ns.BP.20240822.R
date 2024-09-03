#!/share/apps/R-3.6.1/bin/Rscript

# script to take summary file, select most significant gene/annotation results
# then run scoreassoc on second dataset with those pairs

NumBest=100
SummFile="/cluster/project9/bipolargenomes/UKBB/UKBB.BP.annot.ns.20240821/UKBB.BP.annot.ns.20240821.summ.txt"
BestModel="UKBB.BP.forAnnot.20240814"
BestGenesFile="BP.20240822.bestGenes.txt"
TestsFile="/home/rejudcu/UKBB/RAP/pars/UKBB.annot.ns.allTests.20240822.txt"
TopTestsFile="BP.topTests.20240822.txt"
NewModel="UKBB.BP.best.annot.ns.20240822"
bashScript="run.BP.best.annot.ns.20240822.sh"
NewArgFile=sprintf("/home/rejudcu/pars/rsco.%s.rarg",NewModel)

DxCommand="dx cd / ; dx run swiss-army-knife -y  --ignore-reuse --instance-type mem3_hdd2_v2_x4 -imount_inputs=FALSE -iin=/scripts/runMultipleScoreassoc.20240822.sh "


wd="/home/rejudcu/UKBB/BP.20240814"
setwd(wd)

Summary=data.frame(read.table(SummFile,header=TRUE,stringsAsFactors=FALSE))
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
	cmd=sprintf("bash /home/rejudcu/UKBB/RAP/scripts/getAllSaos.20240815.sh %s %s for.%s.lst",BestModel,BestGenesFile,BestModel)
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
sink(bashScript)
for (r in 1:nrow(TopTests)) {
	gene=TopTests[r,1]
	cat(sprintf("Rscript scoreassoc.R --arg-file rsco.%s.rarg --gene %s --testfile %s\n",
		NewModel,gene,TopTests[r,2]))
}
sink()

commStr=sprintf("dx cd /geneLists; dxupload %s",BestGenesFile)
print(commStr)
system(commStr)
commStr=sprintf("dx cd /scripts; dxupload %s",bashScript)
print(commStr)
system(commStr)
cmd=sprintf("bash runMultipleScoreassoc.20240822.sh %s for.%s.lst %s",
		NewModel,NewModel,bashScript)
commStr=sprintf("%s -icmd=\"%s\"",DxCommand,cmd)
print(commStr)
system(commStr)


