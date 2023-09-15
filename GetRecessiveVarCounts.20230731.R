#!/share/apps/R-3.6.1/bin/Rscript

TestName="HT.rec.20230728"

SLPFile="/cluster/project9/bipolargenomes/UKBB/UKBB.HT.rec.20230728/UKBB.HT.rec.20230728.summ.txt"
ResultsFolder="/cluster/project9/bipolargenomes/UKBB/UKBB.HT.rec.20230728/results"
SaoTemplate="/cluster/project9/bipolargenomes/UKBB/UKBB.HT.rec.20230728/results/UKBB.HT.rec.20230728.%s.sao"

VarCountsFile=sprintf("VarCounts.%s.txt",TestName)

SLPs=na.omit(data.frame(read.table(SLPFile,header=TRUE,stringsAsFactors=FALSE)))
print(sprintf("Number of genes with SLPs = %d",nrow(SLPs)))

SLPs=SLPs[SLPs$lrSLP!=0,]
print(sprintf("Number of genes with non-zero SLPs = %d",nrow(SLPs)))

VarCounts=data.frame(matrix(NA, ncol = 3, nrow = nrow(SLPs)),stringsAsFactors=FALSE)
colnames(VarCounts)=c("Gene","NVars","NAlts")
VarCounts$Gene=SLPs$Gene
for (r in 1:nrow(SLPs)) {
	SaoFile=sprintf(SaoTemplate,VarCounts$Gene[r])
	SaoLines=data.frame(read.table(SaoFile,header=FALSE,sep="",fill=TRUE,stringsAsFactors=FALSE))
	if (nrow(SaoLines)<3) {
		next
	}
	SaoLines=SaoLines[-c(1,2),]
	for (rr in 1:nrow(SaoLines)) {
		if (SaoLines[rr,10]=="") {
			break
		}
	}
	rr=rr-1
	SaoLines=SaoLines[1:rr,c(4,6,10,12)]
	SaoLines[]=lapply(SaoLines, as.numeric)
	VarCounts$NVars[r]=rr
	VarCounts$NAlts[r]=sum(SaoLines[,1])+sum(SaoLines[,3])+2*(sum(SaoLines[,2])+sum(SaoLines[,4]))
}
write.table(VarCounts,VarCountsFile,sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

