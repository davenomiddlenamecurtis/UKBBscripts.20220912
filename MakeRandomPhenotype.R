#!/share/apps/R-3.6.1/bin/Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args)<1) {
	args=c("UKBB.HT.txt","22")
}
PhenoFile=args[1]
Phenos=data.frame(read.table(PhenoFile,header=FALSE,stringsAsFactors=FALSE,fill=TRUE)) # will work whether or not there is a header
Phenos=Phenos[Phenos[,2]==0 | Phenos[,2]==1 ,]
Phenos[,2]=sample(as.numeric(Phenos[,2]),nrow(Phenos))
RandomFile=sprintf("%s.random.txt",args[1])
write.table(Phenos,RandomFile,col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
