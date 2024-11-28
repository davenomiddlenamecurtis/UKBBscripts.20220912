#!/share/apps/R-3.6.1/bin/Rscript

# script get correlations between variant scores and principal components, etc

ScoreFile="UKBB.DEMScore.countsOnly.20241108.ABCA7.sco"
PCsFile="ukb23158.common.all.20230806.eigenvec.txt"
SexFile="UKBB.sex.20230807.txt"
APOEFile="APOEStatus.470K.txt"
OutFile="Corr.txt

MainWeights=c("VEPWeight","FivePrime","InDelEtc","IntronicEtc","LOF","ProteinAltering","SpliceRegion","Synonymous","ThreePrime")
ExtraWeightsFile="~/UKBB/RAPfiles/annot/extraWeights.20231106.txt"
ExtraWeights=data.frame(read.table(ExtraWeightsFile,header=TRUE,stringsAsFactors=FALSE,fill=TRUE))[,1]

Scores=data.frame(read.table(ScoresFile,header=FALSE,stringsAsFactors=FALSE)))
colnames(Scores)=c("IID","Pheno",MainWeights,"AM_prediction","AM_Score",ExtraWeights)
PCs=data.frame(read.table(PCsFile,header=TRUE,stringsAsFactors=FALSE)))
Sex=data.frame(read.table(SexFile,header=TRUE,stringsAsFactors=FALSE)))
APOE=data.frame(read.table(APOEFile,header=TRUE,stringsAsFactors=FALSE)))
AllVars=merge(Scores,PCs,by="IID")
AllVars=merge(AllVars,Sex,by="IID")
AllVars=merge(AllVars,APOE,by="IID")

M=cor(AllVars)
corrMatrix=data.frame(matrix(ncol=ncol(AllVars)+1,nrow=length(ncol(AllVars))))
colnames(corrMatrix)=c("Variable",colnames(AllVars))
corrMatrix[,1]=colnames(AllVars)
corrMatrix[,2:(ncol(AllVars)+1)]=M
write.table(corrMatrix,sprintf("corrMatrix.%s.txt",ScoreFile),row.names=FALSE,quote=FALSE,sep="\t")
