LDLimit=2000000

plink="/share/apps/genomics/plink-1.9/plink"
BedFileTemplate="/home/rejudcu/vcf/UKBB.20201103/ukb23155_c%s_b0_v1.bed"
BimFileTemplate="/home/rejudcu/vcf/UKBB.20201103/UKBexomeOQFE_chr%s.bim"
FamFile="/home/rejudcu/vcf/UKBB.20201103/ukb23155_c22_b0_v1_s200632.20210202.fam"
VariantsFile="variants.txt"
PCsFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb23155.common.all.eigenvec"
sexFile="/home/rejudcu/UKBB/UKBB.sex.20201111.txt"

Variants=data.frame(read.table(VariantsFile,header=TRUE,sep="",stringsAsFactors=FALSE))
LastChr=""
LastPos=-(LDLimit+1)
NumSets=0
Sets=data.frame(matrix(ncol=3,nrow=1))
colnames(Sets)=c("Chr","FirstPos","NumVars")
for (r in 1:nrow(Variants)) {
	# loop through, identify groups and save them separately
	if (Variants$chromosome_name[r]!=LastChr | Variants$pos[r]-LastPos>LDLimit) {
		ThisSet=Variants[r,]
		LastChr=Variants$chromosome_name[r]
		LastPos=Variants$pos[r]
		FirstSetPos=LastPos
		rr=2
	} else {
		LastPos=Variants$pos[r]
		ThisSet[rr,]=Variants[r,]
		rr=rr+1
	}
	SetDone=r==nrow(Variants)
	if (!SetDone) {
		SetDone=Variants$chromosome_name[r+1]!=LastChr | Variants$pos[r+1]-LastPos>LDLimit
	}
	if (SetDone & (LastChr!="6" | FirstSetPos !=26409662)) {
		NumSets=NumSets+1
		Sets[NumSets,1]=LastChr # Sets$Chr[NumSets]=LastChr does not work to implicitly add row to table
		Sets$FirstPos[NumSets]=FirstSetPos
		Sets$NumVars[NumSets]=nrow(ThisSet)
		SetFileName=sprintf("set.%s.%d.txt",LastChr,FirstSetPos)
		write.table(ThisSet,SetFileName,sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
	}
}
write.table(Sets,"AllSets.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

# now use plink to extract data for each set
Sets=data.frame(read.table("AllSets.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE))
for (r in 1:nrow(Sets)) {
	SetName=sprintf("%s.%d",Sets$Chr[r],Sets$FirstPos[r])
	if (!file.exists(sprintf("plink.%s.raw",SetName))) {
		SetFileName=sprintf("set.%s.txt",SetName)
		Set=data.frame(read.table(SetFileName,header=TRUE,sep="\t",stringsAsFactors=FALSE))
		Ranges=data.frame(matrix(ncol=4,nrow=nrow(Set)))
		Ranges[,1]=Set$chromosome_name
		Ranges[,2]=Set$pos
		Ranges[,3]=Set$pos
		Ranges[,4]=SetName
		RangeFile=sprintf("range.%s.txt",SetName)
		write.table(Ranges,RangeFile,sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
		BedFile=sprintf(BedFileTemplate,Sets$Chr[r])
		BimFile=sprintf(BimFileTemplate,Sets$Chr[r])
		CommLine=sprintf("%s --bed %s --bim %s --fam %s --extract range %s --set-hh-missing --recodeA --out plink.%s",plink,BedFile,BimFile,FamFile,RangeFile,SetName)
		system(CommLine)
	}
}

PCsTable=data.frame(read.table(PCsFile,header=FALSE,sep="\t"))
colnames(PCsTable)[1:2]=c("FID","IID")
covars=""
for (p in 1:20) {
  covars=sprintf("%sPC%d + ",covars,p)
  colnames(PCsTable)[2+p]=sprintf("PC%d",p) # because first row of PCs file starts with hash sign FID so row is ignored
}
sexTable=data.frame(read.table(sexFile,header=TRUE,sep="\t"))
covars=sprintf("%s Sex ",covars)

ToUse=vector()
PhenoNames=c("Asthma","AllergicRhinitis")
for (p in 1:2) {
	PhenoName=PhenoNames[p]
if (!file.exist(sprintf("AllMultiVar.%s.txt",PhenoName))) {
	FirstMulti=TRUE
	Pheno=data.frame(read.table(sprintf("UKBB.%s.txt",PhenoName),header=TRUE,sep="\t",stringsAsFactors=FALSE))
	colnames(Pheno)=c("IID","Pheno")
	PhenoData=merge(Pheno,sexTable,by="IID")
	PhenoData=merge(PhenoData,PCsTable,by="IID")
	FullModel=sprintf("Pheno ~ %s",covars)
	m=glm(as.formula(FullModel),data=PhenoData,family="binomial")
	LRTest.Coeffs0=summary(m)$coefficients
	LRTest.LL0=logLik(m)
    cat(sprintf("L0 = %f\n",LRTest.LL0))
	print(LRTest.Coeffs0)
	for (r in 1:nrow(Sets)) {
		SetName=sprintf("%s.%d",Sets$Chr[r],Sets$FirstPos[r])
		Genos=data.frame(read.table(sprintf("plink.%s.raw",SetName),header=TRUE,sep="",stringsAsFactors=FALSE))
		Genos[is.na(Genos)]=0
# NB missing genotypes are treated as wildtype - hopefully OK for relatively rare coding variants
		AllData=merge(PhenoData,Genos,by="IID")
		SingleVarResults=data.frame(matrix(ncol=5,nrow=1)) # was ncol(Genos)-6 but there may be variants with all zero genotypes
		colnames(SingleVarResults)=c("Var","Beta","SE","Z","PVal")
		NumToUse=0
		vv=0
		for (v in 7:ncol(Genos)) {
			FullModel=sprintf("Pheno ~ %s + %s",covars,colnames(Genos)[v])
			m=glm(as.formula(FullModel),data=AllData,family="binomial")
			LRTest.Coeffs1=summary(m)$coefficients
			LRTest.LL1=logLik(m)
			print(LRTest.Coeffs1)
			if (nrow(LRTest.Coeffs1)>22) { # can be a variant with all 0 genotypes
			vv=vv+1
			SingleVarResults[vv,1]=rownames(LRTest.Coeffs1)[23]
			SingleVarResults[vv,2:5]=LRTest.Coeffs1[23,]
			if (LRTest.Coeffs1[23,4]<0.05) {
				NumToUse=NumToUse+1
				ToUse[NumToUse]=rownames(LRTest.Coeffs1)[23]
			}
			}
		}
		write.table(SingleVarResults,sprintf("SingleVar.%s.%s.txt",PhenoName,SetName),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
		if (NumToUse>0) {
		FullModel=sprintf("Pheno ~ %s",covars)
		for (v in 1:NumToUse) {
			FullModel=sprintf("%s + %s",FullModel,ToUse[v])
		}
		m=glm(as.formula(FullModel),data=AllData,family="binomial")
		LRTest.Coeffs1=summary(m)$coefficients
		LRTest.LL1=logLik(m)
		print(LRTest.Coeffs1)
		MultiVarResults=data.frame(matrix(ncol=5,nrow=nrow(LRTest.Coeffs1)-22)) 
		# some variants can be dropped if singularities
		# though more likely the problem was due to missing genotypes
		colnames(MultiVarResults)=c("Var","Beta","SE","Z","PVal")
		MultiVarResults[,1]=rownames(LRTest.Coeffs1)[23:nrow(LRTest.Coeffs1)]
		MultiVarResults[,2:5]=LRTest.Coeffs1[23:nrow(LRTest.Coeffs1),]
		write.table(MultiVarResults,sprintf("MultiVar.%s.%s.txt",PhenoName,SetName),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
		}
		
		if (r==1) {
			AllSingleVarResults=SingleVarResults
		} else {
			AllSingleVarResults=rbind(AllSingleVarResults,SingleVarResults)
		}
		if (NumToUse>0) {
			if (FirstMulti) {
				AllMultiVarResults=MultiVarResults
				FirstMulti=FALSE
			} else {
			AllMultiVarResults=rbind(AllMultiVarResults,MultiVarResults)
			}
		}
		
	}
	write.table(AllSingleVarResults,sprintf("AllSingleVar.%s.txt",PhenoName),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
	write.table(AllMultiVarResults,sprintf("AllMultiVar.%s.txt",PhenoName),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
	}
}

Summary=Variants[Variants$chromosome_name!="6",]
Summary=Summary[c("chromosome_name","pos","alleleA","alleleB","Gene","Control.MAF","p.value","Odds.ratio","trait")]
Summary$Var=sprintf("X%s.%d.%s.%s",Summary$chromosome_name,Summary$pos,Summary$alleleA,Summary$alleleB)
NumOrigCols=ncol(Summary)
for (p in 1:2) {
	PhenoName=PhenoNames[p]
	AllSingleVarResults=data.frame(read.table(sprintf("AllSingleVar.%s.txt",PhenoName),header=TRUE,sep="\t",stringsAsFactors=FALSE))
	AllMultiVarResults=data.frame(read.table(sprintf("AllMultiVar.%s.txt",PhenoName),header=TRUE,sep="\t",stringsAsFactors=FALSE))
	AllSingleVarResults$Var=substr( AllSingleVarResults$Var,1,nchar(AllSingleVarResults$Var)-2)
	AllMultiVarResults$Var=substr( AllMultiVarResults$Var,1,nchar(AllMultiVarResults$Var)-2)
	if ("X1.152312600.D.4" %in% AllSingleVarResults$Var) {
		AllSingleVarResults$Var[AllSingleVarResults$Var=="X1.152312600.D.4"]="X1.152312600.CACTG.C"
	}
	if ("X1.152312600.D.4" %in% AllMultiVarResults$Var) {
		AllMultiVarResults$Var[AllMultiVarResults$Var=="X1.152312600.D.4"]="X1.152312600.CACTG.C"
	}
	# X17.39868373.I.3_CTCT is not in variants file, just at same position as an snv
	Summary=merge(Summary,AllSingleVarResults[c("Var","PVal")],by="Var",all.x=TRUE)
	colnames(Summary)[ncol(Summary)]=sprintf("SPV%s",PhenoName)
	Summary=merge(Summary,AllMultiVarResults[c("Var","PVal")],by="Var",all.x=TRUE)
	colnames(Summary)[ncol(Summary)]=sprintf("MPV%s",PhenoName)
}

write.table(Summary,sprintf("Summary.All.20230317.txt",PhenoName),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
