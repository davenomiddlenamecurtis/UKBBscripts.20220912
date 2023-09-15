#!/share/apps/R-3.6.1/bin/Rscript

# R implementation of scoreassoc, using pre-calculated scores on new phenotypes
# library(plyr)

args = commandArgs(trailingOnly=TRUE)

wd="C:/Users/dave/OneDrive/msvc/data/fixTTest"
setwd(wd)
if (length(args)<1) {
  args=c("--arg-file","testADSP.rarg")
}

if (length(args)<1) {
  args=c("--arg-file","~/pars/rsco.UKBB.BMI.20210111.rarg",
     "--geneListFile", "/home/rejudcu/reference/DRDgenes.txt")
}

if (length(args)<1) {
args=c("--IDphenotypefile","/home/rejudcu/UKBB/lipids/UKBB.HL.20201103.txt",
"--dottest","1",
"--varfile","/SAN/ugi/UGIbiobank/data/downloaded/ukb23155.common.all.eigenvec.txt",
"--varfile","/home/rejudcu/UKBB/UKBB.sex.20201111.txt",
"--dolrtest","1","--testfile","/home/rejudcu/pars/justScore.tst",
"--geneListFile","/home/rejudcu/reference/DRDgenes.txt",
"--summaryoutputfile","UKBB.HL.withSex.R.SLPs.txt",
"--outputfilespec","UKBB.HL.withSex.R.GENE.rsao",
"--inputscorefilespec","/cluster/project9/bipolargenomes/UKBB/UKBB.HL.all.20201231/results/UKBB.HL.all.20201231.GENE.sco"
)
}


setClass("parInfo",slots=list(
IDphenotypefile="character",
inputScoreFileSpec="character",
gene="character",
geneListFile="character",
outputFileSpec="character",
summaryOutputFile="character",
doTTest="numeric",doLRTest="numeric",doLinRTest="numeric",isquantitative="numeric",
numVars="numeric",
numVarFiles="numeric",varFiles="vector",
numTestFiles="numeric",testFiles="vector",
numLinTestFiles="numeric",linTestFiles="vector"
))

pars=new("parInfo",numVarFiles=0,numTestFiles=0,numLinTestFiles=0,
doTTest=0,doLRTest=0,doLinRTest=0,isquantitative=0)

a=0
while (TRUE) {
  if (a*2+1>=length(args)) {
    break;
  }
  arg=args[a*2+1]
  if (arg=="--arg-file") {
    if ((a*2+3)<length(args)) {
	  oldArgs=args[(a*2+3):length(args)]
	} else {
	  oldArgs=c("")
	}
	newArgs=data.frame(read.table(args[a*2+2],header=FALSE,sep="",stringsAsFactors=FALSE))
	for (r in 1:nrow(newArgs)) {
	  args[(r-1)*2+1]=newArgs[r,1]
	  args[(r-1)*2+2]=newArgs[r,2]
	}
	args=append(args,oldArgs)
	a=-1
  } else if (arg=="--IDphenotypefile") {
    pars@IDphenotypefile=args[a*2+2]
  } else if (arg=="--isquantitative") {
    pars@isquantitative=as.numeric(args[a*2+2])
  } else if (arg=="--inputscorefilespec") {
    pars@inputScoreFileSpec=args[a*2+2]
  } else if (arg=="--geneListFile") {
    pars@geneListFile=args[a*2+2]
  } else if (arg=="--gene") {
    pars@gene=args[a*2+2]
  } else if (arg=="--outputfilespec") {
    pars@outputFileSpec=args[a*2+2]
  } else if (arg=="--summaryoutputfile") {
    pars@summaryOutputFile=args[a*2+2]
  } else if (arg=="--varfile") {
    pars@numVarFiles=pars@numVarFiles+1
    pars@varFiles[pars@numVarFiles]=args[a*2+2]
  } else if (arg=="--testfile") {
    pars@numTestFiles=pars@numTestFiles+1
    pars@testFiles[pars@numTestFiles]=args[a*2+2]
  } else if (arg=="--lintestfile") {
    pars@numLinTestFiles=pars@numLinTestFiles+1
    pars@linTestFiles[pars@numLinTestFiles]=args[a*2+2]
  } else if (arg=="--dottest") {
    pars@doTTest=as.numeric(args[a*2+2])
  } else if (arg=="--dolrtest") {
    pars@doLRTest=as.numeric(args[a*2+2])
  } else if (arg=="--dolinrtest") {
    pars@doLinRTest=as.numeric(args[a*2+2])
  }  else {
    print(sprintf("Error: Unrecognised command line argument: %s\n",arg))
	q()
  }
  a=a+1
}

phenoTypes=data.frame(read.table(pars@IDphenotypefile,header=FALSE,stringsAsFactors=FALSE,sep="",fill=TRUE))
if (phenoTypes[1,1]=="IID") {
  phenoTypes=phenoTypes[2:nrow(phenoTypes),]
}
colnames(phenoTypes)=c("IID","pheno")
phenoTypes=phenoTypes[complete.cases(phenoTypes),]
phenoTypes$pheno=as.numeric(phenoTypes$pheno)
if (!pars@isquantitative) {
  phenoTypes=phenoTypes[which(phenoTypes$pheno==0 | phenoTypes$pheno==1),]
}

if (pars@numVarFiles>0) {
first=TRUE
for (f in 1:pars@numVarFiles) {
  newVars=data.frame(read.table(pars@varFiles[f],header=TRUE))
  for (c in 1:ncol(newVars)) {
    if (colnames(newVars)[c]=="IID") {
	  break
	}
  }
  newVars=newVars[,c:ncol(newVars)]
  if (first) {
    vars=newVars
	first=FALSE
  } else {
    vars=merge(vars,newVars,by="IID")
  }
}
}

if (length(pars@geneListFile)>0) {
  genes=data.frame(read.table(pars@geneListFile,header=FALSE))[,1]
} else {
  genes=c(pars@gene)
}

if (file.exists(pars@summaryOutputFile)) {
  summary=data.frame(read.table(pars@summaryOutputFile,header=TRUE,stringsAsFactors=FALSE))
} else {
  summary=data.frame(matrix(nrow=0,ncol=1+pars@doTTest+pars@doLRTest+pars@doLinRTest+pars@numTestFiles+pars@numLinTestFiles))
  colnames(summary)[1]="Gene"
}
summaryRow=nrow(summary)+1

inputScoreFileSpec=sub("GENE","%s",pars@inputScoreFileSpec)
if (length(pars@outputFileSpec)>0) {
  outputFileSpec=sub("GENE","%s",pars@outputFileSpec)
}

if (pars@numVarFiles>0) {
allData=merge(phenoTypes,vars,by="IID")
if (pars@numVarFiles>0) {
    first=TRUE
    for (c in 2:ncol(vars)) {
	  if (first) {
        fullModel0=sprintf("pheno ~ %s",colnames(vars)[c])
	    first=FALSE
	  }
	  else {
	    fullModel0=sprintf("%s + %s",fullModel0,colnames(vars)[c])
	  }
	}
	fullModel1=sprintf("%s + score",fullModel0)
} else {
    fullModel0="pheno ~ 1"
    fullModel1="pheno ~ score"
}
} else {
  allData=phenoTypes
    fullModel0="pheno ~ 1"
    fullModel1="pheno ~ score"
}

#reduce allData rows to match subjects we are using
scores=data.frame(read.table(sprintf(inputScoreFileSpec,genes[1]),header=FALSE))
colnames(scores)=c("IID","oldPheno","score")
allData=merge(scores,allData,by="IID")
allData=allData[,c(1,4:ncol(allData))]

if (pars@doLRTest) {
  m=glm(as.formula(fullModel0),data=allData,family="binomial")
  LRTest.Coeffs0=summary(m)$coefficients
  LRTest.LL0=logLik(m)
}

if (pars@doLinRTest) {
  m=glm(as.formula(fullModel0),data=allData)
  LinRTest.Coeffs0=summary(m)$coefficients
  LinRTest.LL0=logLik(m)
}

tests.ndf0=tests.ndf1=tests.model0=tests.model1=tests.LL0=vector()
linTests.ndf0=linTests.ndf1=linTests.model0=linTests.model1=linTests.LL0=vector()

for (testTypes in 1:2) {
if (testTypes==1) {
nTests=pars@numTestFiles
} else {
nTests=pars@numLinTestFiles
}
for (t in 1:nTests) {
  if (nTests==0) {
    break
  }
  ndf0=1
  ndf1=1
  first0=TRUE
  first1=TRUE
  if (testTypes==1) {
    test=read.table(pars@testFiles[t],header=FALSE)
  } else{
    test=read.table(pars@linTestFiles[t],header=FALSE)
  }
  for (r in 1:nrow(test)) {
    if (test[r,3]==1) {
	  ndf0=ndf0+1
	  if (first0) {
	    model0=sprintf("pheno ~ %s",test[r,1])
	    first0=FALSE
	  } else {
	    model0=sprintf("%s + %s",model0,test[r,1])
	  }
	}
    if (test[r,4]==1) {
	  ndf1=ndf1+1
	  if (first1) {
	    model1=sprintf("pheno ~ %s",test[r,1])
	    first1=FALSE
	  } else {
	    model1=sprintf("%s + %s",model0,test[r,1])
	  }
	}
  }
  if (first0) {
    model0="pheno ~ 1"
  }
  if (first1) {
    model1="pheno ~ 1"
  }
  if (testTypes==1) {
    tests.ndf0[t]=ndf0
    tests.model0[t]=model0
    tests.ndf1[t]=ndf1
    tests.model1[t]=model1
  } else {
    linTests.ndf0[t]=ndf0
    linTests.model0[t]=model0
    linTests.ndf1[t]=ndf1
    linTests.model1[t]=model1
  }
}
}

for (gene in genes) {
  if (gene %in% summary$Gene) {
    next
  }
  if (length(pars@outputFileSpec)>0) {
    outFileName=sprintf(outputFileSpec,gene)
	sink(outFileName)
	cat(sprintf("Output for scoreassoc.R, %s\n",outFileName))
  }
  scoresFileName=sprintf(inputScoreFileSpec,gene)
  if (!file.exists(scoresFileName)) {
    cat(sprintf("Scores file %s does not exist\n",scoresFileName))
	sink()
	next
  }
  summary[summaryRow,1]=gene
  summaryCol=2
  scores=data.frame(read.table(scoresFileName,header=FALSE))
  colnames(scores)=c("IID","oldPheno","score")
  testData=merge(scores,allData,by="IID")
  if (pars@doTTest) {
    tt=t.test(testData$score[testData$pheno==0],testData$score[testData$pheno!=0])
	SLP=log10(tt$p.value)*as.numeric(sign(tt$estimate[1]-tt$estimate[2]))
	colnames(summary)[summaryCol]="SLP"
	summary[summaryRow,summaryCol]=SLP
	summaryCol=summaryCol+1
	print(t)
	cat(sprintf("SLP = %f (signed log10(p), positive if cases score higher than controls)\n\n",SLP))
  }
  if (pars@doLRTest) {
    cat(sprintf("L0 = %f\n",LRTest.LL0))
	print(LRTest.Coeffs0)
    m=glm(as.formula(fullModel1),data=testData,family="binomial")
	LL1=logLik(m)
    cat(sprintf("L1 = %f\n",LL1))
	print(summary(m)$coefficients)
	ch2=2*as.numeric(LL1-LRTest.LL0)
	p=pchisq(ch2,1,lower.tail=FALSE)
	SLP=-log10(p)*sign(summary(m)$coefficients[nrow(summary(m)$coefficients),1])
	cat(sprintf("chisq = %f, 1 df, p=%f\nlrSLP = %f (signed log10(p), positive if cases have more variants than controls)\n\n",
	  ch2,p,SLP)) 
	colnames(summary)[summaryCol]="lrSLP"
	summary[summaryRow,summaryCol]=SLP
	summaryCol=summaryCol+1
  }
  if (pars@doLinRTest) {
    cat(sprintf("L0 = %f\n",LinRTest.LL0))
	print(LinRTest.Coeffs0)
    m=glm(as.formula(fullModel1),data=testData)
	LL1=logLik(m)
    cat(sprintf("L1 = %f\n",LL1))
	print(summary(m)$coefficients)
	ch2=2*as.numeric(LL1-LinRTest.LL0)
	p=pchisq(ch2,1,lower.tail=FALSE)
	SLP=-log10(p)*sign(summary(m)$coefficients[nrow(summary(m)$coefficients),1])
	cat(sprintf("chisq = %f, 1 df, p=%f\nlinrSLP = %f (signed log10(p), positive if variant score is positively correlated with phenotype)\n\n",
	  ch2,p,SLP)) 
	colnames(summary)[summaryCol]="linrSLP"
	summary[summaryRow,summaryCol]=SLP
	summaryCol=summaryCol+1
  }
  
  for (testTypes in 1:2) {
    if (testTypes==1) {
      nTests=pars@numTestFiles
    } else {
      nTests=pars@numLinTestFiles
    }
    for (t in 1:nTests) {
      if (nTests==0) {
        break
      }
	  if (testTypes==1) {
	    m0=glm(as.formula(tests.model0[t]),data=testData,family="binomial")
	    m1=glm(as.formula(tests.model1[t]),data=testData,family="binomial")
		df=tests.ndf1[t]-tests.ndf0[t]
	  } else {
	    m0=glm(as.formula(tests.model0[t]),data=testData)
	    m1=glm(as.formula(tests.model1[t]),data=testData)
		df=linTests.ndf1[t]-linTests.ndf0[t]
	  }
	  LL0=logLik(m0)
	  LL1=logLik(m1)
	  coeffs0=summary(m0)$coefficients
	  coeffs1=summary(m1)$coefficients
	  ch2=2*as.numeric(LL1-LL0)
	  p=pchisq(ch2,df,lower.tail=FALSE)
	  SLP=-log10(p)
      cat(sprintf("L0 = %f\n",LL0))
	  print(coeffs0)
      cat(sprintf("L1 = %f\n",LL1))
	  print(coeffs1)

	  if (("score" %in% rownames(coeffs1)) && !("score" %in% rownames(coeffs0))) {
	    rr=match("score",rownames(coeffs1))
		if (coeffs1[rr,1]<0) {
		  SLP=-SLP
		}
     	cat(sprintf("chisq = %f, 1 df, p=%f\ntSLP = %f (signed log10(p), positive if variant score is positively correlated with phenotype)\n\n",
	      ch2,p,SLP)) 
	    colnames(summary)[summaryCol]="tSLP"
	  } else {
     	cat(sprintf("chisq = %f, 1 df, p=%f\ntMLP = %f (minus log10(p))\n\n",
	      ch2,p,SLP)) 
	    colnames(summary)[summaryCol]="tMLP"
	  }
	  summary[summaryRow,summaryCol]=SLP
	  summaryCol=summaryCol+1
	}
  }


  
  if (length(pars@outputFileSpec)>0) {
    sink()
  }
  write.table(summary,pars@summaryOutputFile,row.names=FALSE,quote=FALSE,sep="\t" )
  summaryRow=summaryRow+1
print(summary)
}

