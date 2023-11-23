#!/bin/bash 
# runOneGene.20230804.sh arg-file chromosome gene otherNeededFiles
# script to run on virtual machine in RAP using Swiss Army Knife app

unset DX_WORKSPACE_ID
dx cd $DX_PROJECT_CONTEXT_ID:

dataFilePath="/Bulk/Exome sequences/Population level exome OQFE variants, PLINK format"
dataFileSpec=ukb23155_c
# these will change

testName=$1
argFile=gva.$testName.arg
chr=$2
neededFiles=$3

summFile=$testName.summ.txt


getSLPs='
BEGIN { ORS=""; nSLP=0; } 
{
if ($1 == "SLP" || $1 == "tSLP" || $1 == "tMLP" || $1 == "linrSLP" || $1 == "lrSLP") 
	{
	nSLP=nSLP+1;
	SLPs[nSLP]=$3;
	}
}
END {
    if (nSLP > 0) {
	print gene "\t";
	for (i=1; i<=nSLP; ++i) {
	print SLPs[i] "\t";
	}
	print "\n";
	}
}
'

dx mkdir /results/$testName
dx mkdir /results/$testName/geneResults

mkdir ~/workdir
pushd ~/workdir
dx cd /results/$testName
dx download $summFile
dx cd /bin
dx download scoreassoc
dx download geneVarAssoc
chmod 755 scoreassoc geneVarAssoc
PATH=$PATH:.
dx cd /pars
dx download $argFile # does not handle nested arg files automatically but can be in otherNeededFiles
dx cd /reference38
dx download RefSeq.$chr.txt
dx download AllGenes.$chr.txt
dx cd /annot
dx download annot.$chr.vcf.gz
dx download annot.$chr.vcf.gz.tbi
cd $dataFilePath
dx download $dataFileSpec${chr}'*'

if [ $neededFiles 1= . ] # will probably need at least PCs, sex, phenotype file
then
	cat $neededFiles | while read d f
	do
		dx cd $d
		dx download $f
	done
fi
cat AllGenes.$chr.txt | while read gene
do
	g=
	if [ -z "`grep -w $gene $summFile`" != . ]
	then
		geneVarAssoc --arg-file $argFile --gene $gene
		if [ -e *.$gene.sao ]
		then
			dx cd /results
			for f in *.$gene.s?o
			do
				dx rm $f
				dx upload $f
			done
			awk -v gene=$gene "$getSLPs" *.$gene.sao >> $summFile
			dx cd /results/$testName
			dx rm $summFile
			dx upload $summFile
		else
			echo $gene > $gene.failed.txt
			dx mkdir /results/$testName/failed
			dx cd /results/$testName/failed
			dx upload $gene.failed.txt
		fi
	fi

popd
cp ~/workdir/*.sao .
cp ~/workdir/*.sco .



