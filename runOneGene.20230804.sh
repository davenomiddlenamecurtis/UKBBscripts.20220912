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
gene=$3
neededFiles=$4

dx mkdir /results/$testName
dx mkdir /results/$testName/geneResults

mkdir ~/workdir
pushd ~/workdir

dx cd /bin
dx download scoreassoc
dx download geneVarAssoc
chmod 755 scoreassoc geneVarAssoc
PATH=$PATH:.
dx cd /pars
dx download $argFile # does not handle nested arg files automatically but can be in otherNeededFiles
dx cd /reference38
dx download RefSeq.$chr.txt
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

geneVarAssoc --arg-file $argFile --gene $gene
dx cd /results/$testName/geneResults
for f in *.$gene.s?o
do
	dx rm $f
	dx upload $f
done

popd
cp ~/workdir/*.sao .
cp ~/workdir/*.sco .



