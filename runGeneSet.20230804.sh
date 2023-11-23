#!/bin/bash 

# script to run on virtual machine in RAP using Swiss Army Knife app
# need to have run 
# dx cd /
# before submitting this
# then results will end up in subfolders of /results
set -x 

unset DX_WORKSPACE_ID
dx cd $DX_PROJECT_CONTEXT_ID:

dataFilePath="/Bulk/Exome sequences/Population level exome OQFE variants, PLINK format - final release/"
dataFileSpec=ukb23158_c

testName=$1
argFile=gva.$testName.arg
chr=$2
geneList=$3
# neededFiles="/pars/gva.UKBB.all.20230807.arg /pars/newWeights.20201231.txt /pars/POLYPHENWEIGHTS.20201231.txt /pars/SIFTWEIGHTS.20201231.txt /covars/ukb23158.common.all.20230806.eigenvec.txt /covars/UKBB.sex.20230807.txt /pars/justScore.tst" 
neededFiles="/pars/gva.UKBB.all.20230807.arg /pars/newWeights.20201231.txt /pars/POLYPHENWEIGHTS.20201231.txt /pars/SIFTWEIGHTS.20201231.txt /phenos/UKBB.HT.20230816.txt" 
# omit covariates for t test analysis
# not clear whether phenotype file should be hard wired here or should be on command line

mkdir results # probably in folder ~/out/out
mkdir results/$testName
mkdir results/$testName/geneResults
mkdir results/$testName/failed

mkdir ~/workdir
pushd ~/workdir
dx cd /geneLists
dx download $geneList 
dx cd /pars
dx download $argFile
dx cd /bin
dx download scoreassoc
dx download geneVarAssoc
chmod 755 scoreassoc geneVarAssoc
PATH=$PATH:.
dx cd /reference38
dx download refseqgenes.hg38.20191018.sorted.onePCDHG.txt
dx cd /annot
dx download ukb23158.annot.vcf.gz
dx cd "$dataFilePath"
dx download $dataFileSpec${chr}_'*'
dx download ukb23158_c22_b0_v1.fam
dx cd /annot
dx download ukb23158.annot.vcf.gz.tbi # here to make sure it is newer

for f in $neededFiles
do
	dx download $f
done

cat $geneList | while read gene
do
		geneVarAssoc --arg-file $argFile --gene $gene
		if [ -e *.$gene.sao ]
		then
			pushd
			cp ~/workdir/*.$gene.s?o results/$testName/geneResults
			dx rm /results/$testName/geneResults/'*'.$gene.s'?'o
			pushd
		else
			echo $gene > $gene.failed.txt
			ls -l *.$gene.* >> $gene.failed.txt
			pushd
			cp ~/workdir/$gene.failed.txt results/$testName/failed
			dx rm /results/$testName/failed/$gene.failed.txt
			pushd
		fi
done
popd

# invoke with:
# gene=DNMT3A; chr=2; model=UKBB.HT.20230807; geneList=$gene.txt; testName=$model
# echo $gene > $geneList; dx cd  /geneLists; dxupload $geneList
# dx rm /results/$testName/geneResults/'*'.$gene.s'?'o ; dx rm /results/$testName/failed/$gene.failed.txt
# dx cd / ; dx run swiss-army-knife -y --instance-type mem3_hdd2_v2_x4  -imount_inputs=FALSE -iin="/scripts/runGeneSet.20230804.sh" -icmd="bash runGeneSet.20230804.sh $model $chr $geneList "
# info re instance types is here: https://documentation.dnanexus.com/developer/api/running-analyses/instance-types




