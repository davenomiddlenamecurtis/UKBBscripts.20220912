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
genes=$3
neededFiles="/pars/gva.UKBB.all.20230807.arg /pars/LOF.weights.txt /phenos/UKBB.HT.txt" 

mkdir results # probably in folder ~/out/out
mkdir results/$testName
mkdir results/$testName/geneResults
mkdir results/$testName/failed

mkdir ~/workdir
pushd ~/workdir
dx cd pars
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
dx download $dataFileSpec${chr}'*'
dx download ukb23158_c22_b0_v1.fam
dx cd /annot
dx download ukb23158.annot.vcf.gz.tbi
# this is here to make sure it gets a later timestamp

for f in $neededFiles
do
	dx download $f
done

echo $genes | while read gene
do
		geneVarAssoc --arg-file $argFile --gene $gene --debug 1 --keep-temp-files 1
		if [ -e *.$gene.sao ]
		then
			pushd
			cp ~/workdir/*.$gene.s?o results/$testName/geneResults
			pushd
		else
			echo $gene > $gene.failed.txt
			ls -l *.$gene.* >> $gene.failed.txt
			pushd
			cp ~/workdir/$gene.failed.txt results/$testName/failed
			pushd
		fi
done
popd

# invoke with:
# genes=HERC1; chr=15; model=UKBB.LOF.20230809
# dx rm /results/$testName/geneResults/'*'.$gene.s'?'o ; dx rm /results/$testName/failed/$gene.failed.txt
# dx cd / ; dx run swiss-army-knife -y --instance-type mem3_hdd2_v2_x4  -imount_inputs=FALSE -iin="/scripts/runGeneLOF.20230809.sh"  -icmd="bash runGeneLOF.20230809.sh $model $chr $genes "
# info re instance types is here: https://documentation.dnanexus.com/developer/api/running-analyses/instance-types




