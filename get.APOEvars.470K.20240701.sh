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

genoName=APOE
chr=19


mkdir results # probably in folder ~/out/out
mkdir results/$genoName


mkdir ~/workdir
pushd ~/workdir


dx cd "$dataFilePath"
dx download $dataFileSpec${chr}_'*'

SNPFile=APOEVars.txt
echo "19:44908684:T:C" > $SNPFile
echo "19:44908822:C:T" >> $SNPFile

OutputRoot=APOE470Kgenotypes


plink --bfile ${dataFileSpec}${chr}_b0_v1 --extract $SNPFile --recode A --out $OutputRoot


popd
cd results/$genoName
cp ~/workdir/$OutputRoot.* .



# invoke with:
# dx cd / ; dx run swiss-army-knife -y --instance-type mem3_hdd2_v2_x4  -imount_inputs=FALSE -iin="/scripts/get.APOEvars.470K.20240701.sh" -icmd="bash get.APOEvars.470K.20240701.sh "
# info re instance types is here: https://documentation.dnanexus.com/developer/api/running-analyses/instance-types


