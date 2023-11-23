#!/bin/bash
set -x
unset DX_WORKSPACE_ID
dx cd $DX_PROJECT_CONTEXT_ID:

# dataFileSpec=/mnt/project/localExomes/ukb23158

dataFileSpec=ukb23158
outputFileSpec=ukb23158

mkdir ~/workdir
pushd ~/workdir

dx cd "/localExomes"
dx download $dataFileSpec.'*'.20230806.'*'
ls -l

allChrs="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"

for c in $allChrs
do
  if [ $c != 1 ]
  then
    echo $dataFileSpec.$c.20230806.bed $dataFileSpec.$c.20230806.bim $dataFileSpec.$c.20230806.fam >> $outputFileSpec.common.tomerge.txt
  fi
done

plink --bfile $dataFileSpec.1.20230806 --merge-list $dataFileSpec.common.tomerge.txt --make-bed --out $outputFileSpec.common.all.20230806
ls -lrt
pwd

popd

cp ~/workdir/$outputFileSpec.common.all.20230806.bed .
cp ~/workdir/$outputFileSpec.common.all.20230806.bim .
cp ~/workdir/$outputFileSpec.common.all.20230806.fam .
# cp ~/workdir/$outputFileSpec.common.all* .


# invoke with:
# dx cd /localExomes ; dx run swiss-army-knife -y --instance-type mem3_hdd2_v2_x4  -imount_inputs=FALSE -iin="/scripts/mergeCommonVars.20230806.sh" -icmd="bash mergeCommonVars.20230806.sh "
# info re instance types is here: https://documentation.dnanexus.com/developer/api/running-analyses/instance-types

