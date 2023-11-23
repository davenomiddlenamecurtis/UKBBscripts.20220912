#!/bin/bash
set -x
pwd
echo DX_WORKSPACE_ID is $DX_WORKSPACE_ID
echo DX_PROJECT_CONTEXT_ID is $DX_PROJECT_CONTEXT_ID
dx ls -l
unset DX_WORKSPACE_ID
dx cd $DX_PROJECT_CONTEXT_ID:
dx ls -l

dataFileSpec=ukb23158
mkdir ~/workdir
pushd ~/workdir

# get 20 PCs from exomes using plink
dx cd "/localExomes"
date
dx download $dataFileSpec.common.all.20230806.'*'
date
plink2 --bfile $dataFileSpec.common.all.20230806 --pca 20 approx --out $dataFileSpec.common.all.20230806
date
head -n 1 $dataFileSpec.common.all.20230806.eigenvec | sed s/#// > $dataFileSpec.common.all.20230806.txt
tail -n +2 $dataFileSpec.common.all.20230806.eigenvec >> $dataFileSpec.common.all.20230806.eigenvec.txt

popd
cp ~/workdir/$dataFileSpec.common.all.20230806.eigenvec.txt .

# invoke with:
# dx cd /covars ; dx run swiss-army-knife -y --instance-type mem3_hdd2_v2_x4  -imount_inputs=FALSE -iin="/scripts/get.500Kexomes.PCs.20230806.sh" -icmd="bash get.500Kexomes.PCs.20230806.sh "
# info re instance types is here: https://documentation.dnanexus.com/developer/api/running-analyses/instance-types


