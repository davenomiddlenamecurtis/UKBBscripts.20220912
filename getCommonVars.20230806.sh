#!/bin/bash

# this is all wrong
# files are in /mnt/project but it is very slow to access them this way
# use dx download instead


unset DX_WORKSPACE_ID
dx cd $DX_PROJECT_CONTEXT_ID:

chr=$1

echo I have set chr to be $chr

dataFilePath="/Bulk/Exome sequences/Population level exome OQFE variants, PLINK format - final release/"
dataFileSpec=ukb23158

mkdir ~/workdir
pushd ~/workdir

dx ls
# use dx download instead of following lines
# cp /mnt/project/"$dataFilePath"/${dataFileSpec}_c${chr}_b0_v1.fam .
# ln -s /mnt/project/"$dataFilePath"/${dataFileSpec}_c${chr}_b0_v1.b* .

sed -i 's/redacted/-9/g' ${dataFileSpec}_c${chr}_b0_v1.fam

plink2 -bfile ${dataFileSpec}_c${chr}_b0_v1  --maf 0.1 --make-bed --out $dataFileSpec.${chr}.20230806

popd
dx cd /localExomes
cp ~/workdir/$dataFileSpec.${chr}.20230806.bed .
cp ~/workdir/$dataFileSpec.${chr}.20230806.bim .
cp ~/workdir/$dataFileSpec.${chr}.20230806.fam .


# invoke with:
# dx cd /localExomes ; dx run swiss-army-knife -y --instance-type mem3_hdd2_v2_x4  -imount_inputs=FALSE -iin="/scripts/getCommonVars.20230806.sh" -icmd="bash getCommonVars.20230806.sh $CHR"
# largest file is 83 G
# info re instance types is here: https://documentation.dnanexus.com/developer/api/running-analyses/instance-types

# or:
# allChrs="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"
# for c in $allChrs; do dx cd /localExomes ; dx run swiss-army-knife -y --instance-type mem3_hdd2_v2_x4  -imount_inputs=FALSE -iin="/scripts/getCommonVars.20230806.sh" -icmd="bash /mnt/project/scripts/getCommonVars.20230806.sh $c" ; done
# for c in $allChrs; do dx cd /localExomes ; dx run swiss-army-knife -y --instance-type mem3_hdd2_v2_x4 -icmd="bash /mnt/project/scripts/getCommonVars.20230806.sh $c" ; done
