#!/bin/bash
unset DX_WORKSPACE_ID
dx cd $DX_PROJECT_CONTEXT_ID:

chr=$1

echo I have set chr to be $chr

mkdir ~/workdir
pushd ~/workdir

dx ls
dx cd "/Bulk/Exome sequences/Population level exome OQFE variants, PLINK format"
dx download "ukb23155_c${chr}_b0_v1.bed"
dx download "ukb23155_c${chr}_b0_v1.bim"
dx download "ukb23155_c${chr}_b0_v1.fam"

sed -i 's/redacted/-9/g' ukb23155_c${chr}_b0_v1.fam

dx cd
dx cd localExomes
# above lines do not do anything useful

plink2 -bfile ukb23155_c${chr}_b0_v1  --maf 0.1 --make-bed --out ukb23155.common.c${chr}.20210829

# dx upload "ukb23155.common.c${chr}.20210829.bed"
# dx upload "ukb23155.common.c${chr}.20210829.bim"
# dx upload "ukb23155.common.c${chr}.20210829.fam"

# rm -f *

popd
cp ~/workdir/ukb23155.common.c${chr}.20210829.bed .
cp ~/workdir/ukb23155.common.c${chr}.20210829.bim .
cp ~/workdir/ukb23155.common.c${chr}.20210829.fam .


# invoke with:
# dx cd /localExomes ; dx run swiss-army-knife -y --instance-type mem3_hdd2_v2_x4  -imount_inputs=FALSE -iin="/scripts/getCommonVars.20210829.sh" -icmd="bash getCommonVars.20210829.sh $CHR"
# largest file is 83 G
# info re instance types is here: https://documentation.dnanexus.com/developer/api/running-analyses/instance-types

# or:
# for ( $c=1 ; $c -lt 22 ; $c++) { dx cd /localExomes ; dx run swiss-army-knife -y --instance-type mem3_hdd2_v2_x4  -imount_inputs=FALSE -iin="/scripts/getCommonVars.20210829.sh" -icmd="bash getCommonVars.20210829.sh $c" ; }
