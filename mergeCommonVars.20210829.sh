#!/bin/bash
unset DX_WORKSPACE_ID
dx cd $DX_PROJECT_CONTEXT_ID:

mkdir ~/workdir
pushd ~/workdir

dx cd "/localExomes"
dx download "ukb23155.common.c*.20210829.*"


allChrs="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X"

rm ukb23155.common.tomerge.txt
for c in $allChrs
do
  if [ $c != 1 ]
  then
    echo ukb23155.common.c$c.20210829.bed ukb23155.common.c$c.20210829.bim ukb23155.common.c$c.20210829.fam >> ukb23155.common.tomerge.txt
  fi
done

plink --bfile ukb23155.common.c1.20210829 --merge-list ukb23155.common.tomerge.txt --make-bed --out ukb23155.common.all.20210829

popd
cp ~/workdir/ukb23155.common.all.20210829.bed .
cp ~/workdir/ukb23155.common.all.20210829.bim .
cp ~/workdir/ukb23155.common.all.20210829.fam .


# invoke with:
# dx cd /localExomes ; dx run swiss-army-knife -y --instance-type mem3_hdd2_v2_x4  -imount_inputs=FALSE -iin="/scripts/mergeCommonVars.20210829.sh" -icmd="bash mergeCommonVars.20210829.sh "
# info re instance types is here: https://documentation.dnanexus.com/developer/api/running-analyses/instance-types

