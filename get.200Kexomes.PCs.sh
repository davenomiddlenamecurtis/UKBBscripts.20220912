#!/bin/bash

# get 20 PCs from exomes using plink

allChrs="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X"

rm ukb23155.common.tomerge.txt
for c in $allChrs
do
  plink --maf 0.1 \
  --bed /home/rejudcu/vcf/UKBB.20201103/ukb23155_c${c}_b0_v1.bed \
  --fam /home/rejudcu/vcf/UKBB.20201103/ukb23155_c22_b0_v1_s200632.fam \
  --bim /home/rejudcu/vcf/UKBB.20201103/UKBexomeOQFE_chr${c}.bim \
  --make-bed \
  --out ukb23155.common.20201106.c$c
  if [ $c != 1 ]
  then
    echo ukb23155.common.20201106.c$c.bed ukb23155.common.20201106.c$c.bim ukb23155.common.20201106.c$c.fam >> ukb23155.common.tomerge.txt
  fi
done

plink --bfile ukb23155.common.20201106.c1 --merge-list ukb23155.common.tomerge.txt --make-bed --out ukb23155.common.all

/share/apps/genomics/plink-2.0/bin/plink2 --bfile ukb23155.common.all --pca 20 approx --out ukb23155.common.all

head -n 1 ukb23155.common.all.eigenvec | sed s/#// > ukb23155.common.all.eigenvec.txt
tail -n +2 ukb23155.common.all.eigenvec >> ukb23155.common.all.eigenvec.txt
