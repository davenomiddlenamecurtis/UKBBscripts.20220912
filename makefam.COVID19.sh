#!/bin/bash

tsvFile=covid19_result.20200428.txt
famFile=/cluster/project9/bipolargenomes/UKBB/downloads/ukb51119_cal_chr21_v2_s488264.fam

# in the first instance, just extract the subjects who scored positive to get their MAF]
positiveListFile=positiveList.txt

tail -n +2 $tsvFile | sort | cut -f1,6 | grep -w 1 | while read a b ; do echo $a $a ;done > tempFile.txt
cat tempFile.txt | sort |uniq > $positiveListFile
echo rs12329760 > rs12329760.txt
plink  \
  --bed /cluster/project9/bipolargenomes/UKBB/downloads/ukb_cal_chr21_v2.bed \
  --fam /cluster/project9/bipolargenomes/UKBB/downloads/ukb51119_cal_chr21_v2_s488264.fam \
  --bim /cluster/project9/bipolargenomes/UKBB/downloads/ukb_snp_chr21_v2.bim \
  --keep $positiveListFile \
  --extract rs12329760.txt \
  --recode \
  --out  positive.rs12329760
  
plink  \
  --bed /cluster/project9/bipolargenomes/UKBB/downloads/ukb_cal_chr21_v2.bed \
  --fam /cluster/project9/bipolargenomes/UKBB/downloads/ukb51119_cal_chr21_v2_s488264.fam \
  --bim /cluster/project9/bipolargenomes/UKBB/downloads/ukb_snp_chr21_v2.bim \
  --make-pheno $positiveListFile '*' \
  --extract rs12329760.txt \
  --all-pheno --allow-no-sex --assoc counts \
  --out  assoc.rs12329760

echo 21 41461493 41461693 rs35074065  > rs35074065.txt  
plink  \
  --bed /cluster/project9/bipolargenomes/UKBB/downloads/ukb_cal_chr21_v2.bed \
  --fam /cluster/project9/bipolargenomes/UKBB/downloads/ukb51119_cal_chr21_v2_s488264.fam \
  --bim /cluster/project9/bipolargenomes/UKBB/downloads/ukb_snp_chr21_v2.bim \
  --make-pheno $positiveListFile '*' \
  --extract 'range' rs35074065.txt \
  --all-pheno --allow-no-sex --freq case-control \
  --out  counts.rs35074065 
  
