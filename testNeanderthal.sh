#!/bin/bash

tsvFile=covid19_result.20200428.txt
famFile=/SAN/ugi/UGIbiobank/data/downloaded/ukb51119_cal_chr21_v2_s488264.fam

positiveListFile=positiveList.txt

echo 3 45867532 45867532 rs35044562   > rs35044562.txt  
echo 3 45600000 46300000  LIMD1toCCR3   > NeandHap.txt  
plink  \
  --bed /SAN/ugi/UGIbiobank/data/downloaded/ukb_cal_chr3_v2.bed \
  --fam /SAN/ugi/UGIbiobank/data/downloaded/ukb51119_cal_chr21_v2_s488264.fam \
  --bim /SAN/ugi/UGIbiobank/data/downloaded/ukb_snp_chr3_v2.bim \
  --make-pheno $positiveListFile '*' \
  --extract 'range' NeandHap.txt \
  --all-pheno --allow-no-sex --freq case-control \
  --out  counts.NeandHap  
  
echo 3 45862952 45862952 rs71325088 > rs71325088.txt
plink  \
  --bed /SAN/ugi/UGIbiobank/data/downloaded/ukb_cal_chr3_v2.bed \
  --fam /SAN/ugi/UGIbiobank/data/downloaded/ukb51119_cal_chr21_v2_s488264.fam \
  --bim /SAN/ugi/UGIbiobank/data/downloaded/ukb_snp_chr3_v2.bim \
  --keep $positiveListFile \
  --extract 'range' rs71325088.txt \
  --recode A \
  --out  genos.positive.rs71325088  
  
plink  \
  --bed /SAN/ugi/UGIbiobank/data/downloaded/ukb_cal_chr3_v2.bed \
  --fam /SAN/ugi/UGIbiobank/data/downloaded/ukb51119_cal_chr21_v2_s488264.fam \
  --bim /SAN/ugi/UGIbiobank/data/downloaded/ukb_snp_chr3_v2.bim \
  --remove $positiveListFile \
  --extract 'range' rs71325088.txt \
  --recode A \
  --out  genos.negative.rs71325088  
  
  
