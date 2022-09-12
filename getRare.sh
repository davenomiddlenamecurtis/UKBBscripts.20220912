#!/bin/bash
# get all variants in specified range
MINMAC=100
MAXMAC=200

plink=/share/apps/genomics/plink-1.9/plink 

$plink \
 --bed /SAN/ugi/UGIbiobank/data/downloaded/ukb_efe_chr1_v1.bed \
 --fam /SAN/ugi/UGIbiobank/data/downloaded/ukb51119_efe_chr1_v1_s49953.fam \
 --bim /SAN/ugi/UGIbiobank/data/downloaded/ukb_fe_exm_chrall_v1.bim \
 --memory 2000000 \
 --mac $MINMAC --max-mac $MAXMAC \
 --chr 22 \
 --recode A --out rare.UKBB.22
 
 head rare.UKBB.22.raw | cut -f2-102 > UKBB.rare.22.first100.txt
 
 
 