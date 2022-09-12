#!/bin/bash

removeFile=subjectsToRemove.w51119.20200820.csv
dateStr=20200825

dataFolder=/cluster/project9/bipolargenomes/UKBB/downloads/

famFirst=ukb51119_cal_chr
famLast=_v2_s488264.fam
bedFirst=ukb_cal_chr
bedLast=_v2.bed
bimFirst=ukb_snp_chr
bimLast=_v2.bim
exFam=ukb51119_efe_chr1_v1_s49953.fam
exBed=ukb_efe_chr1_v1.bed
exBim=ukb_fe_exm_chrall_v1.bim
phenoFileRoots="ukb41465.enc_ukb ukb41465.exomes"

snpStr=ukb51119_cal_chr
exomeStr=ukb51119_efe_allchr_s49953


cd $dataFolder

for p in $phenoFileRoots; do grep -v -Ff $removeFile $p.txt > $p.$dateStr.txt; done

bedFile=$exBed
bimFile=$exBim
famFile=$exFam
outRoot=$exomeStr.$dateStr

plink --bed $bedFile --bim $bimFile --fam $famFile --remove-fam $removeFile --make-bed --out $outRoot

allChrs="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X"

for c in $allChrs
do
  bedFile=$bedFirst$c$bedLast
  bimFile=$bimFirst$c$bimLast
  famFile=$famFirst$c$famLast
  outRoot=$exomeStr.$dateStr.$c
  if [ ! -e $snpStr.$dateStr.$c.bed ]
  then
    plink --bed $bedFile --bim $bimFile --fam $famFile --remove-fam $removeFile --make-bed --out $outRoot
  fi
done
