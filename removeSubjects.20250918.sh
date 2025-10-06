#!/bin/bash

removeFile=doNotUse.20250918.txt
dateStr=20250918

dataFolder=/SAN/ugi/UGIbiobank/data/downloaded

famFirst=ukb51119_cal_chr
famLast=_v2_s488264.fam
bedFirst=ukb_cal_chr
bedLast=_v2.bed
bimFirst=ukb_snp_chr
bimLast=_v2.bim
exFam=ukb51119_efe_chr1_v1_s49953.fam
exBed=ukb_efe_chr1_v1.bed
exBim=ukb_fe_exm_chrall_v1.bim
phenoFileRoots="ukb41465.enc_ukb ukb43357.enc_ukb ukb41465.470Kexomes ukb43357.470Kexomes"

bgenFirst=ukb22438_c
bgenLast=_b0_v2
sampleLast=_b0_v2_s487038

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
  outRoot=$snpStr.$dateStr.$c
  if [ ! -e $snpStr.$dateStr.$c.bed ]
  then
    plink --bed $bedFile --bim $bimFile --fam $famFile --remove-fam $removeFile --make-bed --out $outRoot
  fi
  bgenFile=$bgenFirst$c$bgenLast.bgen
  sampleFile=$bgenFirst$c$sampleLast.sample
  outRoot=$bgenFirst$c$sampleLast.$dateStr
  if [ ! -e $outRoot.bgen ]
  then
    /share/apps/genomics/plink-2.0/bin/plink2 --bgen $bgenFile --sample $sampleFile --remove-fam $removeFile --export bgen-1.1 --out $outRoot
  fi
done
