#!/bin/bash

root=ukb23158_500k_OQFE.annotations

allChrs="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23"

for c in $allChrs
do
if [ ! -e $root.$c.annot.vcf ]
then
  head -n 5 $root$c.annot.vcf > $root.$c.annot.vcf
  tail -n +6 $root$c.annot.vcf | sort -k 2n >> $root.$c.annot.vcf
fi
done

if [ ! -e $root.all.annot.vcf ]
then
head -n 5 $root1.annot.vcf >$root.all.annot.vcf
for c in $allChrs
do
  tail -n +6 $root.$c.annot.vcf >> $root.all.annot.vcf
done
fi

for c in $allChrs all
do 
  subComm.sh "bgzip $root.$c.annot.vcf; tabix -p vcf $root.$c.annot.vcf.gz"
done