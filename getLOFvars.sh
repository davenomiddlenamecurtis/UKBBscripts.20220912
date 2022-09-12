#!/bin/bash

# get list of LOF variants which have passed QC

genes="SETD1A CUL1 XPO7 TRIO CACNA1G SP4 GRIA3 GRIN2A  HERC1 RB1CC1"

for g in $genes
do
  echo $g >> LOFs.txt
  cut -f 4 < IDsAndCounts.good.LOF.20210615.$g.txt | sort | uniq -c >> LOFs.txt
  echo >> LOFs.txt
done