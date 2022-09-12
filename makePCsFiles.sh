#!/bin/bash

First20PCs="PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 PC11 PC12 PC13 PC14 PC15 PC16 PC17 PC18 PC19 PC20"
Second20PCs=" PC21 PC22 PC23 PC24 PC25 PC26 PC27 PC28 PC29 PC30 PC31 PC32 PC33 PC34 PC35 PC36 PC37 PC38 PC39 PC40"

File20=ukb51119.20PCs.txt
File40=ukb51119.40PCs.txt
# echo IID FID $First20PCs>$File20
echo IID $First20PCs$Second20PCs>$File40

cut ukb51119_cal_chr9_v2_s488264.fam -f 1 -d " " > ukb51119.IDs.txt

cut ukb_sqc_v2.txt -f 26-65 -d " " | paste ukb51119.IDs.txt - -d " " >> $File40
cut $File40 -f 1-21 -d " " >$File20

