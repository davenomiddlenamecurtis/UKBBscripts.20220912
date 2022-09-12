#!/bin/bash

cut ukb23155_c22_b0_v1_s200632.fam -f 1 -d " " | grep -v "-" | sort > UKBB.exome.IDs.20201103.txt 

head -n 1 ukb41465.txt >ukb41465.exomes.20201103.txt 

tail -n +2  ukb41465.txt | cut -f 2- | join  -t $'\t' - UKBB.exome.IDs.20201103.txt >> ukb41465.exomes.20201103.txt 

head -n 1 ukb43357.txt >ukb43357.exomes.20201103.txt 

tail -n +2  ukb43357.txt | cut -f 2- | join  -t $'\t' - UKBB.exome.IDs.20201103.txt >> ukb43357.exomes.20201103.txt 

# the original files have an unwanted tab before the eid field but not in the first line

# head -n 1 ukb51119.40PCs.txt > UKBB.exome.40PCs.20201103.txt
# grep -f UKBB.exome.IDs.20201103.txt ukb51119.40PCs.txt >> UKBB.exome.40PCs.20201103.txt
# cut UKBB.exome.40PCs.20201103.txt -f 1-21 -d " " > UKBB.exome.20PCs.20201103.txt