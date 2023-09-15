#!/bin/bash

head -n 1 ukb41465.txt > ukb41465.470Kexomes.txt 

tail -n +2  ukb41465.txt | cut -f 2- | join  -t $'\t' - UKBB.470Kexome.IDs.txt >> ukb41465.470Kexomes.txt 

head -n 1 ukb43357.txt > ukb43357.470Kexomes.txt 

tail -n +2  ukb43357.txt | cut -f 2- | join  -t $'\t' - UKBB.470Kexome.IDs.txt >> ukb43357.470Kexomes.txt 

# the original files have an unwanted tab before the eid field but not in the first line
# (they also have an extra tab at the end, which this does not fix)

# head -n 1 ukb51119.40PCs.txt > UKBB.exome.40PCs.20201103.txt
# grep -f UKBB.exome.IDs.20201103.txt ukb51119.40PCs.txt >> UKBB.exome.40PCs.20201103.txt
# cut UKBB.exome.40PCs.20201103.txt -f 1-21 -d " " > UKBB.exome.20PCs.20201103.txt