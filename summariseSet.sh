#!/bin/bash
setName=$1
setFileName=`ls *.$1.*`

nGenes=`head -n 4 $setFileName | tail -n 1 | cut -f 3 -d " "`
MLP=`head -n 6 $setFileName | tail -n 1 | cut -f 2`

echo $setName$'\t'$MLP$'\t'$nGenes
tail -n 8 | sort -k 2 -n -r | head 10

