#!/bin/bash 
# make a script which will extract scores for all genes relevant to a phenotype

refGeneFile=~/reference38/refseqgenes.hg38.20191018.sorted.onePCDHG.txt
scoreFile=/cluster/ref1/UGISharedData/GPN-MSA/scores.tsv.bgz

prologue=~/UKBB/RAPfiles/RAPscripts/makeGeneWESVCF.prologue.20241027.sh

extractGeneCoords=' BEGIN { start=300000000; end=0 } { chr= $3; if ($5<start) start=$5; if ($6>end) end=$6 } END { print chr, start, end }'

rm coords.txt
rm extractedScores.vcf*
for gene in $@
do
  geneArgs=`grep -w $gene $refGeneFile | awk "$extractGeneCoords"`
  coords=($geneArgs)
  chr=${coords[0]/chr/} # removes chr, see https://stackoverflow.com/questions/19551613/modify-the-content-of-variable-using-sed-or-something-similar
  start=${coords[1]}
  end=${coords[2]}

  echo $chr $start $end >> coords.txt
done

sort -k1,1n -k2,2n < coords.txt > coords.sorted.txt

cat coords.sorted.txt | while read chr start end
do
  tabix $scoreFile $chr:$start-$end >> extractedScores.vcf
  bgzip -f extractedScores.vcf
  tabix -p vcf extractedScores.vcf.gz
done




