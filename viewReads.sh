#!/bin/bash

igv=/home/rejudcu/igv/IGV_2.9.4/igv.sh
samtools=/share/apps/genomics/samtools-1.9/bin/samtools
bamReadcount=/home/rejudcu/bam-readcount/build/bin/bam-readcount
targetDir=/home/rejudcu/UKBB/LOF/snapShots
refDir=/home/rejudcu/reference38

ID=$1
chr=$2
pos=$3

if [ ! -e /SAN/ugi/UGIbiobank/data/downloaded/cramFiles/${ID}_23153_0_0.cram ]
then
  pushd /SAN/ugi/UGIbiobank/data/downloaded
  ./ukbfetch -e$ID -d23153_0_0
  ./ukbfetch -e$ID -d23154_0_0
  mv ${ID}_23154_0_0.crai ${ID}_23153_0_0.crai
  mv ${ID}_23153_0_0.* cramFiles
  popd
fi

name=$ID.$chr.$pos

if [ ! -e $targetDir/$name.png ]
then
mkdir $targetDir/$name
echo load /SAN/ugi/UGIbiobank/data/downloaded/cramFiles/${ID}_23153_0_0.cram >$name.inp
echo goto chr$chr:$pos >> $name.inp
echo snapshotDirectory $targetDir/$name >> $name.inp
echo snapshot >> $name.inp
echo exit >> $name.inp

xvfb-run -d $igv -b $name.inp 

mv $targetDir/$name/*.png $targetDir/$name.png
fi

if [ ! -e $targetDir/$name.counts.txt ]
then
  mkdir $targetDir/$name
  $samtools view -b -o $targetDir/$name/$name.bam /SAN/ugi/UGIbiobank/data/downloaded/cramFiles/${ID}_23153_0_0.cram chr$chr:$(( pos - 1000 ))-$(( pos + 1000 )) 
  $samtools index $targetDir/$name/$name.bam 
  $bamReadcount -f $refDir/chr$chr.fa $targetDir/$name/$name.bam chr$chr:$pos-$pos > $targetDir/$name/$name.counts.txt
  read -a words < $targetDir/$name/$name.counts.txt
  counts=""
  for w in {4..9}
  do
    IFS=':', read -a word <<<  ${words[$w]}
	counts="$counts ${word[0]}:${word[5]}/${word[6]}"
  done
  echo $counts > $targetDir/$name.countsByStrand.txt
fi

rm -rf $targetDir/$name
# last line 
