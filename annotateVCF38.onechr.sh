#!/bin/bash
# for single chromosome files e.g. UKBexomeOQFE_chr11.vcf
# I think these can be run in parallel
set -x

root=$1
chr=$2
mult=no

if [ -z "$mult" ]
then
  mult=yes
# annotations for multiple transcripts
fi

if [ -z "$root" ]
then
  echo Usage: $0 vcfRoot chr [ target is vcfRootchr.vcf.gz, can also first export mult=no for only one transcript and X=no if no X chromosome data ]
  exit
fi

if [ ! -e $root$chr.vars.vcf ]
then
      zcat $root$chr.vcf.gz | cut -f1-9 > $root$chr.vars.vcf
fi

PICK=
MULT=.mult

if [ "$mult" == no ]
then
  PICK="--pick_allele_gene"
  MULT=
fi

# export PERL5LIB=/usr/lib/perl64:$PERL5LIB
OLDPERL5LIB=$PERL5LIB
export PERL5LIB=/usr/lib/perl64
chrfile=$root$chr$MULT.annot.vcf
if [ ! -e $chrfile.gz ]
then
# cat $root$chr.vars.vcf | perl /share/apps/ensembl-vep-97/vep \
#     --synonyms ~/vep/chr_synonyms.txt \
# 	--cache --dir /cluster/project9/bipolargenomes/vepcache --merged --port 3337 --force_overwrite \
# 	--sift b --polyphen b --offline --assembly GRCh38 --format vcf \
# 	--fasta /cluster/project9/bipolargenomes/vepcache/homo_sapiens_merged/97_GRCh38 \
# 	--vcf --output_file $chrfile $PICK
  cat $root$chr.vars.vcf | perl /share/apps/ensembl-vep-97/vep \
    --synonyms ~/vep/chr_synonyms.txt \
	--cache --dir /cluster/project9/bipolargenomes/vepcache --merged --force_overwrite \
	--sift b --polyphen b --assembly GRCh38 --format vcf \
	--fasta /cluster/project9/bipolargenomes/vepcache/homo_sapiens_merged/97_GRCh38 \
	--canonical --regulatory \
	--vcf --output_file $chrfile $PICK
  if [ -s ${chrfile}_summary.html ]
  then
    bgzip $chrfile
  else
    echo  ${chrfile}_summary.html is zero length
  fi
fi

root=$1
useX=yes
if [ "$useX" == "yes" ]
then
 allChrs="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X"
else
if [ "$use23" == "yes" ]
then
 allChrs="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23"
else
 allChrs="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"
fi
fi

vcfs=
allDone=yes
for c in $allChrs
do
chrfile=$root$c$MULT.annot.vcf
vcfs="$vcfs $chrfile.gz"
if [ ! -e $chrfile.gz ]
then
  allDone=no
fi
done

if [ $allDone == no ]
then
  exit
fi

# if this is the last file and all others are done then join them

# $bcftools concat -o $root$MULT.annot.vcf $vcfs
PERL5LIB=/usr/lib/perl64:/share/apps/genomics/vcftools-0.1.13/lib/perl5/site_perl/
VCFCONCAT=/share/apps/genomics/vcftools-0.1.13/bin/vcf-concat
$VCFCONCAT $vcfs > $root$MULT.annot.vcf
rm $root$MULT.annot.vcf.gz # just in case
bgzip $root$MULT.annot.vcf
tabix -f -p vcf $root$MULT.annot.vcf.gz
