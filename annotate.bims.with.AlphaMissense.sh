#!/bin/bash
# for single chromosome files e.g. ukb23158_c9_b0_v1.bim
# I think these can be run in parallel
# file names use X but internally use 23

set -x

root=$1
chr=$2
mult=no # set --pick_allele_gene
useX=yes

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

TOVCF='
BEGIN { FS=":"; OFS="\t" }
{ print $1, $2, $1 ":" $2 ":" $3 ":" $4, $3, $4, ".", "." }
'

if [ ! -e $root.$chr.vars.vcf ]
then
	cut -f 2 ukb23158_c${chr}_b0_v1.bim | awk "$TOVCF" - > $root.$chr.vars.vcf
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
export PERL5LIB=/usr/lib/perl64:/share/apps/ensembl-vep-97/plugins/
chrfile=$root.$chr$MULT.AM.annot.vcf
if [ ! -e $chrfile.done ]
then
  perl /share/apps/ensembl-vep-97/vep \
	--input_file $root.$chr.vars.vcf \
    --synonyms ~/vep/chr_synonyms.txt \
	--cache --dir /cluster/project9/bipolargenomes/vepcache --merged --force_overwrite \
	--sift b --polyphen b --assembly GRCh38 --format vcf \
	--fasta /cluster/project9/bipolargenomes/vepcache/homo_sapiens_merged/97_GRCh38 \
	--canonical --regulatory \
	--plugin AlphaMissense,file=/share/ref/VEP/AlphaMissense/AlphaMissense_hg38.tsv.gz,cols=all \
	--vcf --output_file $chrfile $PICK
  if [ -s ${chrfile}_summary.html ]
  then
    bgzip $chrfile
	echo done > $chrfile.done
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
for chr in $allChrs
do
chrfile=$root.$chr$MULT.AM.annot.vcf
vcfs="$vcfs $chrfile.gz"
if [ ! -e $chrfile.done ]
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
$VCFCONCAT $vcfs > $root$MULT.AM.annot.vcf
rm $root$MULT.AM.annot.vcf.gz # just in case
bgzip $root$MULT.AM.annot.vcf
tabix -f -p vcf $root$MULT.AM.annot.vcf.gz
