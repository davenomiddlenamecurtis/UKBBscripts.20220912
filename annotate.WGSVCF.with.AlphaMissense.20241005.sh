#!/bin/bash
# for e.g. single gene files e.g. OMA1.exons.vars.vcf

set -x

root=$1 # e.g. OMA1.exons.vars
mult=no # set --pick_allele_gene

if [ -z "$mult" ]
then
  mult=yes
# annotations for multiple transcripts
fi

if [ -z "$root" ]
then
  echo Usage: $0 vcfRoot [ target is vcfRoot.vcf, can also first export mult=no for only one transcript ]
  exit
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
outfile=$root$MULT.AM.annot.vcf
if [ ! -e $outfile.done ]
then
  perl /share/apps/ensembl-vep-97/vep \
	--input_file $root.vcf \
    --synonyms ~/vep/chr_synonyms.txt \
	--cache --dir /cluster/project9/bipolargenomes/vepcache --merged --force_overwrite \
	--sift b --polyphen b --assembly GRCh38 --format vcf \
	--fasta /cluster/project9/bipolargenomes/vepcache/homo_sapiens_merged/97_GRCh38 \
	--canonical --regulatory \
	--plugin AlphaMissense,file=/share/ref/VEP/AlphaMissense/AlphaMissense_hg38.tsv.gz,cols=all \
	--vcf --output_file $outfile $PICK
  if [ -s ${outfile}_summary.html ]
  then
    bgzip $outfile
	tabix -p vcf $outfile.gz
	echo done > $chrfile.done
  else
    echo  ${chrfile}_summary.html is zero length
  fi
fi
