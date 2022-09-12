#!/bin/bash
# for single chromosome files e.g. UKBexomeOQFE_chr11.vcf
# I think these can be run in parallel
set -x

root=$1
chr=$2
mult=no
filterNonCoding=yes # ignore non-coding transcripts, pseudogenes, etc but keep RNA genes

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

  cat $root$chr.vars.vcf | /share/apps/ensembl-vep-97/vep \
    --synonyms ~/vep/chr_synonyms.txt \
	--cache --dir /cluster/project9/bipolargenomes/vepcache --merged --force_overwrite \
	--sift b --polyphen b --assembly GRCh38 --format vcf \
	--fasta /cluster/project9/bipolargenomes/vepcache/homo_sapiens_merged/97_GRCh38 \
	--canonical --regulatory --use_given_ref \
	--custom /home/rejudcu/loftee/loftee/gerp_conservation_scores.homo_sapiens.GRCh38.bw,GERP,bigwig \
	--vcf --output_file $chrfile $PICK
# --use_given_ref is new, intended to fix e.g. rs1229984
  if [ -s ${chrfile}_summary.html ]
  then
  if [ $filterNonCoding == yes ]
  then
    rm -f $root$chr$MULT.annot.unfiltered.vcf
	mv $chrfile $root$chr$MULT.annot.unfiltered.vcf
	/share/apps/ensembl-vep-97/filter_vep \
	  --input_file $root$chr$MULT.annot.unfiltered.vcf \
	  --format vcf \
	  --filter "biotype in protein_coding,IG_C_gene,IG_D_gene,IG_J_gene,IG_V_gene,lncRNA,miRNA,misc_RNA,Mt_rRNA,Mt_tRNA,ribozyme,rRNA,scaRNA,scRNA,snoRNA,snRNA,sRNA,TEC,TR_C_gene,TR_D_gene,TR_J_gene,TR_V_gene,vault_RNA" \
	  --output_file $chrfile 
  fi
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

# Unfortunately we don't have an easily accessible way to get this in the right format. A rather roundabout way to get this is to copy/paste from BioMart.
# This is the list you'll get from there:
# IG_C_gene , IG_C_pseudogene , IG_D_gene , IG_J_gene , IG_J_pseudogene , IG_pseudogene , 
# IG_V_gene , IG_V_pseudogene , lncRNA , miRNA , misc_RNA , Mt_rRNA , Mt_tRNA , nonsense_mediated_decay , non_stop_decay , 
# polymorphic_pseudogene , processed_pseudogene , processed_transcript , protein_coding , pseudogene , retained_intron , 
# ribozyme , rRNA , rRNA_pseudogene , scaRNA , scRNA , snoRNA , snRNA , sRNA , TEC , transcribed_processed_pseudogene , 
# transcribed_unitary_pseudogene , transcribed_unprocessed_pseudogene , translated_processed_pseudogene , translated_unprocessed_pseudogene , 
# TR_C_gene , TR_D_gene , TR_J_gene , TR_J_pseudogene , TR_V_gene , TR_V_pseudogene , unitary_pseudogene , unprocessed_pseudogene , vault_RNA

