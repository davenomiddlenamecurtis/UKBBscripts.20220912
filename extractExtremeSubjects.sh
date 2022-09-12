#!/bin/bash

# obtain summary clinical pictures of subjects with extreme scores for a given gene

# this will fail if a subject has more than one LOF variant

gene=$1


# exomesFile=/SAN/ugi/UGIbiobank/data/downloaded/ukb41465.exomes.txt
model=LOF.20201227
argFile=~/pars/gva.UKBB.$model.arg
exomesFile=/SAN/ugi/UGIbiobank/data/downloaded/ukb41465.exomes.20201103.txt



saoFile=$model.$gene.sao
if [ ! -e $saoFile ]
then
  geneVarAssoc --arg-file $argFile --gene $gene --test-name $model
fi
# ID pheno score
# cat $saoFile | awk -v t=$threshold '{ if ((t>0 && $3>t ) || (t<0 && $3<t )) print $1 }' > IDs.$gene.txt
cat $saoFile | awk '
{ 
if ($1=="Controls") exit; 
if ($17 != "" && $17!="comment")
 { 
 split($1,l,":");
 split(l[2],p,"-");
 split(p[1],pp,"."); 
 split(p[2],a,","); 
 if ($2>$6)  print l[1] ":" pp[1],a[1],a[2],$1,"outputAlt"
 else print l[1] ":" pp[1],a[1],a[2],$1,"outputRef"
 } 
}' > $model.$gene.vars.txt
# with split loci may be . in position

rm IDs.$model.$gene.txt
rm IDsAndVars.$model.$gene.txt
rm IDsAndAnnots.$model.$gene.txt
if [ ! -s $model.$gene.vars.txt ]
then
  exit
fi
cat $model.$gene.vars.txt | while read a b c d e
do
  if [ ! -e altSubs.$model.$gene.$a.$b.$c.txt ]
  then
    if [ $e == outputRef ]
	then
      showAltSubs --arg-file $argFile --position $a --ref-all $b --alt-all $c --test-name altSubs.$model.$gene.$a.$b.$c --output-ref 1
	else
      showAltSubs --arg-file $argFile --position $a --ref-all $b --alt-all $c --test-name altSubs.$model.$gene.$a.$b.$c
	fi
  fi
  cat altSubs.$model.$gene.$a.$b.$c.txt | awk -v v=$d '{ print $1, v }' >> IDsAndAnnots.$model.$gene.txt
  cat altSubs.$model.$gene.$a.$b.$c.txt | while read s rest
  do
   echo $s >> IDs.$model.$gene.txt
   echo $s$'\t'$a-$b >> IDsAndVars.$model.$gene.txt
  done
done 
cat IDsAndAnnots.$model.$gene.txt | sort > IDsAndAnnots.sorted.$model.$gene.txt
cat IDs.$model.$gene.txt  | sort > IDs.sorted.$model.$gene.txt
echo eid$'\t'extra > extraField.txt
head -n 1 $exomesFile | join -t  $'\t' - extraField.txt > $model.$gene.phenos.txt
# this is because the exomes file has got a tab appended to each line after the header
# and otherwise fread will not work
tail -n +2 $exomesFile | join -t  $'\t' IDs.sorted.$model.$gene.txt - >> $model.$gene.phenos.txt 
# grep -f  IDs.sorted.$model.$gene.txt $exomesFile >> $model.$gene.phenos.txt # this is not really safe

# then R can read in the two files and get the column headings right
# can use IDs file to get a phenotype of who has any LOF variants for other analyses




 