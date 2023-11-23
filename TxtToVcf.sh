root=ukb23158_500k_OQFE.annotations

TOVCF='
BEGIN { FS=":"; OFS="\t" }
{ print $1, $2, $1 ":" $2 ":" $3 ":" $4, $3, $4, ".", "." }
'

 zcat $root.txt.gz | cut -f 1 -d " " | awk "$TOVCF" - > $root.vcf

allChrs="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X"
allChrs="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23"

for chr in $allChrs
do
 awk -v chr=$chr '{ if ($1==chr) {print $0 } }' $root.vcf > $root$chr.vars.vcf
 subComm.sh ~/UKBB/UKBBscripts.20220912/annotateVCF38.onechr.20210225.sh $root $chr
done
