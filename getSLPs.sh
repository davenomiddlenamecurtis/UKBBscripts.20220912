#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l tmem=1G,h_vmem=1G
#$ -l h_rt=20:0:0
#$ -V
#$ -R y
#$ -e /home/rejudcu/tmp
#$ -o /home/rejudcu/tmp

# model="ct08 ct08.rare"
# model=ExAC.ct08.rare
# model=1000G.ct08.rare
# disease=BP
# disease=WKS
# disease=IBDAJ
# disease="MIGen" # this MUST be the disease because otherwise copyVCF does not work
# disease=ADSP2
# model=common.withAPOE

disease=UKBB
# model="HL.all.20201103 HL.all.20201103.XGenes HL.all.20201103.notXGenes"
model="sex.all.20201111 HL.withSex.20201207 Depn.withSex.20201207 LOAD.withSex.20201208"
model="sex.all.20201111.XGenes sex.all.20201111.notXGenes"
# disease="UCLEx.Prionb2"
model="alcHigh.withSex.20201229"
model=HL.all.20201231
model=HL.withSex.20210102
model="alcProb.withSex.20210105"

if [ .$disease == . ]
then
  disease=Ashkenazi_ASJ
fi

if [ -z "$model" -o -z "$disease" ]
then
	echo need to set disease and model
	exit
fi

for d in $disease
do
for m in $model
do 

testName=$d.$m

if [ .$testName == . ]
then
	echo Need to set testName
	exit
fi
# testName=BPrare
 
workFolder=/cluster/project9/bipolargenomes/$d/$testName
# workFolder=/cluster/project9/bipolargenomes/SCHEMA/results/$d/$testName
resultsFolder=$workFolder/results
summFile=$workFolder/$testName.summ.txt

# NB must never have two tabs at start of line because will break you out of quoted segment

getSLPs='
BEGIN { ORS=""; nSLP=0; } 
{
if ($1 == "SLP" || $1 == "tSLP" || $1 == "tMLP" || $1 == "linrSLP" || $1 == "lrSLP") 
	{
	nSLP=nSLP+1;
	SLPs[nSLP]=$3;
	}
}
END {
    if (nSLP > 0) {
	print gene "\t";
	for (i=1; i<=nSLP; ++i) {
	print SLPs[i] "\t";
	}
	print "\n";
	}
}
'

# echo Gene$'\t'SLPD$'\t'SLPR$'\t'SLPHA$'\t'SLPHO > $summFile
# echo Gene$'\t'SLP$'\t'tSLPscore$'\t'tSLPscorePC$'\t'tSLPscorePCPRS$'\t'tSLPscorePCCNV$'\t'tSLPscoreALL> $summFile
echo Gene$'\t'tSLP$'\t'lrSLPwithPCs$'\t'lrSLP> $summFile
# echo Gene$'\t'SLP$'\t'tSLP$'\t'tSLPPC> $summFile

find  $resultsFolder -name '*.sao' | while read resultsFile
	do
	gene=${resultsFile%.sao}
	gene=${gene#*$m.}
	tail -n 100 $resultsFile | grep LP | awk -v gene=$gene "$getSLPs" >> $summFile # all valid lines contain LP
#	awk -v gene=$gene "$getSLPs" $resultsFile >> $summFile
	done
done
done
