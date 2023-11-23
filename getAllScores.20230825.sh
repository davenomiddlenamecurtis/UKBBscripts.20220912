#!/bin/bash 
set -x
# get scores for all genes in a list using a model and file list of extra needed files

if [ .%1 == . ]
then
	echo Usage: $0  model geneListFile otherNeededFileList
	exit
fi

testName=$1
geneList=$2
otherNeededFileList=$3

allChrs="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X"

resultsDir=/cluster/project9/bipolargenomes/UKBB/$testName/results
failedDir=/cluster/project9/bipolargenomes/UKBB/$testName/failed
mkdir /cluster/project9/bipolargenomes/UKBB/$testName
mkdir $resultsDir
mkdir $failedDir

if [ ! -e $resultsDir/$testName.scoresSoFar.txt ]
then
	echo GenesWithScores > $resultsDir/$testName.scoresSoFar.txt
fi

if [ ! -e $failedDir/$testName.failedSoFar.txt ]
then
	echo GenesWhichFailed > $failedDir/$testName.failedSoFar.txt
fi

pushd $resultsDir
dx mkdir /results
dx mkdir /results/$testName/
dx mkdir /results/$testName/geneResults/
dx mkdir /results/$testName/failed/
# important that these folders exist so that dx cd will be successful

dx cd /results/$testName/geneResults/
dx ls  > $testName.onDx.txt
grep '.sco' $testName.onDx.txt | while read f
do
	dx download $f -f
	if [ -e $f -a -s $f ]
	then
		g=${f%.sco}
		h=${g##*.}
		echo $h >> $testName.scoresSoFar.txt
		dx rm $f
	fi
done
# this will only work if the gva output files are named as expected
grep '.sao' $testName.onDx.txt | while read f
do
	dx download $f -f
	if [ -e $f -a -s $f ]
	then
		dx rm $f
	fi
done
cd $failedDir
dx cd /results/$testName/failed
dx ls  > $testName.onDx.txt
cat $testName.onDx.txt | while read f
do
	dx download $f -f 
	if [ -e $f -a -s $f ]
	then
		g=${f%.failed.txt}
		echo $g >> $testName.failedSoFar.txt
		dx rm $f
	fi
done
popd

cat $geneList | grep -w -v -f $resultsDir/$testName.scoresSoFar.txt | grep -w -v -f $failedDir/$testName.failedSoFar.txt >$testName.leftToDo.txt

for chr in $allChrs 
do
	rm $testName.toDo.$chr.txt
	grep -w -f ~/UKBB/RAP/reference38/AllGenes.$chr.txt $testName.leftToDo.txt > $testName.toDo.$chr.txt
done

for chr in $allChrs 
do
	if [ -e $testName.toDo.$chr.txt -a -s $testName.toDo.$chr.txt ]
	then
		dx cd /geneLists
		dx rm $testName.toDo.$chr.txt
		dx upload $testName.toDo.$chr.txt
		dx cd /
		dx run swiss-army-knife --name "$testName.toDo.$chr" -y --instance-type mem3_hdd2_v2_x4  -imount_inputs=FALSE -iin="/scripts/runGeneSetWithFileList.20230825.sh" -icmd="bash runGeneSetWithFileList.20230825.sh $testName $chr $testName.toDo.$chr.txt $otherNeededFileList"
	fi
done 

# invoke with e.g. 
# bash ~/UKBB/RAP/scripts/getAllScores.20230816.sh UKBB.HT.20230807 HT.genes.20230816.txt /phenos/UKBB.HT.20230816.txt
