#!/bin/bash
if [ .$2 == . ]
then
  echo Usage: $0 tsvFile date
  exit
fi

# covid19_result.20200428.txt

tsvFile=$1
date=$2	
outFile=testResult.$date.txt
testFile=hadTest.$date.txt
positiveFile=testedPositive.$date.txt
testPerformedFile=testPerformed.$date.txt

tail -n +2 $tsvFile | sort | cut -f1,6 > rawResults.txt

rm -f tempfile.txt 
cat rawResults.txt | while read i p ; do echo $i 1 >> tempfile.txt ; done
cat tempfile.txt | sort | uniq > $testFile 

rm -f tempfile.txt 
cat rawResults.txt | while read i p ; do if [ $p == 1 ] ; then echo $i $p >> tempfile.txt; fi; done
cat tempfile.txt | sort | uniq > $positiveFile

echo IID testResult > $outFile
cat $positiveFile >> $outFile
tail -n +2 /home/rejudcu/vcf/UKBB/UKBB.BMI.txt | join -v1 - $positiveFile > notPositive.txt
cat notPositive.txt | while read i p; do echo $i 0 >> $outFile; done

echo IID hadTest > $testPerformedFile
cat $testFile >> $testPerformedFile
tail -n +2 /home/rejudcu/vcf/UKBB/UKBB.BMI.txt | join -v1 - $testFile > notTested.txt
cat notTested.txt| while read i p; do echo $i 0 >> $testPerformedFile; done
