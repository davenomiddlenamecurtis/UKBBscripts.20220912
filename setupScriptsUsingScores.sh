#!/bin/bash
# DC script to set up GVA analyses, one script per gene

# geneList=/home/rejudcu/reference/allGenes140817.onePCDHG.txt
# geneList=/home/rejudcu/reference38/allGenes.20191018.onePCDHG.txt
geneList=/home/rejudcu/reference38/genes.with.vars.20201225.txt
# geneList=/home/rejudcu/reference/DRDgenes.txt
# disease=MPexomes
# model=bp1.myWeights

# disease=ADSP2
# model=common.withAPOE

disease=UKBB
sourceModel=HL.all.20201103
model="sex.all.20201111 HL.withSex.20201207 Depn.withSex.20201207 LOAD.withSex.20201208"
model="alcHigh.withSex.20201229"
model=HL.withSex.20201207

refdir=reference38

homeFolder=/cluster/project9/bipolargenomes
argFolder=/home/rejudcu/pars
softwareFolder=/home/rejudcu/bin
dataHome=/home/rejudcu

if [ -z "$disease" -o -z "$model" ]
then
	echo Error in $0: must set environment variables disease and model
	exit
fi

someGenesLeft=no

for d in $disease
do
for m in $model
do

testName=$d.$m
argFile=$argFolder/sco.$testName.args
# argFileTxt=`cat $argFile`


# workFolder=/cluster/project8/bipolargenomes/GVA

workFolder=$homeFolder/$d/$testName
inputFolder=$homeFolder/$d/$d.$sourceModel/results
mkdir $homeFolder/$d
mkdir $workFolder

nSplits=100

splitScript=$workFolder/scripts/split${nSplits}s.sh
scriptName=$testName.runSplit${nSplits}.sh
mainSplitScript=$workFolder/scripts/$scriptName

qdel $testName.'runSplit*'

nhours=4
vmem=6 
memory=2
queue=queue6
scratch=0

# UKBB analyses were running out of memory with vmem=6 with geneVarAssoc
# but maybe these will need much less - nope, failed with 2
vmem=8

if [ ! -e $workFolder ]; then mkdir $workFolder; fi;
wastebin=$workFolder/wastebin
if [ ! -e $wastebin ]; then mkdir $wastebin; fi
if [ ! -e $workFolder/results ]; then mkdir $workFolder/results; fi;
if [ -e $workFolder/error ]; then mv $workFolder/error $wastebin/error; ( rm -r $wastebin/error & ) ; fi;
mkdir $workFolder/error
if [ -e $workFolder/scripts ]; then mv $workFolder/scripts $wastebin/scripts; (rm -r $wastebin/scripts & ); fi;
mkdir $workFolder/scripts; 
if [ -e $workFolder/temp ]; then mv $workFolder/temp $wastebin/temp; (rm -r $wastebin/temp & ); fi;
mkdir $workFolder/temp; 

cat $geneList | while read geneName
    do
    outFile=$workFolder/results/$testName.$geneName.sao
    if [ ! -e $outFile ]
    then 
		shellScript=$workFolder/scripts/runScoreassoc.$testName.$geneName.sh
		inputScoreFile=$inputFolder/$d.$sourceModel.$geneName.sco
		commLine="scoreassoc --argfile $argFile --inputscorefile $inputScoreFile --outfile $testName.$geneName.sao"
		echo "rm $testName.$geneName.sao
		pwd
		if [ -e $inputScoreFile ]
		then
		  echo Running:
		  echo $commLine
		  $commLine 
		  echo finished running scoreassoc
		else
		  echo No score file $inputScoreFile > $testName.$geneName.sao
		fi
		cp $testName.$geneName.sao $outFile 
		if [ ! -s $outFile ] ; then rm -f $outFile ; fi
		" >> $shellScript
    fi
    done
	
nScriptsWritten=`find $workFolder/scripts -name 'runScoreassoc.*.sh' | wc -l`
if [ $nScriptsWritten -lt $nSplits ]
then 
	nSplits=$nScriptsWritten
fi 

if [ -e  $mainSplitScript ] ; then rm  $mainSplitScript; fi

echo "
#!/bin/bash
#$ -S /bin/bash
#$ -e $workFolder/error
#$ -o $workFolder/error
#$ -l tscr=${scratch}G
#$ -l tmem=${vmem}G,h_vmem=${vmem}G
#$ -l h_rt=${nhours}:0:0
#$ -t 1-$nSplits
#$ -V
#$ -R y

# I used to have #$ -cwd but I am going to try just omitting it as sometimes cannot cd to it
# If that does not work may try -wd /scratch0

date
echo bash $splitScript \$SGE_TASK_ID
bash -x $splitScript \$SGE_TASK_ID
date
" > $mainSplitScript

echo "
#!/bin/bash
set +e
#  was exiting after running just one, possibly because no proper exit code from script
# this should switch off errexit
echo Running \$0 with argument \$1
$workFolder/temp
cd $workFolder/temp
myDir=\$RANDOM
mkdir \$myDir
cd \$myDir # this is all so I can have local vcf and reference folders so par files will work with this and with scratch0
mkdir temp
cd temp # so relative paths will work OK
n=1
find $workFolder/scripts -name 'runScoreassoc*sh' | while read f
do
if [ .\$n == .\$1 ]
then
	echo running source \$f # try using source $f instead of bash $f
	source \$f > \$f.log  2>&1
	cp \$f.log $workFolder/error
	echo finished running source \$f
	rm *
fi
if [ \$n -eq $nSplits ]
then
	n=1
else
	n=\$(( \$n + 1 ))
fi
done
cd ../..
rm -r \$myDir
" > $splitScript

count=`find $workFolder/scripts -name 'runScoreassoc*sh' | wc -l`

if [ $count -gt 0 ]
then
	echo wrote $count scripts
	echo qsub -N $scriptName $mainSplitScript
	pushd $workFolder
# reason for this is that I would get Eqw with qstat -j error message: error: can't chdir to /home/rejudcu/tmp: No such file or directory 
	qsub -N $scriptName $mainSplitScript
	popd
	someGenesLeft=yes
else
	echo No genes left to do for $testName	
fi

done
done
# the model etc. are currently hard-coded so exports below are irrelevant
logFile=${0##*/}
if [ $someGenesLeft = yes ]
then
	echo will schedule script to run again
	echo "export disease=\"$disease\"; export model=\"$model\"; export geneList=$geneList; bash $0 &> $workFolder/$logFile.log" | at now + $nhours hours
else
	echo date > $workFolder/$logFile.log
	echo All results files written OK >> $workFolder/$logFile.log
fi
