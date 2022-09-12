#!/bin/bash
# DC script to set up GVA analyses, one script per gene

# geneList=/home/rejudcu/reference/allGenes140817.onePCDHG.txt
# geneList=/home/rejudcu/reference38/allGenes.20201224.onePCDHG.txt
geneList=/home/rejudcu/reference38/genes.with.vars.20201225.txt
# ONLY FUNCTIONAL VARIANTS - MAY NEED TO CHANGE THIS !!!
geneList=/home/rejudcu/reference38/genes.with.vars.func.20210103.txt
# geneList=/home/rejudcu/reference/DRDgenes.txt
# disease=MPexomes
# model=bp1.myWeights

# disease=ADSP2
# model=common.withAPOE

disease=UKBB
# model=T2D.func
model=HL.all.20201231

refdir=reference38

if [ -z $geneList ]
then
# geneList=/home/rejudcu/SSSDNMclinical/notNeuroGenes.lst
# geneList=/home/rejudcu/SSSDNMclinical/dominantGenes.lst
# geneList=/home/rejudcu/SSSDNMclinical/recessiveGenes.lst
echo geneList is not set
exit
fi
# geneList=/home/rejudcu/reference/DRDgenes.txt
# geneList=/home/rejudcu/tmp/FAM21EP.lst

# disease="UCLEx.Prionb2"
# model="ExAC.ct08.rare"
# model="ct08.cleaned"
# must be in this order or else qdel will delete all the ct08 jobs
if [ -z "$disease" ]
then
  disease=ADSP
fi
if [ -z "$model" ]
then
   model=all
#   model="codeVars.Dam codeVars.Dis"
#   model="codeVars.Dam.NotNeuro codeVars.Dis.NotNeuro"
#  model="codeRecDam codeRecDis"
fi

nhours=10 # till this script runs again, was 6 but this just interrupted running scripts
joblength=8 # time for job to run before it gets terminated, was 4 then 6
queue=queue6
scratch=100

# Trying back to 6 now that memory for subject loci is allocated dynamically
# UKBB analyses were running out of memory with vmem=6
# vmem=8
# I think one or two genes did not work with vmem=8
# yup, calloc() failed
# Assertion `sub[s]=(subject *)calloc(1,sizeof(subject))' failed.
# vmem=12
# Still not enough
# scoreassoc: ../src/scoreassoc.cpp:49: int main(int, char**): Assertion `sub[s]=(subject *)calloc(1,sizeof(subject))' failed.
# vmem=16

# actually, qstat -j says maxvmem is never over 1G
# vmem=2
# vmem=6 # just for last few genes 
# vmem=24 # 12 did not work for a few, this worked for all but KIAA1109
# vmem=32 # last three - of queues forever go back to 24 and find out what went wrong
# vmem=48 # just for TTN - failed with 34,000 loci, just said segmentation fault 
 vmem=60 # just for KIAA1109 and TTN


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
argFile=$argFolder/gva.$testName.arg


# workFolder=/cluster/project8/bipolargenomes/GVA

workFolder=$homeFolder/$d/$testName
mkdir $homeFolder/$d
mkdir $workFolder

nSplits=100

splitScript=$workFolder/scripts/split${nSplits}s.sh
scriptName=$testName.runSplit${nSplits}.sh
mainSplitScript=$workFolder/scripts/$scriptName

qdel $testName.'runSplit*'


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
		scoreFile=$workFolder/results/$testName.$geneName.sco
		shellScript=$workFolder/scripts/runGVA.$testName.$geneName.sh
		elogFile=$workFolder/results/$testName.$geneName.elog
# I may add an exclusion log file so I can find which variants failed which conditions
		echo "export LMOD_SH_DBG_ON=1
# always use scratch0 now there is a flatfile access - the cd /scratch0/uniquedir is done by the calling script
#		cd /scratch0
#		df /scratch0
		mkdir $geneName
		cd $geneName
		rm gva.$geneName.*
		hostname
		pwd
		commLine=\"geneVarAssoc --arg-file $argFile --gene $geneName --keep-temp-files 1\" 
		echo Running:
		echo \$commLine
		\$commLine 
		echo finished running geneVarAssoc
		cp *.$geneName.sco $scoreFile 
		cp *.$geneName.sao $outFile 
		cp *.$geneName.elog $elogFile
		nSLPs=\`grep SLP $outFile | wc -l\`
		if [ \$nSLPs == 0 ]
		then
		  ls -l
		  echo no result so will try to run scoreassoc on its own
		  cat *.$geneName.sh
		  bash -x *.$geneName.sh
		  cp *.$geneName.sco $scoreFile 
		  cp *.$geneName.sao $outFile 
		  cp *.$geneName.elog $elogFile
		fi
		nSLPs=\`grep SLP $outFile | wc -l\`
		if [ \$nSLPs == 0  ] ; then ls -l;  hostname; df /scratch0; echo still no output so deleting files; rm -f $outFile $scoreFile $elogFile; fi
		cd ..
		rm -rf $geneName
		" >> $shellScript
    fi
    done

# was \$commLine > gva.$testName.$geneName.elog 
		
nScriptsWritten=`find $workFolder/scripts -name 'runGVA.*.sh' | wc -l`
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
#$ -l h_rt=${joblength}:0:0
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
# cd $workFolder/temp
cd /scratch0
myDir=\$RANDOM
mkdir \$myDir
cd \$myDir # this is all so I can have local vcf and reference folders so par files will work with this and with scratch0
mkdir vcf
mkdir vcf/$disease
cd vcf/$disease
ln -s $dataHome/vcf/$disease/* .
cd ../..
mkdir $refdir
cd $refdir
ln -s $dataHome/$refdir/* .
cd ..
# mkdir temp
# cd temp # so relative paths will work OK
# the called script will make a directory called geneName and cd into it
n=1
find $workFolder/scripts -name 'runGVA*sh' | while read f
do
if [ .\$n == .\$1 ]
then
	echo running source \$f # try using source $f instead of bash $f
	source \$f
	echo finished running source \$f
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

count=`find $workFolder/scripts -name 'runGVA*sh' | wc -l`

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

logFile=${0##*/}
if [ $someGenesLeft = yes ]
then
	echo will schedule script to run again
	echo "export disease=\"$disease\"; export model=\"$model\"; export geneList=$geneList; bash $0 &> $workFolder/$logFile.log" | at now + $nhours hours
else
	echo date > $workFolder/$logFile.log
	echo All results files written OK >> $workFolder/$logFile.log
fi
