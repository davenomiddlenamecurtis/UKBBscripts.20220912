#!/bin/bash
 
cmd=$1
args=($1)
echo $cmd
echo $args

echo dx run swiss-army-knife -y -imount_inputs=FALSE -iin=${args[0]} -iin=$2 -iin=/scripts/DownloadNeededFiles.sh -icmd="bash DownloadNeededFiles.sh; bash $1" $3 $4 $5
# or might need:
# echo dx run swiss-army-knife -y -imount_inputs=FALSE -iin=${args[0]} -iin=$2 -iin=/scripts/DownloadNeededFiles.sh -icmd=\"bash DownloadNeededFiles.sh; bash $1\" $3 $4 $5

# RunOnRap.sh "ThisScriptFile.sh arg1 arg2" TheseAreNeededFiles.txt -these-flags-get-appended -such-as --instance-type mem3_hdd2_v2_x4  