#/bin/bash

badIDList=/home/rejudcu/UKBB/removedSubs.20210202.txt

if [ .$1 == . ]
then
 echo $0 usage:
 $0 fileToReplaceBadIDsIn
fi

fileToReplaceBadIDsIn=$1

badIDStart=20

cat $badIDList | while read ID
do
  sed -i s/$ID/-$badIDStart/g $fileToReplaceBadIDsIn
  badIDStart=$(( badIDStart + 1 ))
done

