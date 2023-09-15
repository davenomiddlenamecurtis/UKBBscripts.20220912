#!/bin/bash
# note that the column number to provide is one higher than that given in http://www.davecurtis.net/UKBB/ukb41465.html

tableFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb41465.470Kexomes.txt"
 
if [ .$2 == . ]
then
  echo Usage:
  echo $0 varName colNum
  exit
fi

echo IID$'\t'$1 > UKBB.$1.txt
tail -n +2 $tableFile | cut -f 1,$2 >> UKBB.$1.txt
