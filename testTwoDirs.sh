#!/bin/bash
unset DX_WORKSPACE_ID
dx cd $DX_PROJECT_CONTEXT_ID:

mkdir results
mkdir results/dir1
mkdir results/dir2

echo $1 > results/dir1/$1.txt
echo $2 > results/dir2/$2.txt
