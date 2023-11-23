#!/bin/bash 

# script to run on virtual machine in RAP using Swiss Army Knife app

unset DX_WORKSPACE_ID
dx cd $DX_PROJECT_CONTEXT_ID:

annotFilePath="/Bulk/Exome sequences/Population level exome OQFE variants, PLINK format - final release/helper_files/"
annnotFile=ukb23158_500k_OQFE.annotations.txt.gz

mkdir workDir
pushd workDir
dx cd "$annotFilePath"
dx download $annnotFile
dx cd /annot
dx upload $annnotFile

popd
echo Transfer completed > ranProg.txt

