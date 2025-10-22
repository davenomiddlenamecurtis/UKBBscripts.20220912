#!/share/apps/R-3.6.1/bin/Rscript

# script to run intVarAssoc on a number of intervals

# Tasks
# Create a VCF for each interval
# Produce annotations with GPN for that interval and upload them
# Run intVarAssoc on each interval 

wd="/home/rejudcu/UKBB/LDLRpromoter.20250814"
Model="UKBB.HL.nonCoding.20250814"
Model="UKBB.HL.nonCoding.analysePrevious.20251017"
Model="UKBB.HL.nonCoding.analyseWithPrevious.20251017"
IntervalsFile="LDLRnonCoding.20251017.txt"

Model="UKBB.HL.nonCoding.analysePrevious.20251017"
IntervalsFile="LDLRallNonCoding.20251019.txt"

Model="UKBB.HL.nonCoding.analyseWithPrevious.20251017"
IntervalsFile="LDLRnonCoding.20251017.txt"

setwd(wd)
args = commandArgs(trailingOnly=TRUE)
if (length(args)>=2) {
	Model=args[1]
	IntervalsFile=args[2]
}


if (!file.exists(sprintf("/cluster/project9/bipolargenomes/UKBB/%s/results",Model))) {
	system(sprintf("mkdir /cluster/project9/bipolargenomes/UKBB/%s",Model))
	system(sprintf("mkdir /cluster/project9/bipolargenomes/UKBB/%s/results",Model))
}

GPNScoreFile="/cluster/ref1/UGISharedData/GPN-MSA/scores.tsv.bgz"

Intervals=data.frame(read.table(IntervalsFile,header=FALSE,stringsAsFactors=FALSE))

AllDone=TRUE
for (r in 1:nrow(Intervals)) {
	Name=Intervals[r,1]
	VCF=sprintf("%s.vcf.gz",Name)
	system(sprintf("rm %s.exists.txt",VCF))
	system(sprintf("dx cd /WGS/VCFs;f=`dx ls %s`; if [ .$f == .%s ]; then echo got it > %s.exists.txt; fi",VCF,VCF,VCF))
	if (!file.exists(sprintf("%s.exists.txt",VCF))) {
		AllDone=FALSE
		system(sprintf("bash /home/rejudcu/UKBB/RAPfiles/RAPscripts/makeExtractedVCF.20250919.sh %s chr%s %d %d",Intervals[r,1],Intervals[r,2],Intervals[r,3],Intervals[r,4]))
		ScriptName=sprintf("extract.%s.sh",Name)
		system(sprintf("dx cd /scripts ; dxupload %s",ScriptName))
		CommStr=sprintf("dx cd /WGS/VCFs ; dx run swiss-army-knife -y --ignore-reuse --instance-type mem3_hdd2_v2_x4  -imount_inputs=FALSE -iin=/scripts/%s -icmd=\"bash %s\"",ScriptName,ScriptName)
		print(CommStr)
		system(CommStr)
	}
}

if (AllDone!=TRUE) {
	quit()
}

for (r in 1:nrow(Intervals)) {
	Name=Intervals[r,1]
	AnnotVCF=sprintf("%s.GPN.annot.vcf.gz",Name)
	system(sprintf("rm %s.exists.txt",AnnotVCF))
	system(sprintf("dx cd /annot;f=`dx ls %s`; if [ .$f == .%s ]; then echo got it > %s.exists.txt; fi",AnnotVCF,AnnotVCF,AnnotVCF))
	if (!file.exists(sprintf("%s.exists.txt",AnnotVCF))) {
		AllDone=FALSE
		system(sprintf("rm %s*",AnnotVCF))
		system(sprintf("tabix %s %s:%d-%d > %s.GPN.annot.vcf",GPNScoreFile,Intervals[r,2],Intervals[r,3],Intervals[r,4],Name))
		system(sprintf("bgzip %s.GPN.annot.vcf",Name))
		system(sprintf("tabix -p vcf %s",AnnotVCF))
		system(sprintf("dx cd /annot; dxupload %s",AnnotVCF))
		system(sprintf("dx cd /annot; dxupload %s.tbi",AnnotVCF))
	}
}

for (r in 1:nrow(Intervals)) {
	Name=Intervals[r,1]
	IntFile=sprintf("%s.int",Name)
	cat(sprintf("%s:%s-%s\n",Intervals[r,2],Intervals[r,3],Intervals[r,4]),file=IntFile,append=FALSE)
	system(sprintf("dx cd /intervals; dxupload %s",IntFile))
}

AllDone=TRUE
for (r in 1:nrow(Intervals)) {
	Name=Intervals[r,1]
	SaoFile=sprintf("%s.%s.sao",Model,Name)
	if (!file.exists(sprintf("/cluster/project9/bipolargenomes/UKBB/%s/results/%s",Model,SaoFile))) {
		CommStr=sprintf("dx cd /results/%s/intervalResults/;pushd /cluster/project9/bipolargenomes/UKBB/%s/results;dx download %s;popd",Model,Model,SaoFile)
		print(CommStr)
		system(CommStr)
		if (!file.exists(sprintf("/cluster/project9/bipolargenomes/UKBB/%s/results/%s",Model,SaoFile))) {
			AllDone=FALSE
			CommStr=sprintf("cp /home/rejudcu/UKBB/RAPfiles/fileLists/for.%s.lst for.%s.lst; echo /reference38/chr%s.fa>> for.%s.lst;dx cd /fileLists; dxupload for.%s.lst",Model,Name,Intervals[r,2],Name,Name)
			# I think this is only for CpG islands
			print(CommStr)
			system(CommStr)
			CommStr=sprintf("dx cd /;dx run swiss-army-knife -y --ignore-reuse --instance-type mem3_ssd2_v2_x4 -imount_inputs=FALSE -iin=\"/scripts/runOneIntervalWithFileList.20250418.sh\" -icmd=\"bash runOneIntervalWithFileList.20250418.sh %s %s for.%s.lst\" ",Model,Name,Name)
			print(CommStr)
			system(CommStr)
		}
	}
}

if (AllDone!=TRUE) {
	quit()
}

for (r in 1:nrow(Intervals)) {
	Name=Intervals[r,1]
	SaoFile=sprintf("%s.%s.sao",Model,Name)
	CommStr=sprintf("cp /cluster/project9/bipolargenomes/UKBB/%s/results/%s .",Model,SaoFile)
	system(CommStr)
}
