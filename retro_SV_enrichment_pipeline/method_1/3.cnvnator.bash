#!/bin/bash
#PBS -P ku3
#PBS -q express
#PBS -l other=gdata1a
#PBS -l walltime=18:00:00
#PBS -l mem=32GB
#PBS -l ncpus=1
## For licensed software, you have to specify it to get the job running. For unlicensed software, you should also specify it to help us analyse the software usage on our system.
#PBS -l software=cnvnator
## The job will be executed from current working directory instead of home.
#PBS -l wd 

# link dependencies:
export ROOTSYS=/home/913/jt3341/local/lib/root
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${ROOTSYS}/lib
export PATH=$PATH:$ROOTSYS/bin

# directory hierarchy:
homeDir="/g/data1a/ku3/jt3341/"
projectDir="$homeDir/projects/hgsoc_repeats/genome/"
resultsDir="$projectDir/results/"
mkdir -p $resultsDir

inDir="$resultsDir/bwa/"
outDir="$resultsDir/cnvnator/"
mkdir -p $outDir

scriptsDir="$projectDir/scripts"
genomeDir="$homeDir/genomes/GRCh37/"

inFile="$inDir/AOCS-076-1-3.bam"
sampleName=`basename $inFile | sed s/.bam//`

# extract read mapping from bam files:
#cnvnator -root $outDir/$sampleName.out.root -tree $inFile

# generate histogram:
#cnvnator -root $outDir/$sampleName.out.root -his 100

# calculate statistics:
#cnvnator -root $outDir/$sampleName.out.root -stat 100

# rd signal partitioning:
#cnvnator -root $outDir/$sampleName.out.root -partition 100

# call CNVs:
cnvnator -root $outDir/$sampleName.out.root -call 100 > $outDir/AOCS-076-1-3_CNV_calls.txt


