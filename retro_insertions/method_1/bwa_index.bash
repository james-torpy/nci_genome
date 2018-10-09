#!/bin/bash

#!/bin/bash
#PBS -P ku3
#PBS -q express
#PBS -l walltime=10:00:00
#PBS -l mem=64GB
#PBS -l ncpus=1
## For licensed software, you have to specify it to get the job running. For unlicensed software, you should also specify it to help us analyse the software usage on our system.
#PBS -l software=bwa
## The job will be executed from current working directory instead of home.
#PBS -l wd

# took 01:06:37 with 1 CPU, 128 GB #

module load BWA/0.7.17
module load gcc/4.9.0
module load samtools/0.1.18

homeDir="/g/data1a/ku3/jt3341/"
genomeDir="$homeDir/genomes/GRCh37/"

bwa index $genomeDir/GRCh37.fa

