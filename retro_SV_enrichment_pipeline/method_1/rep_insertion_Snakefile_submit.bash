#!/bin/bash
#PBS -P ku3
#PBS -q normal
#PBS -l other=gdata1a
#PBS -l walltime=24:00:00
#PBS -l mem=128GB
#PBS -l ncpus=16
## For licensed software, you have to specify it to get the job running. For unlicensed software, you should also specify it to help us analyse the software usage on our system.
#PBS -l software=snakemake
## The job will be executed from current working directory instead of home.
#PBS -l wd

# activate snakemake environment
conda activate p3.6.3env

# link dependencies:
export ROOTSYS=/home/913/jt3341/local/lib/root
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${ROOTSYS}/lib
export PATH=$PATH:$ROOTSYS/bin

# submit snakemake job:
snakemake -s /g/data1a/ku3/jt3341/projects/hgsoc_repeats/genome/rep_insertion_Snakefile --unlock
snakemake -s /g/data1a/ku3/jt3341/projects/hgsoc_repeats/genome/rep_insertion_Snakefile -k --cores=15

# for cnvnator, 1 NCPU, 8 GB per job
