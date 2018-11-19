#!/bin/bash 
#PBS -P ku3
#PBS -q express
#PBS -l ncpus=16
#PBS -l mem=128GB
#PBS -l walltime=24:00:00
#PBS -l other=snakemake
#PBS -l other=gdata1a
#PBS -l wd

project_dir="/g/data1a/ku3/jt3341/projects/hgsoc_repeats/genome"
log_dir="$project_dir/logs/"

conda activate p3.6.3env
cd $project_dir
#snakemake -s $project_dir/retro_SV_enrichment_Snakefile --cores 16 -k
snakemake -s retro_SV_enrichment_Snakefile_subs -k \
	-j 16 --cluster "qsub -P ku3 -q express -l ncpus=16 -l mem=32GB \
	-l walltime=00:05:00"