#!/bin/bash 
#PBS -P ku3
#PBS -q express
#PBS -l ncpus=1
#PBS -l mem=2GB
#PBS -l walltime=24:00:00
#PBS -l other=snakemake
#PBS -l other=gdata1a
#PBS -l wd

project_dir="/g/data1a/ku3/jt3341/projects/hgsoc_repeats/genome/"

echo "Initiating retrotransposon SV enrichment snakemake pipeline..."
echo -e
conda activate p3.6.3env
cd $project_dir
snakemake -s retro_SV_enrichment_Snakefile -k --cluster "qsub -P ku3 -q \
	express -l ncpus=8 -l mem=64GB -l walltime=24:00:00"  -j 8 --jobs 4