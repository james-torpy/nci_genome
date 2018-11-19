#!/bin/bash 
#PBS -P ku3
#PBS -q copyq
#PBS -l ncpus=1
#PBS -l mem=2GB
#PBS -l walltime=10:00:00
#PBS -l other=mdss
#PBS -l other=gdata1a
#PBS -l wd

project_dir="/g/data1a/ku3/jt3341/projects"
md5_dir="$project_dir/hgsoc_repeats/genome/raw_files"

echo "md5 checksumming ${f}"
echo -e

md5sum ${f} >> $md5_dir/AOCS_genome_bams.md5

echo "Checksum complete! Gzipping ${f}"
echo -e

gzip ${f}

