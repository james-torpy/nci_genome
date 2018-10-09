#!/bin/bash
project_dir="/g/data1a/ku3/jt3341/projects/hgsoc_repeats/genome"
script_dir="$project_dir/scripts"
bam_dir="$project_dir/raw_files/bowtell/AOCS_genome_bams/mdss"


for f in $bam_dir/*.bam
do
  
  echo "Gzipping $f..."
  echo -e

  qsub -v f=$f "$script_dir/gzip_bams.bash"

done