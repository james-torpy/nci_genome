#!/bin/bash

script_dir="/g/data1a/ku3/jt3341/projects/hgsoc_repeats/genome/scripts"
mdss_dir="jt3341/projects/hgsoc_repeats/genome/raw_files/bowtell/"


for f in $(mdss ls $mdss_dir)
do
  
  echo "Getting $f..."
  echo -e

  qsub -v f=$f "$script_dir/mdss_get_bams.bash"

done