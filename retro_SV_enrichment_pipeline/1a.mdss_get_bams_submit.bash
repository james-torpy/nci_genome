#!/bin/bash

script_dir="/g/data1a/ku3/jt3341/projects/hgsoc_repeats/genome/scripts/retro_SV_enrichment_pipeline/"
mdss_dir="jt3341/projects/hgsoc_repeats/genome/raw_files/"

fils=( "AOCS-001-1.bam.gz" "AOCS-001-5.bam.gz" ) #\
	#"AOCS-005sub-1.bam.gz" "AOCS-005sub-5.bam.gz" )
#fils=( "AOCS-034-1.bam.gz" "AOCS-034-5.bam.gz" \
#	"AOCS-055-1.bam.gz" "AOCS-055-5.bam.gz" )
#fils=( $(mdss ls $mdss_dir/*-[1,5]-* | head -20) )

for f in ${fils[@]}
do
  
  echo "Getting $f..."
  echo -e

  qsub -v f=$f "$script_dir/1b.mdss_get_bams.bash"

done