#!/bin/bash 
#PBS -P ku3
#PBS -q copyq
#PBS -l ncpus=1
#PBS -l mem=2GB
#PBS -l walltime=10:00:00
#PBS -l other=mdss
#PBS -l other=gdata1a
#PBS -l wd

mdss_dir="jt3341/projects/hgsoc_repeats/genome/raw_files/bowtell/"

project_dir="/g/data1a/ku3/jt3341/projects/"
md5_dir="$project_dir/hgsoc_repeats/genome/raw_files/bowtell/AOCS_genome_bams/"
out_dir="$md5_dir/mdss/"
mkdir -p $out_dir

mdss_dir="jt3341/projects/hgsoc_repeats/genome/raw_files/bowtell/"

echo "Getting ${f} from $mdss_dir"
echo -e

mdss get $mdss_dir/${f} $out_dir

if [ -f $out_dir/${f} ]; then
	
	echo "${f} got!"
	echo -e

	echo "md5 checksumming"
	echo -e

	md5sum $out_dir/${f} >> $md5_dir/AOCS_genome_bams.md5

	echo "Checksum complete!"
	echo -e

	touch "$out_path/${f}_got_and_checksummed"

fi

