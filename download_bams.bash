#!/bin/bash 
#PBS -P ku3
#PBS -q copyq
#PBS -l ncpus=1
#PBS -l mem=2GB
#PBS -l walltime=10:00:00
#PBS -l other=mdss
#PBS -l other=gdata1a
#PBS -l wd

# change download links and num for each batch:
script_dir="/g/data1a/ku3/jt3341/projects/hgsoc_repeats/genome/scripts"

out_path="/g/data1a/ku3/jt3341/projects/hgsoc_repeats/genome/raw_files/bowtell"
out_dir="$out_path/AOCS_genome_bams"
mkdir -p $out_dir

mdss_dir="jt3341/projects/hgsoc_repeats/genome/raw_files/bowtell/"

echo "Downloading ${id} via ${lnk}"
echo -e

curl -d from-share=true -o $out_dir/${id} ${lnk}

if [ -f $out_dir/${id} ]; then
	
	echo "${id} downloaded!"
	echo -e

	echo "md5 checksumming"
	echo -e

	md5sum $out_dir/${id} >> $out_dir/AOCS_genome_bams.md5

	echo "Checksum complete!"
	echo -e

	touch "$out_path/${id}_downloaded_and_checksummed"

fi

