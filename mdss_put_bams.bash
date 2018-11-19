#!/bin/bash 
#PBS -P ku3
#PBS -q copyq
#PBS -l ncpus=1
#PBS -l mem=2GB
#PBS -l walltime=10:00:00
#PBS -l other=mdss
#PBS -l other=gdata1a
#PBS -l wd

mdss_dir="jt3341/projects/hgsoc_repeats/genome/raw_files/"

echo "Putting ${f} into $mdss_dir"
echo -e

mdss put ${f} $mdss_dir

id=$(echo ${f} | sed "s/^.*AOCS/AOCS/")

# check if file transferred, and if so delete:
chek=""
chek=$(mdss ls $mdss_dir/$id)
if [ ${#chek} != 0 ]; then
	rm ${f}
fi
