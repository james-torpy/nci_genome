#!/bin/bash 
#PBS -P ku3
#PBS -q copyq
#PBS -l ncpus=1
#PBS -l mem=2GB
#PBS -l walltime=10:00:00
#PBS -l other=samtools
#PBS -l other=gdata1a
#PBS -l wd

mdss_dir="jt3341/projects/hgsoc_repeats/genome/raw_files/"

project_dir="/g/data1a/ku3/jt3341/projects//hgsoc_repeats/genome/"
script_dir="$project_dir/scripts/retro_SV_enrichment_pipeline/"
bam_dir="$project_dir/raw_files/bams/"
check_dir="$bam_dir/check/"
mkdir -p $check_dir

echo "Getting ${f} from $mdss_dir"
echo -e

# define bam name:
bam=$(echo ${f} | sed "s/\.gz//")

# if bam not in bam_dir, get bam from MDSS:
if [ ! -f $bam_dir/${f} ]; then
	mdss get $mdss_dir/${f} $bam_dir
fi

if [ -f $bam_dir/${f} ]; then
	echo "${f} got!"
	echo -e
	touch $check_dir/${f}.got
fi

# if all 4 files have been transferred, process using snakemake pipeline:
wait

# if all 4 files have been transferred, process using snakemake pipeline:
fileno=$(ls $check_dir/*got | wc -l)
jorbz=( $(qstat -s | grep jt3341 | grep 2a.retro ) )
if [ "$fileno" = 4 ] && [ -z ${jorbz+x} ]; then
	echo "Initiating retrotransposon SV enrichment snakemake pipeline..."
	echo -e
	conda activate p3.6.3env
	cd $project_dir
	#snakemake -s $project_dir/retro_SV_enrichment_Snakefile --cores 16 -k
	snakemake -s retro_SV_enrichment_Snakefile -k --cluster "qsub -P ku3 -q \
	express -l ncpus=16 -l mem=128GB -l walltime=24:00:00" -j 16
fi




