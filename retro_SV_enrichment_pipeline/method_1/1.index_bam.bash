#!/bin/bash
#PBS -l ncpus=1
#PBS -l mem=64GB
#PBS -l jobfs=30GB
#PBS -l walltime=20:00:00
#PBS -P ku3
#PBS -l wd
#PBS -l other=gdata1a
#PBS -q normal
#PBS -N samtools
set -euo pipefail
module load samtools/1.9
samtools index /g/data1a/ku3/jt3341/projects/hgsoc_repeats/genome/results/bwa/AOCS-076-1-3_chr1.bam
