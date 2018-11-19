#!/bin/bash 
#PBS -P ku3
#PBS -q express
#PBS -l ncpus=8
#PBS -l mem=64GB
#PBS -l walltime=24:00:00
#PBS -l other=bwa
#PBS -l other=gdata1a
#PBS -l wd

bwa index -p /g/data1a/ku3/jt3341/projects/hgsoc_repeats/genome/refs/hg38_ercc/hg38_ercc.bwa \
/g/data1a/ku3/jt3341/projects/hgsoc_repeats/genome/refs/hg38_ercc/hg38_ercc.fa

