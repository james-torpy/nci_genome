#!/bin/bash

sample="AOCS-152-1-X"

# transfer sample from mdss:
cd /g/data1a/ku3/jt3341/projects/hgsoc_repeats/genome/results/bwa/
#mdss get jt3341/projects/hgsoc_repeats/genome/raw_files/bowtell/$sample.bam

# run snakemake submit script:
qsub /g/data1a/ku3/jt3341/projects/hgsoc_repeats/genome/scripts/retro_insertions/rep_insertion_Snakemake_submit.bash
