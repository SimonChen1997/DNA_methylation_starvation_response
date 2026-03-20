#!/bin/bash -l
#SBATCH --job-name="lacidophilus_sliding_window_bed"
#SBATCH --partition=general
#SBATCH -o lacidophilus_sliding_window_bed.o
#SBATCH -e lacidophilus_sliding_window_bed.e

#########################################################
module load bedtools/2.30.0-gcc-10.3.0

############################################################################
## extract the chromosome name and genome length columns
awk 'BEGIN{IFS=OFS="\t"}{print$1, $2}' $ref_dir/GCF_034298135.1_ASM3429813v1_genomic.fna.fai > $ref_dir/GCF_034298135.1_ASM3429813v1_refseq_chrids_size.tsv

## use bed tools to crate windows of each chromosome with 10kb sliding windows
bedtools makewindows -g $ref_dir/GCF_034298135.1_ASM3429813v1_refseq_chrids_size.tsv -w 10000 > $ref_dir/GCF_034298135.1_ASM3429813v1_windows_10kb.bed
