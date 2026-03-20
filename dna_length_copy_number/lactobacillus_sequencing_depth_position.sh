#!/bin/bash -l
#SBATCH --job-name="lactobacillus_seq_depth_position"
#SBATCH --partition=general
#SBATCH --array=1-3
#SBATCH -o lactobacillus_seq_depth_position.o
#SBATCH -e lactobacillus_seq_depth_position.e

##############################################################
module load samtools

##############################################################
### samtools to calculate the sequencing depth of each position
samtools depth --threads 5 $minimap2_primary_fastq_100x/lacidophilus_exp_${SLURM_ARRAY_TASK_ID}_100x_primary.sorted.bam > $samtools_depth/lacidophilus_exp_${SLURM_ARRAY_TASK_ID}_100x_primary.depth
samtools depth --threads 5 $minimap2_primary_fastq_100x/lacidophilus_sta_${SLURM_ARRAY_TASK_ID}_100x_primary.sorted.bam > $samtools_depth/lacidophilus_sta_${SLURM_ARRAY_TASK_ID}_100x_primary.depth
