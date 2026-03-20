#!/bin/bash -l
#SBATCH --job-name="ecoli_seq_depth_position"
#SBATCH --partition=general
#SBATCH --array=1-5
#SBATCH -o ecoli_seq_depth_position.o
#SBATCH -e ecoli_seq_depth_position.e

##############################################################
module load samtools

##############################################################
### samtools to calculate the sequencing depth of each position
samtools depth --threads 5 $minimap2_primary_60x/ecoli_exp_${SLURM_ARRAY_TASK_ID}_60x_primary.sorted.bam > $samtools_depth/ecoli_exp_${SLURM_ARRAY_TASK_ID}_60x_primary.depth
samtools depth --threads 5 $minimap2_primary_60x/ecoli_sta_${SLURM_ARRAY_TASK_ID}_60x_primary.sorted.bam > $samtools_depth/ecoli_sta_${SLURM_ARRAY_TASK_ID}_60x_primary.depth
