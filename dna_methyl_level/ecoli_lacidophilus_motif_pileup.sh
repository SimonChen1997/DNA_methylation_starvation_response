#!/bin/bash -l
#SBATCH --job-name="ecoli_lacidophilus_motif_pileup"
#SBATCH --partition=general
#SBATCH --array=1-3
#SBATCH -o ecoli_lacidophilus_motif_pileup.o
#SBATCH -e ecoli_lacidophilus_motif_pileup.e

#########################################################
module load anaconda3
module load samtools

#########################################################
### get modification at base level
source activate ont-modkit

########## e.coli
## targeting m6a
modkit pileup $minimap2_primary_m6a_60x/ecoli_exp_${SLURM_ARRAY_TASK_ID}_primary_m6a.sorted.bam $ecoli_motif_focus_pileup/ecoli_exp_${SLURM_ARRAY_TASK_ID}_m6a_gatc.bed --motif GATC 1 --ref $ecoli_ref
modkit pileup $minimap2_primary_m6a_60x/ecoli_sta_${SLURM_ARRAY_TASK_ID}_primary_m6a.sorted.bam $ecoli_motif_focus_pileup/ecoli_sta_${SLURM_ARRAY_TASK_ID}_m6a_gatc.bed --motif GATC 1 --ref $ecoli_ref

## targeting m4c
modkit pileup $minimap2_primary_m4c_60x/ecoli_exp_${SLURM_ARRAY_TASK_ID}_primary_m4c.sorted.bam $ecoli_motif_focus_pileup/ecoli_exp_${SLURM_ARRAY_TASK_ID}_m4c_garcntc.bed --motif GATCNTC 6 --ref $ecoli_ref
modkit pileup $minimap2_primary_m4c_60x/ecoli_sta_${SLURM_ARRAY_TASK_ID}_primary_m4c.sorted.bam $ecoli_motif_focus_pileup/ecoli_sta_${SLURM_ARRAY_TASK_ID}_m4c_garcntc.bed --motif GATCNTC 6 --ref $ecoli_ref

## targeting m5c
modkit pileup $minimap2_primary_m5c_60x/ecoli_exp_${SLURM_ARRAY_TASK_ID}_primary_m5c.sorted.bam $ecoli_motif_focus_pileup/ecoli_exp_${SLURM_ARRAY_TASK_ID}_m5c_ccwgg.bed --motif CCWGG 1 --ref $ecoli_ref
modkit pileup $minimap2_primary_m5c_60x/ecoli_sta_${SLURM_ARRAY_TASK_ID}_primary_m5c.sorted.bam $ecoli_motif_focus_pileup/ecoli_sta_${SLURM_ARRAY_TASK_ID}_m5c_ccwgg.bed --motif CCWGG 1 --ref $ecoli_ref


########## l.acidophilus
## targeting m6a
modkit pileup $minimap2_primary_m6a_100x/lacidophilus_exp_${SLURM_ARRAY_TASK_ID}_primary_m6a.sorted.bam \
    $lacidophilus_motif_focus_pileup/lacidophilus_exp_${SLURM_ARRAY_TASK_ID}_m6a_gcatc.bed --motif GCATC 2 --ref $lacidophilus_ref

modkit pileup $minimap2_primary_m6a_100x/lacidophilus_sta_${SLURM_ARRAY_TASK_ID}_primary_m6a.sorted.bam \
    $lacidophilus_motif_focus_pileup/lacidophilus_sta_${SLURM_ARRAY_TASK_ID}_m6a_gcatc.bed --motif GCATC 2 --ref $lacidophilus_ref