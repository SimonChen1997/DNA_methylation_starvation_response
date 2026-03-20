#!/bin/bash -l
#SBATCH --job-name="ecoli_replication_start"
#SBATCH --partition=general
#SBATCH --array=1-5
#SBATCH -o ecoli_replication_start.o
#SBATCH -e ecoli_replication_start.e

#########################################################
module load anaconda3

#########################################################
### extract the mapping bam file
module load samtools
samtools view -F 0x900 -F 16 -F 4 -hSb $minimap2_primary_60x/ecoli_exp_${SLURM_ARRAY_TASK_ID}_60x_primary.bam > $forward_reverse_strand_bam/ecoli_exp_${SLURM_ARRAY_TASK_ID}_60x_plus.bam
samtools view -f 16 -hSb $minimap2_primary_60x/ecoli_exp_${SLURM_ARRAY_TASK_ID}_60x_primary.bam > $forward_reverse_strand_bam/ecoli_exp_${SLURM_ARRAY_TASK_ID}_60x_minus.bam

samtools view -F 0x900 -F 16 -F 4 -hSb $minimap2_primary_60x/ecoli_sta_${SLURM_ARRAY_TASK_ID}_60x_primary.bam > $forward_reverse_strand_bam/ecoli_sta_${SLURM_ARRAY_TASK_ID}_60x_plus.bam
samtools view -f 16 -hSb $minimap2_primary_60x/ecoli_sta_${SLURM_ARRAY_TASK_ID}_60x_primary.bam > $forward_reverse_strand_bam/ecoli_sta_${SLURM_ARRAY_TASK_ID}_60x_minus.bam

#########################################################
### convert the mapping bam to bed
module load bedtools
bedtools bamtobed -i $forward_reverse_strand_bam/ecoli_exp_${SLURM_ARRAY_TASK_ID}_60x_plus.bam > $forward_reverse_strand_bed/ecoli_exp_${SLURM_ARRAY_TASK_ID}_60x_plus.bed
bedtools bamtobed -i $forward_reverse_strand_bam/ecoli_exp_${SLURM_ARRAY_TASK_ID}_60x_minus.bam > $forward_reverse_strand_bed/ecoli_exp_${SLURM_ARRAY_TASK_ID}_60x_minus.bed

bedtools bamtobed -i $forward_reverse_strand_bam/ecoli_sta_${SLURM_ARRAY_TASK_ID}_60x_plus.bam > $forward_reverse_strand_bed/ecoli_sta_${SLURM_ARRAY_TASK_ID}_60x_plus.bed
bedtools bamtobed -i $forward_reverse_strand_bam/ecoli_sta_${SLURM_ARRAY_TASK_ID}_60x_minus.bam > $forward_reverse_strand_bed/ecoli_sta_${SLURM_ARRAY_TASK_ID}_60x_minus.bed

#########################################################################################
### subset the bed file to tsv
awk 'BEGIN{OFS=IFS="\t"}{if ($6=="+") print $1,$2+1,$3+1,$6; else if ($6=="-") print $1,$3+1,$2+1,$6}' \
    $forward_reverse_strand_bed/ecoli_exp_${SLURM_ARRAY_TASK_ID}_60x_plus.bed > $subset_strand_bed/ecoli_exp_${SLURM_ARRAY_TASK_ID}_60x_plus.tsv
awk 'BEGIN{OFS=IFS="\t"}{if ($6=="+") print $1,$2+1,$3+1,$6; else if ($6=="-") print $1,$3+1,$2+1,$6}' \
    $forward_reverse_strand_bed/ecoli_exp_${SLURM_ARRAY_TASK_ID}_60x_minus.bed > $subset_strand_bed/ecoli_exp_${SLURM_ARRAY_TASK_ID}_60x_minus.tsv

awk 'BEGIN{OFS=IFS="\t"}{if ($6=="+") print $1,$2+1,$3+1,$6; else if ($6=="-") print $1,$3+1,$2+1,$6}' \
    $forward_reverse_strand_bed/ecoli_sta_${SLURM_ARRAY_TASK_ID}_60x_plus.bed > $subset_strand_bed/ecoli_sta_${SLURM_ARRAY_TASK_ID}_60x_plus.tsv
awk 'BEGIN{OFS=IFS="\t"}{if ($6=="+") print $1,$2+1,$3+1,$6; else if ($6=="-") print $1,$3+1,$2+1,$6}' \
    $forward_reverse_strand_bed/ecoli_sta_${SLURM_ARRAY_TASK_ID}_60x_minus.bed > $subset_strand_bed/ecoli_sta_${SLURM_ARRAY_TASK_ID}_60x_minus.tsv
