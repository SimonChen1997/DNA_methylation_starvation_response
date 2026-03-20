#!/bin/bash -l
#SBATCH --job-name="ecoli_rna_seq"
#SBATCH --partition=general
#SBATCH --array=1-3
#SBATCH -o ecoli_rna_seq.o
#SBATCH -e ecoli_rna_seq.e

##############################################################
module load anaconda3

##############################################################
### fastqc for quality control check
source activate fastqc

for i in {1,2};do
	fastqc $fastq/EE_${SLURM_ARRAY_TASK_ID}_R${i}.fastq -q -o $fastqc
	fastqc $fastq/ES_${SLURM_ARRAY_TASK_ID}_R${i}.fastq -q -o $fastqc
done

##############################################################
### fastp for adapter cutting
source activate fastp

fastp --thread 8 -i $fastq/EE_${SLURM_ARRAY_TASK_ID}_R1.fastq -I $fastq/EE_${SLURM_ARRAY_TASK_ID}_R2.fastq -o $fastp/EE_${SLURM_ARRAY_TASK_ID}_R1.fastq -O $fastp/EE_${SLURM_ARRAY_TASK_ID}_R2.fastq
fastp --thread 8 -i $fastq/ES_${SLURM_ARRAY_TASK_ID}_R1.fastq -I $fastq/ES_${SLURM_ARRAY_TASK_ID}_R2.fastq -o $fastp/ES_${SLURM_ARRAY_TASK_ID}_R1.fastq -O $fastp/ES_${SLURM_ARRAY_TASK_ID}_R2.fastq

##############################################################
### use STAR to build genome index
source activate star

STAR --runThreadN 2 --runMode genomeGenerate --genomeSAindexNbases 10 --genomeDir $ref_index --genomeFastaFiles $ref --sjdbGTFfile $gtf_ref --sjdbOverhang 150-1 

##############################################################
### use STAR to for alignment
source activate star

for i in {1,2};do
	STAR --runThreadN 8 --genomeDir $ref_index --readFilesIn $fastp/EE_${SLURM_ARRAY_TASK_ID}_R1.fastq $fastp/EE_${SLURM_ARRAY_TASK_ID}_R2.fastq --outFileNamePrefix $star/EE_${SLURM_ARRAY_TASK_ID}_
	STAR --runThreadN 8 --genomeDir $ref_index --readFilesIn $fastp/ES_${SLURM_ARRAY_TASK_ID}_R1.fastq $fastp/ES_${SLURM_ARRAY_TASK_ID}_R2.fastq --outFileNamePrefix $star/ES_${SLURM_ARRAY_TASK_ID}_
done

##############################################################
### htseq for raw counts
source activate htseq

htseq-count -n 8 -f sam -s yes -t CDS -m intersection-strict $star/EE_${SLURM_ARRAY_TASK_ID}_Aligned.out.sam $gtf_ref > $htseq/EE_${SLURM_ARRAY_TASK_ID}_gene_counts.tsv
htseq-count -n 8 -f sam -s yes -t CDS -m intersection-strict $star/ES_${SLURM_ARRAY_TASK_ID}_Aligned.out.sam $gtf_ref > $htseq/ES_${SLURM_ARRAY_TASK_ID}_gene_counts.tsv