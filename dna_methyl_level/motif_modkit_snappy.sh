#!/bin/bash -l
#SBATCH --job-name="motif_modkit_snappy"
#SBATCH --partition=general
#SBATCH --array=xxx
#SBATCH -o motif_modkit_snappy.o
#SBATCH -e motif_modkit_snappy.e

#########################################################
module load anaconda3/2023.09-0

#########################################################
### get modification at base level
source activate ont-modkit

#### strain
## targeting m6a
modkit motif search -i $modkit_bed/barcode${SLURM_ARRAY_TASK_ID}_m6a.bed -r $ref \
    --min-coverage 20 --low-thresh 0.01 --min-sites 110 --high-thresh 0.7 --min-frac-mod 0.7 --skip-search --min-log-odds 1.6 \
    -o $motif_modkit/barcode${SLURM_ARRAY_TASK_ID}_m6a.tsv --threads 5 --log $motif_modkit/barcode${SLURM_ARRAY_TASK_ID}_m6a_log.txt --force-override-spec

## targeting m4C
modkit motif search -i $modkit_bed/barcode${SLURM_ARRAY_TASK_ID}_m4c.bed -r $ref \
    --min-coverage 20 --low-thresh 0.01 --min-sites 110 --high-thresh 0.7 --min-frac-mod 0.7 --skip-search --min-log-odds 1.6 \
    -o $motif_modkit/barcode${SLURM_ARRAY_TASK_ID}_m4c.tsv --threads 5 --log $motif_modkit/barcode${SLURM_ARRAY_TASK_ID}_m4c_log.txt --force-override-spec

## targeting m5c
modkit motif search -i $modkit_bed/barcode${SLURM_ARRAY_TASK_ID}_m5c.bed -r $ref \
    --min-coverage 20 --low-thresh 0.01 --min-sites 110 --high-thresh 0.7 --min-frac-mod 0.7 --skip-search --min-log-odds 1.6 \
    -o $motif_modkit/barcode${SLURM_ARRAY_TASK_ID}_m5c.tsv --threads 5 --log $motif_modkit/barcode${SLURM_ARRAY_TASK_ID}_m5c_log.txt --force-override-spec

#########################################################
source activate snappy
#### strain
## targeting m6a
snappy -mk_bed $modkit_bed/barcode${SLURM_ARRAY_TASK_ID}_m6a.bed -genome $ref -outdir $motif_snappy/barcode${SLURM_ARRAY_TASK_ID}_m6a

## targeting m4c
snappy -mk_bed $modkit_bed/barcode${SLURM_ARRAY_TASK_ID}_m4c.bed -genome $ref -outdir $motif_snappy/barcode${SLURM_ARRAY_TASK_ID}_m4c

## targeting m5c
snappy -mk_bed $modkit_bed/barcode${SLURM_ARRAY_TASK_ID}_m5c.bed -genome $ref -outdir $motif_snappy/barcode${SLURM_ARRAY_TASK_ID}_m5c