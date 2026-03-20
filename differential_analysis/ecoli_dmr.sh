#!/bin/bash -l
#SBATCH --job-name="ecoli_dmr"
#SBATCH --partition=general
#SBATCH -o ecoli_dmr.o
#SBATCH -e ecoli_dmr.e

#########################################################
module load anaconda3
module load samtools
module load htslib/1.15.1-gcc-11.3.0

#########################################################
### dmr analysis m6a
source activate ont-modkit

modkit dmr pair \
    -a $modkit_pileup_primary/ecoli_exp_1_pileup_m6a.bed.gz \
    -a $modkit_pileup_primary/ecoli_exp_2_pileup_m6a.bed.gz \
    -a $modkit_pileup_primary/ecoli_exp_3_pileup_m6a.bed.gz \
    -b $modkit_pileup_primary/ecoli_sta_1_pileup_m6a.bed.gz \
    -b $modkit_pileup_primary/ecoli_sta_2_pileup_m6a.bed.gz \
    -b $modkit_pileup_primary/ecoli_sta_3_pileup_m6a.bed.gz \
    --min-valid-coverage 20 \
    --ref $ref \
    --out-path $modkit_dmr_file_m6a \
    --header \
    --base A \
    -t 4

### dmr analysis m4c
modkit dmr pair \
    -a $modkit_pileup_primary/ecoli_exp_1_pileup_m4c.bed.gz \
    -a $modkit_pileup_primary/ecoli_exp_2_pileup_m4c.bed.gz \
    -a $modkit_pileup_primary/ecoli_exp_3_pileup_m4c.bed.gz \
    -b $modkit_pileup_primary/ecoli_sta_1_pileup_m4c.bed.gz \
    -b $modkit_pileup_primary/ecoli_sta_2_pileup_m4c.bed.gz \
    -b $modkit_pileup_primary/ecoli_sta_3_pileup_m4c.bed.gz \
    --min-valid-coverage 20 \
    --ref $ref \
    --out-path $modkit_dmr_file_m4c \
    --header \
    --base C \
    -t 4

### dmr analysis m5c
modkit dmr pair \
    -a $modkit_pileup_primary/ecoli_exp_1_pileup_m5c.bed.gz \
    -a $modkit_pileup_primary/ecoli_exp_2_pileup_m5c.bed.gz \
    -a $modkit_pileup_primary/ecoli_exp_3_pileup_m5c.bed.gz \
    -b $modkit_pileup_primary/ecoli_sta_1_pileup_m5c.bed.gz \
    -b $modkit_pileup_primary/ecoli_sta_2_pileup_m5c.bed.gz \
    -b $modkit_pileup_primary/ecoli_sta_3_pileup_m5c.bed.gz \
    --min-valid-coverage 20 \
    --ref $ref \
    --out-path $modkit_dmr_file_m5c \
    --header \
    --base C \
    -t 4
