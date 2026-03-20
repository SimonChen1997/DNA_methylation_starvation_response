#!/bin/bash -l
#SBATCH --job-name="ecoli_dmr_motif_focus"
#SBATCH --partition=general
#SBATCH -o ecoli_dmr_motif_focus.o
#SBATCH -e ecoli_dmr_motif_focus.e

#########################################################
module load anaconda3
module load samtools
module load htslib/1.15.1-gcc-11.3.0

#########################################################
### compress and index the pileup bed file

for file in $modkit_pileup_primary/*.bed;do bgzip -k $file;done
for file in $modkit_pileup_primary/*.bed.gz;do tabix -p bed $file;done

#########################################################
### dmr analysis m6a
source activate ont-modkit

modkit dmr pair \
    -a $modkit_pileup_primary/ecoli_exp_1_m6a_gatc.bed.gz \
    -a $modkit_pileup_primary/ecoli_exp_2_m6a_gatc.bed.gz \
    -a $modkit_pileup_primary/ecoli_exp_3_m6a_gatc.bed.gz \
    -b $modkit_pileup_primary/ecoli_sta_1_m6a_gatc.bed.gz \
    -b $modkit_pileup_primary/ecoli_sta_2_m6a_gatc.bed.gz \
    -b $modkit_pileup_primary/ecoli_sta_3_m6a_gatc.bed.gz \
    --min-valid-coverage 20 \
    --max-coverages 39 38 \
    --ref $ref \
    --out-path $modkit_dmr_file_m6a \
    --header \
    --base A \
    -t 4

### dmr analysis m4c
modkit dmr pair \
    -a $modkit_pileup_primary/ecoli_exp_1_m4c_garcntc.bed.gz \
    -a $modkit_pileup_primary/ecoli_exp_2_m4c_garcntc.bed.gz \
    -a $modkit_pileup_primary/ecoli_exp_3_m4c_garcntc.bed.gz \
    -b $modkit_pileup_primary/ecoli_sta_1_m4c_garcntc.bed.gz \
    -b $modkit_pileup_primary/ecoli_sta_2_m4c_garcntc.bed.gz \
    -b $modkit_pileup_primary/ecoli_sta_3_m4c_garcntc.bed.gz \
    --min-valid-coverage 20 \
    --max-coverages 40 39 \
    --ref $ref \
    --out-path $modkit_dmr_file_m4c \
    --header \
    --base C \
    -t 4

### dmr analysis m5c
modkit dmr pair \
    -a $modkit_pileup_primary/ecoli_exp_1_m5c_ccwgg.bed.gz \
    -a $modkit_pileup_primary/ecoli_exp_2_m5c_ccwgg.bed.gz \
    -a $modkit_pileup_primary/ecoli_exp_3_m5c_ccwgg.bed.gz \
    -b $modkit_pileup_primary/ecoli_sta_1_m5c_ccwgg.bed.gz \
    -b $modkit_pileup_primary/ecoli_sta_2_m5c_ccwgg.bed.gz \
    -b $modkit_pileup_primary/ecoli_sta_3_m5c_ccwgg.bed.gz \
    --min-valid-coverage 20 \
    --max-coverages 41 40 \
    --ref $ref \
    --out-path $modkit_dmr_file_m5c \
    --header \
    --base C \
    -t 4