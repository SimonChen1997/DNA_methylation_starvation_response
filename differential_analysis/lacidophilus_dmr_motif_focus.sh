#!/bin/bash -l
#SBATCH --job-name="lacidophilus_dmr_motif_focus"
#SBATCH --partition=general
#SBATCH -o lacidophilus_dmr_motif_focus.o
#SBATCH -e lacidophilus_dmr_motif_focus.e

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
    -a $modkit_pileup_primary/lacidophilus_exp_1_m6a_gcatc.bed.gz \
    -a $modkit_pileup_primary/lacidophilus_exp_2_m6a_gcatc.bed.gz \
    -a $modkit_pileup_primary/lacidophilus_exp_3_m6a_gcatc.bed.gz \
    -b $modkit_pileup_primary/lacidophilus_sta_1_m6a_gcatc.bed.gz \
    -b $modkit_pileup_primary/lacidophilus_sta_2_m6a_gcatc.bed.gz \
    -b $modkit_pileup_primary/lacidophilus_sta_3_m6a_gcatc.bed.gz \
    --min-valid-coverage 20 \
    --max-coverages 46 54 \
    --ref $ref \
    --out-path $modkit_dmr_file_m6a \
    --header \
    --base A \
    -t 4