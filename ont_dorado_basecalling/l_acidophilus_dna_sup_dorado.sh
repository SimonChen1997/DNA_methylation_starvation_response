#!/bin/bash -l
#SBATCH --job-name="l_acidophilus_dna_sup_dorado"
#SBATCH --qos=gpu
#SBATCH --partition=gpu_cuda
#SBATCH --gres=gpu:h100:2
#SBATCH -o l_acidophilus_dna_sup_dorado.o
#SBATCH -e l_acidophilus_dna_sup_dorado.e

#########################################################
module load cuda
module load anaconda3

#########################################################
### basecalling

$dorado basecaller --no-trim --recursive $model $input_scratch_pod5 --modified-bases 6mA 4mC_5mC \
    --kit-name SQK-NBD114-96 \
    > $output_scratch_bam/lacidophilus_3_rep_methyl_sup.bam

$dorado demux --output-dir $output_scratch_demultiplex \
    --kit-name SQK-NBD114-96 $output_scratch_bam/lacidophilus_3_rep_methyl_sup.bam