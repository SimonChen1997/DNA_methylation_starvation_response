#!/bin/bash -l
#SBATCH --job-name="ecoli_dna_sup_dorado_2rep"
#SBATCH --partition=gpu_cuda
#SBATCH --gres=gpu:h100:2
#SBATCH -o ecoli_dna_sup_dorado_2rep.o
#SBATCH -e ecoli_dna_sup_dorado_2rep.e

#########################################################
module load cuda
module load anaconda3

#########################################################
### basecalling

$dorado basecaller --no-trim --recursive $model $input_scratch --modified-bases 6mA 4mC_5mC \
    --kit-name SQK-NBD114-96 \
    > $output_scratch_bam/ecoli_dna_exp_sta_6mA_4mC_5mC_sup_2rep.bam

$dorado demux --output-dir $output_scratch_demultiplex \
    --kit-name SQK-NBD114-96 $output_scratch_bam/ecoli_dna_exp_sta_6mA_4mC_5mC_sup_2rep.bam
