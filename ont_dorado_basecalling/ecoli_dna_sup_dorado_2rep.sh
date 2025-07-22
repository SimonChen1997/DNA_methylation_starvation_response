#!/bin/bash -l
#SBATCH --job-name="ecoli_dna_sup_dorado_2rep"
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=10G
#SBATCH --partition=gpu_cuda
#SBATCH --gres=gpu:h100:2
#SBATCH -o ecoli_dna_sup_dorado_2rep.o
#SBATCH -e ecoli_dna_sup_dorado_2rep.e

#########################################################

model=/path/dna_r10.4.1_e8.2_400bps_sup@v5.0.0

input_scratch=/path

output_scratch_bam=/path
output_scratch_demultiplex=/path

dorado=/path/dorado-0.8.0-linux-x64/bin/dorado

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
