#!/bin/bash -l
#SBATCH --job-name="l_acidophilus_dna_sup_dorado"
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=10G
#SBATCH --qos=gpu
#SBATCH --partition=gpu_cuda
#SBATCH --gres=gpu:h100:2
#SBATCH -o l_acidophilus_dna_sup_dorado.o
#SBATCH -e l_acidophilus_dna_sup_dorado.e

#########################################################

model=/path/dna_r10.4.1_e8.2_400bps_sup@v5.0.0

input_scratch_pod5=/path

output_scratch_bam=/path
output_scratch_demultiplex=/path

dorado=/scratch/project/genoepic_rumen/dorado-0.8.0-linux-x64/bin/dorado

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