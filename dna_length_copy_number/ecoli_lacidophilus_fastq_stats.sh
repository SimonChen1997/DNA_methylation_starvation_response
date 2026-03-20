#!/bin/bash -l
#SBATCH --job-name="ecoli_lacidophilus_fastq_stats"
#SBATCH --partition=general
#SBATCH -o ecoli_lacidophilus_fastq_stats.o
#SBATCH -e ecoli_lacidophilus_fastq_stats.e
#SBATCH --account=a_eross

##############################################################
module load anaconda3

##############################################################
### ecoli
source activate seqkit

echo -e "sample_id\taverage_length\tn50\tgenus\tphase\treplicate_id" > $ecoli_fastq_stats/ecoli_lacidophilus_primary_fastq_stat.tsv

for file in $ecoli_fastq_60x/*.fastq;do
	name_1=${file%%.*}
	name_2=${name_1//_60x}
	name_3=${name_2##*/}
	genus="Escherichia"
	
	if [[ ${name_3} == *exp* ]]; then
        phase="exponential"
    elif [[ ${name_3} == *sta* ]]; then
        phase="stationary"
    fi
    
    if [[ ${name_3} == *_1* ]]; then
        replicate="replicate_1"
    elif [[ ${name_3} == *_2* ]]; then
        replicate="replicate_2"
    elif [[ ${name_3} == *_3* ]]; then
        replicate="replicate_3"
    elif [[ ${name_3} == *_4* ]]; then
        replicate="replicate_4"
    elif [[ ${name_3} == *_5* ]]; then
        replicate="replicate_5"
    fi
    
	length=$(awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}' $file)
	n50_value=$(seqkit stats $file -a | awk 'NR==2 {print $13}')
	echo -e "$name_3\t$length\t$n50_value\t$genus\t$phase\t$replicate" >> $ecoli_fastq_stats/ecoli_lacidophilus_primary_fastq_stat.tsv
done

### lacidophilus
for file in $lacidophilus_fastq_100x/*.fastq;do
	name_1=${file%%.*}
	name_2=${name_1//_100x}
	name_3=${name_2##*/}
	genus="Lactobacillus"
	
	if [[ ${name_3} == *exp* ]]; then
        phase="exponential"
    elif [[ ${name_3} == *sta* ]]; then
        phase="stationary"
    fi
    
    if [[ ${name_3} == *_1* ]]; then
        replicate="replicate_1"
    elif [[ ${name_3} == *_2* ]]; then
        replicate="replicate_2"
    elif [[ ${name_3} == *_3* ]]; then
        replicate="replicate_3"
    fi
    
	length=$(awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}' $file)
	n50_value=$(seqkit stats $file -a | awk 'NR==2 {print $13}')
	echo -e "$name_3\t$length\t$n50_value\t$genus\t$phase\t$replicate" >> $ecoli_fastq_stats/ecoli_lacidophilus_primary_fastq_stat.tsv
done

