#!/bin/bash -l
#SBATCH --job-name="ecoli_seq_depth_determination"
#SBATCH --partition=general
#SBATCH --array=2,4,5
#SBATCH -o ecoli_seq_depth_determination.o
#SBATCH -e ecoli_seq_depth_determination.e

#########################################################
module load anaconda3
module load samtools

#########################################################
### extract fastq with methylation tag
samtools fastq -T ML,MM,MN $output_scratch_demultiplex/ecoli_exp_${SLURM_ARRAY_TASK_ID}.bam > $fastq_sup/ecoli_exp_${SLURM_ARRAY_TASK_ID}.fastq
samtools fastq -T ML,MM,MN $output_scratch_demultiplex/ecoli_sta_${SLURM_ARRAY_TASK_ID}.bam > $fastq_sup/ecoli_sta_${SLURM_ARRAY_TASK_ID}.fastq

#########################################################
### first round mapping to e.coli reference genome
source activate minimap2

minimap2 -y -ax lr:hq $ref $fastq_sup/ecoli_exp_${SLURM_ARRAY_TASK_ID}.fastq -t 5 -o $minimap2/ecoli_exp_${SLURM_ARRAY_TASK_ID}_align.bam
minimap2 -y -ax lr:hq $ref $fastq_sup/ecoli_sta_${SLURM_ARRAY_TASK_ID}.fastq -t 5 -o $minimap2/ecoli_sta_${SLURM_ARRAY_TASK_ID}_align.bam

### index the sam file and convert to bam file
samtools view -bhS $minimap2/ecoli_exp_${SLURM_ARRAY_TASK_ID}_align.bam | samtools sort -T $minimap2/ecoli_exp_${SLURM_ARRAY_TASK_ID}_align.sorted -o $minimap2/ecoli_exp_${SLURM_ARRAY_TASK_ID}_align.sorted.bam
samtools index $minimap2/ecoli_exp_${SLURM_ARRAY_TASK_ID}_align.sorted.bam

samtools view -bhS $minimap2/ecoli_sta_${SLURM_ARRAY_TASK_ID}_align.bam | samtools sort -T $minimap2/ecoli_sta_${SLURM_ARRAY_TASK_ID}_align.sorted -o $minimap2/ecoli_sta_${SLURM_ARRAY_TASK_ID}_align.sorted.bam
samtools index $minimap2/ecoli_sta_${SLURM_ARRAY_TASK_ID}_align.sorted.bam

#########################################################
### extract only the primary mapped reads
samtools view -F 0x900 -F 4 -bhS $minimap2/ecoli_exp_${SLURM_ARRAY_TASK_ID}_align.bam > $minimap2_primary/ecoli_exp_${SLURM_ARRAY_TASK_ID}_primary.bam
samtools view -F 0x900 -F 4 -bhS $minimap2/ecoli_sta_${SLURM_ARRAY_TASK_ID}_align.bam > $minimap2_primary/ecoli_sta_${SLURM_ARRAY_TASK_ID}_primary.bam

#########################################################
### extract fastq with methylation tag from primary bam
samtools fastq -T ML,MM,MN $minimap2_primary/ecoli_exp_${SLURM_ARRAY_TASK_ID}_primary.bam > $fastq_primary/ecoli_exp_${SLURM_ARRAY_TASK_ID}_primary.fastq
samtools fastq -T ML,MM,MN $minimap2_primary/ecoli_sta_${SLURM_ARRAY_TASK_ID}_primary.bam > $fastq_primary/ecoli_sta_${SLURM_ARRAY_TASK_ID}_primary.fastq

echo "ecoli_exp_${SLURM_ARRAY_TASK_ID}_primary.fastq is done"
echo "ecoli_sta_${SLURM_ARRAY_TASK_ID}_primary.fastq is done"

#########################################################
### remove short reads
source activate nanopack

NanoFilt --length 250 $fastq_primary/ecoli_exp_${SLURM_ARRAY_TASK_ID}_primary.fastq > $nanofilt/ecoli_exp_${SLURM_ARRAY_TASK_ID}_filt.fastq
NanoFilt --length 250 $fastq_primary/ecoli_sta_${SLURM_ARRAY_TASK_ID}_primary.fastq > $nanofilt/ecoli_sta_${SLURM_ARRAY_TASK_ID}_filt.fastq

echo "ecoli_exp_${SLURM_ARRAY_TASK_ID}_filt.fastq is done"
echo "ecoli_sta_${SLURM_ARRAY_TASK_ID}_filt.fastq is done"

#########################################################
### subsample fastq to different coverage
source activate rasusa

for i in $(seq 5.5 5.5 550); do
	coverage=$(awk "BEGIN {print $i / 5.5}")
	rasusa reads -s 1727 --bases ${i}MB $nanofilt/ecoli_exp_${SLURM_ARRAY_TASK_ID}_filt.fastq -o $fastq_subsample/ecoli_exp_${SLURM_ARRAY_TASK_ID}_${coverage}x.fastq
	rasusa reads -s 1727 --bases ${i}MB $nanofilt/ecoli_sta_${SLURM_ARRAY_TASK_ID}_filt.fastq -o $fastq_subsample/ecoli_sta_${SLURM_ARRAY_TASK_ID}_${coverage}x.fastq
done

#########################################################
### map to e.coli reference genome
source activate minimap2

for i in {1,10,20,30,40,50,60,70,80,90,100};do
	minimap2 -y -ax lr:hq $ref $fastq_subsample/ecoli_exp_${SLURM_ARRAY_TASK_ID}_${i}x.fastq -t 5 -o $minimap2_subsample/ecoli_exp_${SLURM_ARRAY_TASK_ID}_${i}x_align.bam
	minimap2 -y -ax lr:hq $ref $fastq_subsample/ecoli_sta_${SLURM_ARRAY_TASK_ID}_${i}x.fastq -t 5 -o $minimap2_subsample/ecoli_sta_${SLURM_ARRAY_TASK_ID}_${i}x_align.bam
done

### index the sam file and convert to bam file

for i in {1,10,20,30,40,50,60,70,80,90,100};do
	samtools view -bhS $minimap2_subsample/ecoli_exp_${SLURM_ARRAY_TASK_ID}_${i}x_align.bam | samtools sort -T $minimap2_subsample/ecoli_exp_${SLURM_ARRAY_TASK_ID}_${i}x_align.sorted -o $minimap2_subsample/ecoli_exp_${SLURM_ARRAY_TASK_ID}_${i}x_align.sorted.bam
	samtools index $minimap2_subsample/ecoli_exp_${SLURM_ARRAY_TASK_ID}_${i}x_align.sorted.bam
	samtools view -bhS $minimap2_subsample/ecoli_sta_${SLURM_ARRAY_TASK_ID}_${i}x_align.bam | samtools sort -T $minimap2_subsample/ecoli_sta_${SLURM_ARRAY_TASK_ID}_${i}x_align.sorted -o $minimap2_subsample/ecoli_sta_${SLURM_ARRAY_TASK_ID}_${i}x_align.sorted.bam
	samtools index $minimap2_subsample/ecoli_sta_${SLURM_ARRAY_TASK_ID}_${i}x_align.sorted.bam
done

#########################################################
### extract only the primary mapped reads

for i in {1,10,20,30,40,50,60,70,80,90,100};do
	samtools view -F 0x900 -F 4 -bhS $minimap2_subsample/ecoli_exp_${SLURM_ARRAY_TASK_ID}_${i}x_align.bam > $minimap2_primary_subsample/ecoli_exp_${SLURM_ARRAY_TASK_ID}_${i}x_primary.bam
	samtools view -F 0x900 -F 4 -bhS $minimap2_subsample/ecoli_sta_${SLURM_ARRAY_TASK_ID}_${i}x_align.bam > $minimap2_primary_subsample/ecoli_sta_${SLURM_ARRAY_TASK_ID}_${i}x_primary.bam
done

#########################################################
### sort and index bam file
for i in {1,10,20,30,40,50,60,70,80,90,100};do
	samtools view -bhS $minimap2_primary_subsample/ecoli_exp_${SLURM_ARRAY_TASK_ID}_${i}x_primary.bam | \
	samtools sort -T $minimap2_primary_subsample/ecoli_exp_${SLURM_ARRAY_TASK_ID}_${i}x_primary.sorted -o $minimap2_primary_subsample/ecoli_exp_${SLURM_ARRAY_TASK_ID}_${i}x_primary.sorted.bam
	samtools index $minimap2_primary_subsample/ecoli_exp_${SLURM_ARRAY_TASK_ID}_${i}x_primary.sorted.bam
	
	samtools view -bhS $minimap2_primary_subsample/ecoli_sta_${SLURM_ARRAY_TASK_ID}_${i}x_primary.bam | \
	samtools sort -T $minimap2_primary_subsample/ecoli_sta_${SLURM_ARRAY_TASK_ID}_${i}x_primary.sorted -o $minimap2_primary_subsample/ecoli_sta_${SLURM_ARRAY_TASK_ID}_${i}x_primary.sorted.bam
	samtools index $minimap2_primary_subsample/ecoli_sta_${SLURM_ARRAY_TASK_ID}_${i}x_primary.sorted.bam
done

#########################################################
### adjust the modification information in bam file
source activate ont-modkit

for i in {1,10,20,30,40,50,60,70,80,90,100};do
    ## targeting m6A
    modkit adjust-mods $minimap2_primary_subsample/ecoli_exp_${SLURM_ARRAY_TASK_ID}_${i}x_primary.sorted.bam stdout --ignore m | modkit adjust-mods stdin $minimap2_primary_m6a_subsample/ecoli_exp_${SLURM_ARRAY_TASK_ID}_${i}x_primary_m6a.bam --ignore 21839
    modkit adjust-mods $minimap2_primary_subsample/ecoli_sta_${SLURM_ARRAY_TASK_ID}_${i}x_primary.sorted.bam stdout --ignore m | modkit adjust-mods stdin $minimap2_primary_m6a_subsample/ecoli_sta_${SLURM_ARRAY_TASK_ID}_${i}x_primary_m6a.bam --ignore 21839

    ## targeting m4C
    modkit adjust-mods $minimap2_primary_subsample/ecoli_exp_${SLURM_ARRAY_TASK_ID}_${i}x_primary.sorted.bam stdout --ignore m | modkit adjust-mods stdin $minimap2_primary_m4c_subsample/ecoli_exp_${SLURM_ARRAY_TASK_ID}_${i}x_primary_m4c.bam --ignore a
    modkit adjust-mods $minimap2_primary_subsample/ecoli_sta_${SLURM_ARRAY_TASK_ID}_${i}x_primary.sorted.bam stdout --ignore m | modkit adjust-mods stdin $minimap2_primary_m4c_subsample/ecoli_sta_${SLURM_ARRAY_TASK_ID}_${i}x_primary_m4c.bam --ignore a

    ## targeting m5C
    modkit adjust-mods $minimap2_primary_subsample/ecoli_exp_${SLURM_ARRAY_TASK_ID}_${i}x_primary.sorted.bam stdout --ignore 21839 | modkit adjust-mods stdin $minimap2_primary_m5c_subsample/ecoli_exp_${SLURM_ARRAY_TASK_ID}_${i}x_primary_m5c.bam --ignore a
    modkit adjust-mods $minimap2_primary_subsample/ecoli_sta_${SLURM_ARRAY_TASK_ID}_${i}x_primary.sorted.bam stdout --ignore 21839 | modkit adjust-mods stdin $minimap2_primary_m5c_subsample/ecoli_sta_${SLURM_ARRAY_TASK_ID}_${i}x_primary_m5c.bam --ignore a
done

#########################################################
### sort and index bam file
for i in {1,10,20,30,40,50,60,70,80,90,100};do
	## targeting m6a
	samtools view -bhS $minimap2_primary_m6a_subsample/ecoli_exp_${SLURM_ARRAY_TASK_ID}_${i}x_primary_m6a.bam | \
	samtools sort -T $minimap2_primary_m6a_subsample/ecoli_exp_${SLURM_ARRAY_TASK_ID}_${i}x_primary_m6a.sorted -o $minimap2_primary_m6a_subsample/ecoli_exp_${SLURM_ARRAY_TASK_ID}_${i}x_primary_m6a.sorted.bam
	
	samtools view -bhS $minimap2_primary_m6a_subsample/ecoli_sta_${SLURM_ARRAY_TASK_ID}_${i}x_primary_m6a.bam | \
	samtools sort -T $minimap2_primary_m6a_subsample/ecoli_sta_${SLURM_ARRAY_TASK_ID}_${i}x_primary_m6a.sorted -o $minimap2_primary_m6a_subsample/ecoli_sta_${SLURM_ARRAY_TASK_ID}_${i}x_primary_m6a.sorted.bam
	
	samtools index $minimap2_primary_m6a_subsample/ecoli_exp_${SLURM_ARRAY_TASK_ID}_${i}x_primary_m6a.sorted.bam
	samtools index $minimap2_primary_m6a_subsample/ecoli_sta_${SLURM_ARRAY_TASK_ID}_${i}x_primary_m6a.sorted.bam
	
	## targeting m4c
	samtools view -bhS $minimap2_primary_m4c_subsample/ecoli_exp_${SLURM_ARRAY_TASK_ID}_${i}x_primary_m4c.bam | \
	samtools sort -T $minimap2_primary_m4c_subsample/ecoli_exp_${SLURM_ARRAY_TASK_ID}_${i}x_primary_m4c.sorted -o $minimap2_primary_m4c_subsample/ecoli_exp_${SLURM_ARRAY_TASK_ID}_${i}x_primary_m4c.sorted.bam
	
	samtools view -bhS $minimap2_primary_m4c_subsample/ecoli_sta_${SLURM_ARRAY_TASK_ID}_${i}x_primary_m4c.bam | \
	samtools sort -T $minimap2_primary_m4c_subsample/ecoli_sta_${SLURM_ARRAY_TASK_ID}_${i}x_primary_m4c.sorted -o $minimap2_primary_m4c_subsample/ecoli_sta_${SLURM_ARRAY_TASK_ID}_${i}x_primary_m4c.sorted.bam
	
	samtools index $minimap2_primary_m4c_subsample/ecoli_exp_${SLURM_ARRAY_TASK_ID}_${i}x_primary_m4c.sorted.bam
	samtools index $minimap2_primary_m4c_subsample/ecoli_sta_${SLURM_ARRAY_TASK_ID}_${i}x_primary_m4c.sorted.bam
	
	## targeting m5c
	samtools view -bhS $minimap2_primary_m5c_subsample/ecoli_exp_${SLURM_ARRAY_TASK_ID}_${i}x_primary_m5c.bam | \
	samtools sort -T $minimap2_primary_m5c_subsample/ecoli_exp_${SLURM_ARRAY_TASK_ID}_${i}x_primary_m5c.sorted -o $minimap2_primary_m5c_subsample/ecoli_exp_${SLURM_ARRAY_TASK_ID}_${i}x_primary_m5c.sorted.bam
	
	samtools view -bhS $minimap2_primary_m5c_subsample/ecoli_sta_${SLURM_ARRAY_TASK_ID}_${i}x_primary_m5c.bam | \
	samtools sort -T $minimap2_primary_m5c_subsample/ecoli_sta_${SLURM_ARRAY_TASK_ID}_${i}x_primary_m5c.sorted -o $minimap2_primary_m5c_subsample/ecoli_sta_${SLURM_ARRAY_TASK_ID}_${i}x_primary_m5c.sorted.bam
	
	samtools index $minimap2_primary_m5c_subsample/ecoli_exp_${SLURM_ARRAY_TASK_ID}_${i}x_primary_m5c.sorted.bam
	samtools index $minimap2_primary_m5c_subsample/ecoli_sta_${SLURM_ARRAY_TASK_ID}_${i}x_primary_m5c.sorted.bam
done

#########################################################
### adjust the modification information in bam file
source activate ont-modkit

for i in {1,10,20,30,40,50,60,70,80,90,100};do
    ## targeting m6A
    modkit pileup $minimap2_primary_m6a_subsample/ecoli_exp_${SLURM_ARRAY_TASK_ID}_${i}x_primary_m6a.sorted.bam $modkit_pileup_primary_m6a/ecoli_exp_${SLURM_ARRAY_TASK_ID}_${i}x_primary_m6a.bed --motif A 0 --ref $ref
    modkit pileup $minimap2_primary_m6a_subsample/ecoli_sta_${SLURM_ARRAY_TASK_ID}_${i}x_primary_m6a.sorted.bam $modkit_pileup_primary_m6a/ecoli_sta_${SLURM_ARRAY_TASK_ID}_${i}x_primary_m6a.bed --motif A 0 --ref $ref

    ## targeting m4C
    modkit pileup $minimap2_primary_m4c_subsample/ecoli_exp_${SLURM_ARRAY_TASK_ID}_${i}x_primary_m4c.sorted.bam $modkit_pileup_primary_m4c/ecoli_exp_${SLURM_ARRAY_TASK_ID}_${i}x_primary_m4c.bed --motif C 0 --ref $ref
    modkit pileup $minimap2_primary_m4c_subsample/ecoli_sta_${SLURM_ARRAY_TASK_ID}_${i}x_primary_m4c.sorted.bam $modkit_pileup_primary_m4c/ecoli_sta_${SLURM_ARRAY_TASK_ID}_${i}x_primary_m4c.bed --motif C 0 --ref $ref

    ## targeting m5C
    modkit pileup $minimap2_primary_m5c_subsample/ecoli_exp_${SLURM_ARRAY_TASK_ID}_${i}x_primary_m5c.sorted.bam $modkit_pileup_primary_m5c/ecoli_exp_${SLURM_ARRAY_TASK_ID}_${i}x_primary_m5c.bed --motif C 0 --ref $ref
    modkit pileup $minimap2_primary_m5c_subsample/ecoli_sta_${SLURM_ARRAY_TASK_ID}_${i}x_primary_m5c.sorted.bam $modkit_pileup_primary_m5c/ecoli_sta_${SLURM_ARRAY_TASK_ID}_${i}x_primary_m5c.bed --motif C 0 --ref $ref
done

