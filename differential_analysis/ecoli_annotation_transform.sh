#!/bin/bash -l
#SBATCH --job-name="ecoli_annotation_transform"
#SBATCH --partition=general
#SBATCH -o ecoli_annotation_transform.o
#SBATCH -e ecoli_annotation_transform.e

#########################################################
module load anaconda3

############################################################################
### clean the gff file and get gene promoter file
awk 'BEGIN{IFS=OFS="\t"}{print $0}' $annotation_dir/GCF_000008865.2_ASM886v2_genomic.gff | grep -v "^#" | awk 'BEGIN{IFS=OFS="\t"} $3!="region" && $4!="1" {print $0}' > $annotation_dir/GCF_000008865.2_ASM886v2_genomic_clean.gff

## use agat to transform gff file to tsv file
source activate agat
agat_convert_sp_gff2tsv.pl --gff $annotation_dir/GCF_000008865.2_ASM886v2_genomic_clean.gff -o $annotation_dir/GCF_000008865.2_ASM886v2_genomic_clean.tsv

## subset to gene region
awk -F"\t" '$3=="gene"&& $2=="RefSeq" {print $1,$3,$4,$5,$7,$9,$15,$16}' $annotation_dir/GCF_000008865.2_ASM886v2_genomic_clean.tsv |\
awk 'BEGIN{IFS=OFS="\t"} {if ($5=="1") {$5="+"} else {$5="-"}; print $0}' > $annotation_dir/GCF_000008865.2_ASM886v2_genomic_gene.tsv

## subset to promoter region
awk 'BEGIN{IFS=OFS="\t"} {if($5=="+") {$4=$3-7; $3=$3-35} else {$3=$4+7; $4=$4+35}; print $0}' $annotation_dir/GCF_000008865.2_ASM886v2_genomic_gene.tsv \
> $annotation_dir/GCF_000008865.2_ASM886v2_genomic_gene_extended_promoter.tsv
