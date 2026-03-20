#!/bin/bash -l

#########################################################
source activate mijamp

for code in m6a m4c m5c; do
    for file in $bam_dir/$code/*.sorted.bam; do
        python $preprocess -b $file -g $ecoli_ref -t 2 -o $preprocess_bam/${file%%.*}
        python $motif -f $preprocess_bam/${file%%.*}
    done
done