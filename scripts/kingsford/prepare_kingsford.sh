#!/usr/bin/bash

DATA_DIR="$HOME/metagenome/benchmark_kingsford/data_fasta"


for cutoff in {1,3,10,20,50}; do
  files_list="kingsford_${cutoff}.txt"
  bsub -J "filter[1-$(cat $files_list | wc -l)]%100" \
    -W 72:00 -n 10 -R "rusage[mem=3400]" \
    "../count_kmers.sh ${DATA_DIR}/\"\$(sed -n \${LSB_JOBINDEX}p $files_list)\".fasta.gz 20 $cutoff"
done

