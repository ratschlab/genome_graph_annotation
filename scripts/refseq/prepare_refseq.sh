#!/usr/bin/bash

dataset_name="refseq_complete_family_16"
DATA_DIR="$HOME/metagenome/$dataset_name"

if [ ! -d $DATA_DIR ]; then
  echo "ERROR: bad directory"
  exit 1
fi

files_list="${dataset_name}.txt"

ls -1S "$DATA_DIR"/*.gz > $files_list

bsub -J "filter[1-$(cat $files_list | wc -l)]%500" \
  -W 24:00 -n 10 -R "rusage[mem=3400]" \
  "./count_kmers.sh \"\$(sed -n \${LSB_JOBINDEX}p $files_list)\" 16 1"

