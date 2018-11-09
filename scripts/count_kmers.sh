#!/usr/bin/bash

KMC="$(dirname ${BASH_SOURCE[0]})/../build/KMC/kmc"


if [ $# -ne 3 ]; then
    echo "Usage $0 <fasta.gz> <k> <cutoff>"
    exit 1
fi

FILE="$1"
K="$2"
cutoff="$3"
num_threads=10


if [ ! -f $FILE ]; then
    echo "File does not exist"
    exit 1
fi

if [ ${FILE: -4} == '.bz2' ]; then
    echo "Error: gzip compressed file expected"
    exit 1
fi

mkdir -p "$FILE.cache"
/usr/bin/time -v $KMC -k$K -m10 -ci$cutoff -fm -t$num_threads $FILE $FILE.kmc $FILE.cache
rm -r "$FILE.cache"

