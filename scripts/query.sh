#/usr/bin/env bash

file="$(dirname ${BASH_SOURCE[0]})/../tests/data/transcripts_1000.fa"

if [ -f annograph ]; then
  exe="./annograph"
else
  exe="$(dirname ${BASH_SOURCE[0]})/../build/annograph"
fi


if [ $# -ne 3 ]; then
  echo "Usage: $0 <graph> <annotation> <type>"
  exit 1
fi


$exe classify -v -i $1 -a $2 --anno-type $3 -o out.$1.$2.$3.tsv <(cat $file) 2>&1 | tee out.$1.$2.$3.log
echo "Done"
