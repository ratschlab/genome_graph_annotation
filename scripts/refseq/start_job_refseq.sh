cat ./refseq_complete_family_15_kmc.txt \
    | ~/annotation_schemes/build_release/annograph build -k 15 --complete -o ~/big_graph/binrel_annotations/refseq_k15

mkdir ~/big_graph/binrel_annotations/refseq_k15_annotation

bsub -J "annotate_refseq[1-3173]%200" -W 20:00 -n 1 -R "rusage[mem=5000] span[hosts=1]" \
    "file=\"\$(sed -n \${LSB_JOBINDEX}p refseq_complete_family_15_kmc.txt)\"; x=\$(basename \${file%.*.*.*.*}); gtime -v ~/annotation_schemes/build_release/annograph annotate -v -p 2 --anno-label \${x} --kmc --complete -i ~/big_graph/binrel_annotations/refseq_k15 -o ~/big_graph/binrel_annotations/refseq_k15_annotation/annotation_refseq_family_k15_\${x} \$file"

bsub -J assemble_hash_refseq_15 -W 72:00 -n 1 -R "rusage[mem=300000] span[hosts=1]" \
    "cat refseq_complete_family_15_kmc.txt | gtime -v ../../build_release/annograph build -v -k 15 --kmc -o ~/big_graph/refseq/hash_refseq_k15 2>&1 | tee ~/big_graph/refseq/log_assemble_hash_refseq_k15.txt"
