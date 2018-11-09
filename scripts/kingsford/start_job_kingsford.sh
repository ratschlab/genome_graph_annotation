bsub -J assemble_hash_king -W 72:00 -n 1 -R "rusage[mem=250000] span[hosts=1]" \
    "cat kingsford_kmc.txt | gtime -v ../../build_release/annograph build -v -k 20 -c --kmc -o ~/big_graph/kingsford/kingsford_k20_canonical 2>&1 | tee ~/big_graph/kingsford/log_assemble_hash_kingsford_k20.txt"

for l in {1..6}; do bsub -J "annotate_kingsford_canonical_$l" -W 50:00 -n 1 -R "rusage[mem=300000] span[hosts=1]" \
    "cat ~/annotation_schemes/scripts/kingsford/kingsford_kmc_$l.txt | gtime -v ~/annotation_schemes/build_annotate_king/annograph annotate -v --kmc --anno-filename -i ~/big_graph/kingsford/kingsford_k20_canonical.orhashdbg -o ~/big_graph/kingsford/kingsford_k20_canonical_$l | tee ~/big_graph/kingsford/log_annotate_kingsford_can_k20_$l.txt"; \
done
