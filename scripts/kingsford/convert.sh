#!/usr/bin/bash


annotation="$HOME/metagenome/benchmark_kingsford_results/binrel_annotations/kingsford_k20_canonical"
exe="../../build/annograph"


#/usr/bin/time -v ./annograph transform_anno --anno-type flat -o $annotation $annotation
#./annograph stats --anno-type flat -a $annotation
#echo ""

#/usr/bin/time -v ./annograph transform_anno --anno-type bin_rel_wt_sdsl -o $annotation $annotation
#./annograph stats --anno-type bin_rel_wt_sdsl -a $annotation
#echo ""

#/usr/bin/time -v ./annograph transform_anno --anno-type bin_rel_wt -o $annotation $annotation
#./annograph stats --anno-type bin_rel_wt -a $annotation
#echo ""

#/usr/bin/time -v $exe transform_anno --anno-type brwt -o $annotation $annotation
#$exe stats --anno-type brwt -a $annotation
#echo ""

/usr/bin/time -v $exe transform_anno --anno-type brwt --greedy -o ${annotation}_greedy $annotation
$exe stats --anno-type brwt -a $annotation
echo ""

#/usr/bin/time -v $exe transform_anno --anno-type rbfish -o $annotation $annotation
#$exe stats --anno-type rbfish -a $annotation
#echo ""

