#!/bin/bash
# Author: MDL, 2021
#
# Run: nohup bash phylostr_funct_summary_v2.sh A > phylostr_funct_summary_Athal.out &
 

SPECIES=$1
for c in {0..8} # 0, 4, 8
do
	time python3 phylostratum_function_enrichment_analysis_v2.py -s $SPECIES -C -c 0.$c -p > phylo_funct_summary_$SPECIES"_c_0_"$c"_go.out"
done

echo "Done."

