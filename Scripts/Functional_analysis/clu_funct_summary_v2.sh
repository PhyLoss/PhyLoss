#!/bin/bash
# Author: MDL, 2021-2022
# 
# Run: nohup bash clu_funct_summary_v2.sh H > clu_funct_summary_Hsap.out &

SPECIES=$1
for c in {0..8}
do
	time python3 cluster_function_enrichment_summary_go_v4.py -s $SPECIES -c 0.$c > clu_funct_summary_$SPECIES"_c_0_"$c"_go.out"
done

echo "Done."
