#!/bin/bash
# Author: MDL, 2021-2022
# # Run: nohup bash clu_funct_analysis_v2.sh H > clu_funct_analysis_Hsap.out &

SPECIES=$1
for c in {0..8}
do
	#time python3 cluster_function_enrichment_analysis_v4.py -s H -G -c 0.8 -S > clu_funct_analysis_HS_2022_c_0_8_go.out
	time python3 cluster_function_enrichment_analysis_v4.py -s $SPECIES -C -c 0.$c > clu_funct_analysis_$SPECIES"_c_0_"$c"_go.out"
done

echo "Done."
