#!/bin/bash
# Author: MDL, 2021-2022
# Run: nohup bash parseClusters_NEW.sh PATH_TO_ANALYSIS_FOLDER Hsap ./Clusters_2022 39 > parseClusters_Hsap.out &

INPUT_DIR=$1
SPECIES=$2
OUTPUT_DIR=$3
PS=$4

for c in {0..8}
do
	echo python3 parseClusters_for_FA_NEW.py -i $INPUT_DIR/results_0_$c/$SPECIES/dbAllPlus_clu_0_$c"_"$SPECIES"_out.txt" -o $OUTPUT_DIR/c_0_$c/$SPECIES"_2022_clusters.txt" -p $PS >> parseClu_for_FA_$SPECIES".out"
	time python3 parseClusters_for_FA_NEW.py -i $INPUT_DIR/results_0_$c/$SPECIES/dbAllPlus_clu_0_$c"_"$SPECIES"_out.txt" -o $OUTPUT_DIR/c_0_$c/$SPECIES"_2022_clusters.txt" -p $PS >> parseClu_for_FA_$SPECIES".out"
done
echo "Done."
