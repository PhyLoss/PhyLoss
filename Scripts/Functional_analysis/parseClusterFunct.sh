#!/bin/bash
# Author: MDL, 2023
# bash parseClusterFunct.sh 39 Hsap PATH/Data/ClusterAnalysis/COG PATH/Data/COGs.txt
#
#
MAX_PS=$1
SPECIES=$2
OUTPUT_DIR=$3
FUNCT_FILE=$4

for c in {0..8}
do
	C_DIR_2=0_$c

	time python3 parseClusterFunct.py -P $1 -S $2 -i $3"/0_"$c/res_dbAll_H_sapiens_9606_0_$c"_clusters_funct.txt" -c $4 -o $3/"0_"$c > cog_out.txt
done
echo "Done."
