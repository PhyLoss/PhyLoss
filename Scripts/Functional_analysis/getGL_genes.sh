#!/bin/bash
# Author: MDL, 2021-2022
# Run: bash getGL_genes.sh 9606 PATH_TO_TSV_FOLDER PATH_TO_PARENTS_FOLDER PATH_TO_LCA_FILE PATH_TO_OUTPUT_FOLDER H_sapiens_9606 > getGL_genes_H_sapiens_9606.out &
SPECIES=$1
TSV_DIR=$2
PARENTS_DIR=$3
LCA_FILE=$4
OUTPUT_DIR=$5
SPECIES_NAME=$6

for c in {0..8}
do
	C_DIR=results_0_$c
	C_DIR_2=0_$c
	GENES_FILE=$OUTPUT_DIR/$C_DIR_2/res_dbAll_$SPECIES_NAME"_0_"$c"_genes.txt"
	SUMMARY_FILE=$OUTPUT_DIR/$C_DIR_2/summary_GL_$SPECIES_NAME"_0_"$c".txt"
	echo $GENES_FILE $SUMMARY_FILE
	echo python3 getGL_genes.py -f $SPECIES -i $TSV_DIR/$C_DIR/db_clu_all.tsv -p $PARENTS_DIR/ -l $LCA_FILE -g $GENES_FILE -s $SUMMARY_FILE > out.txt
	time python3 getGL_genes.py -f $SPECIES -i $TSV_DIR/$C_DIR/db_clu_all.tsv -p $PARENTS_DIR/ -l $LCA_FILE -g $GENES_FILE -s $SUMMARY_FILE > out.txt
done
echo "Done."
