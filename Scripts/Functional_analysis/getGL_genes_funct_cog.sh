#!/bin/bash
# Author: MDL, 2023
# bash getGL_genes_funct.sh 9606 PATH/Data H_sapiens_9606
#
#
SPECIES=$1
WORKING_DIR=$2
TSV_DIR=$WORKING_DIR/tsv_new
PARENTS_DIR=$WORKING_DIR/Parents
LCA_FILE=$WORKING_DIR/allLCA.txt
OUTPUT_DIR=$WORKING_DIR/ClusterAnalysis/COG
EMAPPER_DIR=$WORKING_DIR/eggnogResultsNew
SPECIES_NAME=$3


for c in {0..0}
do
	C_DIR=results_0_$c
	C_DIR_2=0_$c
	GENES_FILE=$OUTPUT_DIR/$C_DIR_2/res_dbAll_$SPECIES_NAME"_0_"$c"_genes_funct.txt"
	OUTPUT_FILE=$OUTPUT_DIR/$C_DIR_2/res_dbAll_$SPECIES_NAME"_0_"$c"_clusters_funct.txt"
	SUMMARY_FILE=$OUTPUT_DIR/$C_DIR_2/summary_GL_$SPECIES_NAME"_0_"$c"_funct.txt"
	echo $GENES_FILE $SUMMARY_FILE
	echo python3 getGL_genes_funct.py -f $SPECIES -i $TSV_DIR/$C_DIR/db_clu_all.tsv -p $PARENTS_DIR/ -l $LCA_FILE -e $EMAPPER_DIR -C -g $GENES_FILE -s $SUMMARY_FILE -o $OUTPUT_FILE > out_funct.txt
	time python3 getGL_genes_funct.py -f $SPECIES -i $TSV_DIR/$C_DIR/db_clu_all.tsv -p $PARENTS_DIR/ -l $LCA_FILE -e $EMAPPER_DIR -C -g $GENES_FILE -s $SUMMARY_FILE -o $OUTPUT_FILE > out_funct.txt
done
echo "Done."
