#!/bin/bash
# Author: Mirjana Domazet-Loso, May 2020
# Filter all fasta files
# For each proteome:
#	--> remove all proteins with empty protein sequences 
#	--> in case of multiple splicing variants, keep the protein with the longest sequence, and other variants (proteins with the same name) delete
#   --> rename all proteins to using into pgi number (using ordinal numbers + genome ID)
# Use clean.awk + filter.awk
#RUN: bash filter_E.sh taxIDNames.txt path_to_original_fasta_files path_to_edited_fasta_files > filter_E.log

declare -A array
# Idea was to improve this by reading list of names + taxIDs from a file into array
FILE=$1
# e.g. array[Vitis_vinifera]=29760 from input "29760	Vitis_vinifera"

while read line; do
	cnt=0
	line=${line##+([[:space:]])}          # Trim leading spaces (including tab, EOL)
	line=${line%%+([[:space:]])}          # Trim trailing spaces
    for word in $line; do
		#echo -e $word
		if [ $cnt -eq 0 ]; then
			taxID=$word
		else
			name=$word			
			array[$name]=$taxID
			echo $name ${array["$name"]} # prints space between $name and array elem
		fi
		let 'cnt = cnt + 1'
    done
done < $FILE

#done < taxID_names_Dmel.txt

# .fa, .faa, or .fasta files
#FILES=/home/mdloso/PROJECTS/ProteinsT/Bacillus2019/new_genomes_ff/*.fa*
OUTPUT_DIR=$3
cd $2
FILES=*.fa*

for f in $FILES
do

	# remove prefix ending with / (i.e. remove the directory name)
	tmp=${f#*/}   
	# stripping the longest match .* from back
	key=${tmp%%[.]*}
	#id=${array["$key"]}
	#key=Abeoforma_whisleri # != key1
	id=${array[$key]}
	echo -e "Key:" $key "Id:" $id
	
	outputFile=$OUTPUT_DIR/$id.ff	
	#echo -e $outputFile
	listFile=$OUTPUT_DIR/temp_$id.list	 
	
	#cp $f All_ff/$key.ff
	echo -e $f "\t" $listFile "\t" $outputFile		
	awk -v t=$id -f clean.awk $f > $listFile # get multiple splicings	
	awk -v t=$id -f filter_E.awk $listFile $f > $outputFile
	#exit
done

echo Done.
