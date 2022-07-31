#!/bin/bash
# Author: Mirjana Domazet-Lošo, 2018-2022
#
# For each genome:
#	--> compute parental nodes up to the root for each species
# Use treeTaxId.awk and auxiliary file nodes.txt
#RUN: bash computeParents.sh nodes.txt Parents > computeParents.log

PARENTS=$2
rm -rf $PARENTS
mkdir $PARENTS

FILES=All_ff/*.ff
FILE_NODES=$1

for f in $FILES
do

	# remove prefix ending with / (i.e. remove the directory name)
	tmp=${f#*/}   
	# stripping the longest match .* from back
	id=${tmp%%[.]*}	
	tFileName=$PARENTS/$id'parents.txt'
	
	# find all t's parents up to the root
	# e.g. awk -v t=7227 -f treeTaxId.awk nodes.dmp.fmt.new.sync > Parents/7227parents.txt
	echo -e $id "\t" $tFileName
	awk -v t=$id -f treeTaxId.awk $FILE_NODES > $tFileName
done

echo Done.
