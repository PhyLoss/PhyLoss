# Author: MDL, 20 July, 2018
# For each genome:
#	--> get the list of genes with the longest sequence (multiple splicing!) from the .list file and store it in an array
#	--> filter genes from .fa file: get only those in the array and transform their id using pdi number
#
# Example: awk -v t=9606 -f filter_E.awk path_to_edited_fasta_files/temp_9606.list path_to_original_fasta_files/Homo_sapiens.fa > path_to_edited_fasta_files/7227.ff
BEGIN{
	# get taxId as program input -v
	cnt = 0;
	taxId = t;
}
FNR == NR {
	list[$1]=$1;
	next; # is this necessary?
}
/>/{ #print $0, a[$1]	
	if ($1 in list) {
		++ cnt;
		proteinId = taxId "" cnt;
		proteinFullId = sprintf(">pid|%019d|tx|%d|\t%s", proteinId, taxId, substr($1, 2));
		print proteinFullId
		printSeq = 1;
	}
	else {
		printSeq = 0;
	}
}
!/>/{
	if (printSeq) print;
}
