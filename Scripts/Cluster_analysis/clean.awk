# Author: MDL, 20 July, 2018
# For each genome in Eukaryota_ensembl:
#	--> in case of multiple splicing variants, keep the protein with the longest sequence, and other variants (proteins with the same name) delete
#	--> here: multiple variants are found and deleted using filter.awk
# Example: awk -v t=9606 -f clean_E.awk Homo_sapiens.fa > path_to_edited_fasta_files/temp_7227.list
#
BEGIN{
	prev = 0; # previous gene
}
/>/{
	# store the length of the previous sequence/protein
	if (!prev) {
		prev = 1;
	}
	else {
		if (!arrLen[geneId] || seqLen > arrLen[geneId]) {
			arrLen[geneId] = seqLen;
			arrHeaders[geneId] = geneHeader;			
		}	
		++ numCopies[geneId];		
	}
	
	# extract gene id; 
	# e.g. extract ENSG00000283951.1 from a header: ">ENSP00000491436.1 pep chromosome:GRCh38:CHR_HSCHR19KIR_7191059-1_CTG3_1:54850452:54867240:1 gene:ENSG00000283951.1 transcript:ENST00000640782.1 gene_biotype:protein_coding ..."
	startPos = index($0, "gene:") + 5; # + 5 to skip "gene:"
	endPos = index($0, "transcript:") - 2;		
	geneId = substr($0, startPos, endPos - startPos + 1);
	
	# in case geneId is empty string, than set geneID to $1; this can happen, if the file's source is NCBI, since they don't have the above naming convention
	if (!length(geneId)) geneId = $1;
	seqLen = 0;	
	geneHeader = $1;
}
!/>/{
	# print sequence
	seqLen += length($0);
}
END{
	# last gene
	if (!arrLen[geneId] || seqLen > arrLen[geneId]) {
		arrLen[geneId] = seqLen;
		arrHeaders[geneId] = geneHeader;			
	}	
	++ numCopies[geneId];		
	
	# proteins longer than maxLen, should be deleted from the final versions of proteomes (i.e. excluded from .ff files)
	maxLen = 10000;
	
	# print results
	for (geneId in arrHeaders) {
		if (arrLen[geneId] >= maxLen) { # print list of these proteins to a separate file, e.g. report_TooLongProteins.txt
			print "taxID = " t "\t" arrHeaders[geneId] "\t" geneId "\tLength = " arrLen[geneId] >> "./report_TooLongProteins.txt";
		}
		else {
			print arrHeaders[geneId] "\t" geneId "\t" numCopies[geneId];
		}
	}
}
