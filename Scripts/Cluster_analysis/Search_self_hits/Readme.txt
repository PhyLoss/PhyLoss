# Korina Šimičević and Mirjana Domazet-Lošo, 2023/2024

######### Note #########
# Required Python modules: pandas, matplotlib

####################
1. Find all genes in the DB
	
	>python3 find_all_genes.py -d DB_PATH -o GENES_PATH
	
	--> generates a list of all genes within the database (from the DB folder)
	Input:
		DB_PATH --> folder containing taxId.faa files
	Output:
		GENES_PATH --> list of all genes ids saved in a text file, e.g. PATH/genes.csv
	
2. Running mmseqs easy-search + extracting self-hits
	
	>qsub esearch_all.pbs 
	
	--> runs an all-vs-all mmseq search for each .faa file in the database +
		extracts self-hits for each taxId to a file of the form "s_hits_taxId.txt", e.g. "s_hits_3981.txt" for taxId=3981
		using get_self_hits.py)
		Output:
			s_hits_*.txt files
	
3. 	Concatenate all s_hits_*.txt files to a single file 
	
	>cat PATH/s_hits_*.txt > PATH/s_hits_all.txt
	
	
4. Find genes without self-hits
	
	>python3 find_genes_with_no_self_hits.py -i HITS_PATH -g GENES_PATH -o MISSES_PATH
	
	--> generates a list of all the proteins lacking a self-hit from the mmseqs search.

	Input e.g.: 
		GENES_PATH= "genes.csv"
		HITS_PATH = "s_hits_all.txt"
	Output e.g.:
		MISSES_PATH = "self_misses.tsv"

5. Filter DB, i.e. 
	for each taxId.faa file: delete proteins without self-hits
	
	> cp -r DB_2024 DB_2024_UPD
	> python3 filter_db.py -i FILTER_PATH/self_misses.tsv -d DB_2024_UPD
	
	--> from each taxId.faa in the database folder delete protein sequences that lack a self-hit
