# Author: Mirjana Domazet-Loso
# 2022 - 2023
# Contact: Mirjana.Domazet-Loso@fer.hr
#------------------------------------------
#
# Reference:
# 	Mirjana Domazet-Lošo, Tin Široki, Korina Šimičević, Tomislav Domazet-Lošo, 2023, 
# 		Macroevolutionary dynamics of gene family gain and loss along multicellular eukaryotic lineages
# 
# ------------------------------------------
# Code contributors:
#	Mirjana Domazet-Lošo, Tin Široki and Korina Šimičević
#	University of Zagreb Faculty of Electrical Engineering and Computing, Zagreb, Croatia

-- PREREQUISITES --
	- installed programs: 
		python 3

----------- DATA -----------
# create a folder "Data" in your working directory, e.g.
	> cd YOUR_PATH
	> mkdir Data
	> mkdir Data/DB	

***** DATA REQUIRED FOR THE CLUSTER ANALYSIS PIPELINE *****

(1) a database that contains protein sequences in fasta format, 
	i.e. a list of fasta files (a single file for each species) in the format taxId.faa ("taxID" is replaced with a true taxID of the species)
	We provide the script "convertAll.py" to convert fasta files to the taxId.faa format (see below).
	--> save in YOUR_PATH/Data/DB
		
		(*) We used the following format for fasta files (required in subsequent scripts):
			- a separate file for each species named taxID.faa ("taxID" is replaced with the species true taxID) 
				e. g. for H. sapiens (taxID = 9606): 9606.faa
			- each .faa file contains protein sequences for the corresponding species

			- header lines in our .faa files have the following format, e.g.: 
			>pid|0000000000000012601|tx|1260|	KXA09781
			
			where	
				pid|0000000000000012601|tx|1260| represents a protein name (taxID/tx = 1260; proteinID/pid = 0000000000000012601)
				KXA09781 represents the original protein name		
		
		Notes:
		(*) To get the above format for protein sequences (i.e. pid|0000000000000012601|tx|1260|), 
			you may use the script "convertAll.py" that generates the .faa files,
			whose inputs are the original fasta files and a file containg all genomes' taxID-s e.g. taxIDNames.txt (see (2))
			e.g.
			> python3 time python3 convertAll.py -t PATH_TO/taxIDNames.txt -i PATH_TO_INPUT_FOLDER [-r PATH_TO_REPORT_FOLDER] -o PATH_TO_OUTPUT_FOLDER
			
			- the script "convertAll.py" calls "clean.py" for each taxId
			- the script "clean.py":
			a) deletes protein seqeunces shorter than MIN_LEN (10) and longer than MAX_LEN (10000) a.a.
			b) deletes protein seqeunces containing "*" in the middle of the sequence
			c) in case of multiple splicing variants, it keeps the longest sequence, 
				and other variants (protein sequences with the same name) are deleted
		
		(**) We additionally removed the protein sequences without self-hits.
				To find the protein sequences without self-hits, we compared each sequence to itself using MMseqs2 easy-search with the following options:
				--exhaustive-search -e 10e-3 --max-seqs 10000000
				(MMseqs2; Steinegger M and Soeding J. 2017 [1])
		
		(***) Finally, we applied the BUSCO analysis (BUSCO v5.4.3 analysis [2]) to our taxId.faa fasta files and replaced or removed low-quality files.		
		
(2) taxIDNames.txt
	(see an example file in PhyLoss/ExampleData)
	--> save in YOUR_PATH/Data
		
		- the list of all species in the database in the format: 
			taxID	species_name
		E.g.
			9606	Homo_sapiens
			7227	Drosophila_melanogaster
			
(4) names.txt 
	(see an example file in PhyLoss/ExampleData)
	--> save in YOUR_PATH/Data
		- the list of all nodes in the input phylogeny: from the species/leaves in the tree up to the root (e.g. "Cellular organism")
		
		- the required format: 
			nodeID	node_name
		E.g.
			9606	Homo_sapiens
			7227	Drosophila_melanogaster
			131567	Cellular_organisms

(5) nodes.txt
	(see an example file in PhyLoss/ExampleData)
	--> save in YOUR_PATH/Data

		- the list of all nodes' pairs in the input phylogeny in the format:
			nodeID	parent_nodeID
		E.g. 
			9606	131567
		where		
			9606	Homo_sapiens
			131567	Cellular_organisms





***** DATA REQUIRED FOR THE FUNCTIONAL ANALYSIS PIPELINE *****

(5) The clusters computed using MMseqs2 [1] and the input data (1) - (4)
	The clusters are computed and analyzed using scripts listed in the section "CLUSTER ANALYSIS PIPELINE".

(6) The eggnog mapper files --> save in YOUR_PATH/Data/eggnogResults
		(see an example file in PhyLoss/ExampleData/chunk_sample.emapper.annotations)
		(eggNOG-mapper v2; Cantalapiedra et al. 2021 [3])
		
		Note:
		(*) the eggnog mapper functional annotations (COG and GO annotations) should be precomputed
			--> each file (emapper.annotations) contains list of genes withe their GO and COG annotations
			- for more information, see "FUNCTIONAL ANALYSIS PIPELINE"
		
		(**) We used eggNOG-mapper emapper.py with the options:
				-m diamond --no_annot --no_file_comments
		
(7) go-basic.obo --> save in YOUR_PATH/Data
		- the GO-basic functions file
		(see http://geneontology.org/docs/download-ontology/ [4])

(8) COGs.txt --> save in YOUR_PATH/Data
		- the COG functional categories file
		(see http://clovr.org/docs/clusters-of-orthologous-groups-cogs/ [5])





----------- CLUSTER ANALYSIS (CA) PIPELINE -----------

PREREQUISITES:
	- install programs: 
		python 3.6 or later (http://www.python.org/getit/)
	
	- for MMseqs2 cluster analysis:
		install program MMseqs2 (https://github.com/soedinglab/MMseqs2) [1]

REQUIREMENTS: 
	- list of SCRIPTS files (see the folder "PhyLoss/Scripts")
		- getAllLCA.py
		- get_PS_short_names.py
		- getGL.py
		- submit_cluster_c.sge
		- treeTaxId_all.py
		- upd_NodesWith1Child.py		

(CA-1) Concatenate all DB fasta files to a single faa file "db_all.faa" (in the subolder "DB")	
	> cat YOUR_PATH/Data/DB/*.faa > YOUR_PATH/Data/DB/db_all.faa


(CA-2) Position to the folder "YOUR_PATH/Data"
	> cd YOUR_PATH/Data


(CA-3) 
	a) Run the mmseqs createdb [1] program to create the "db_all" from "db_all.faa"
	> mmseqs createdb db_all.faa db_all
	
	b) Run the mmseqs cluster [1] program 
	(we used the following options: -e 0.001 --max-seqs 400 --cluster-mode 1 --cov-mode 0; -c 0.0 to 0.8)
	> bash submit_cluster_c.sge &
	Alternatively: run the script using qsub
	
	Output: 
		- the list of subfolders named "YOUR_PATH/Data/results_0_0" to "YOUR_PATH/Data/results_0_8"
		each containing the mmseqs output tsv file "db_clu_all.tsv"
	
	
(CA-4) If there are parental nodes with only one leaf, then nodes, names and taxa file should be updated using upd_NodesWith1Child.py:  
	--> generate UPD files (where the parents/inner nodes are connected to their only child)
	(required for the later gene families gain/loss analysis)
	
	> python3 upd_NodesWith1Child.py -i YOUR_PATH/Data/nodes.txt \
										-j YOUR_PATH/Data/names.txt \
										-t YOUR_PATH/Data/taxIDNames.txt \
										-o YOUR_PATH/Data/nodes_with1Child.txt \
										-n YOUR_PATH/Data/nodes_UPD.txt \ 
										-m YOUR_PATH/Data/names_UPD.txt \
										-s YOUR_PATH/Data/taxIDNames_UPD.txt
	INPUT FILES:
		nodes.txt
		names.txt
		taxIDNames.txt
	
	OUTPUT FILES: 
		taxIDNames_UPD.txt
		names_UPD.txt
		nodes_UPD.txt
		nodes_with1Child.txt


(CA-5a) Create a list of parents for each focal species (in "Parents" subdirectory);
	The "parents" are parental nodes from the focal species up to the root (i.e. "Cellular organisms") 
	(required for the later gene families gain/loss analysis)
	
	> mkdir YOUR_PATH/Data/Parents
	> time python3 treeTaxId_all.py -i YOUR_PATH/Data/names_UPD.txt -j YOUR_PATH/Data/nodes_UPD.txt \
									-t YOUR_PATH/Data/taxIDNames_UPD.txt -o YOUR_PATH/Data/Parents/
	
	INPUT FILES:
		taxIDNames_UPD.txt
		names_UPD.txt
		nodes_UPD.txt

	OUTPUT FOLDER: 
		the subfolder "Parents" comprising a list of parent files 
		--> a file for each taxID, e.g. "Parents/7227parents.txt"

	
(CA-5b) Concatenate parents' files into a single file ("allParents.txt")
	> cat YOUR_PATH/Data/Parents/*parents.txt > YOUR_PATH/Data/allParents.txt
	

(CA-6) Compute the LCA (Last Common Ancestor) matrix for all taxID (n x n; in our case: 667 x 667); 
	note that the taxon's names are truncated to 10 char-s:
	
	> python3 getAllLCA.py -i YOUR_PATH/Data/allParents.txt -o YOUR_PATH/Data/allLCA.txt  

	E.g. output (allLCA.txt):
									880073     526226     526227    1924735	... 
								Caldithrix Gordonia_b Meiothermu Salpingoec 
			880073 (Caldithrix)         -1          2          2     131567 
			526226 (Gordonia_b)          2         -1          2     131567 
			526227 (Meiothermu)          2          2         -1     131567 
			...

(CA-7) Compute list of phylostrata (PS) names starting with (1) "Cellular org." (use existing abbrev. for PS 1 to 4)
	e.g.
	> python3 get_PS_short_names.py -i YOUR_PATH/Data/listMetazoa.txt \
									-p YOUR_PATH/Data/Parents/ \
									-o YOUR_PATH/Data/PS_names_short/ 
	INPUT FILES:
		YOUR_PATH/Data/Parents/ (step (CA-5a))
		YOUR_PATH/Data/listMetazoa.txt
		(a list of species with their names, e.g.
			9606	Homo_sapiens)

	OUTPUT FILE:
		YOUR_PATH/Data/PS_names_short/ (step (CA-5a))		
		--> the output files in the folder will be named after the taxon they refer to in the form:
			"Homo_sapiens_9606.txt" (name + "_" + taxId + ".txt")

(CA-8) Compute the gene families' gain/loss from the perspective of a given focal species/taxID 
	(e.g. if taxID = 9606, then compute the gene families' gain/loss for H. sap.)
	We apply the Dollo's parsimony --> the specific rules are listed in the script "getGL.py".
	
	e. g. for H. sapiens and c = 0.0 (i.e. the input subfolder "YOUR_PATH/Data/results_0_0")
	> python3 getGL.py -f 9606 -p YOUR_PATH/Data/Parents/ -i YOUR_PATH/Data/results_0_0/db_clu_all.tsv \
								-l YOUR_PATH/Data/allLCA.txt -o YOUR_PATH/Data/res_dbAll_H_sap_0_0.txt -s YOUR_PATH/Data/summary_GL_H_sap_0_0.txt

	INPUT FILES:
		YOUR_PATH/Data/Parents/ (step (CA-5a))
		YOUR_PATH/Data/results_0_0/db_clu_all.tsv  (mmseqs cluster results; step (CA-3))
		YOUR_PATH/Data/allLCA.txt (step (CA-6))
		
	OUTPUT FILES:
		(i) YOUR_PATH/Data/summary_GL_H_sap_0_0.txt
		(ii) YOUR_PATH/Data/res_dbAll_H_sap_0_0.txt
		
		-----		
		(i) a summary output file for each focal species and each c (e.g. summary_GL_H_sap_0_0.txt ):
		-------------	
			ps	GF_Gained	GF_Lost
			1	8393	0
			2	193	0
		...
		-------------
		
		(ii¸) a report file containing cluster analysis for clusters involved in gene families gain/loss for a focal species and c
		-------------
		6	pid|0000000001022416913|tx|10224|	8	15	19
		...
		-------------
		where
			6 is the ordinal number of cluster among all clusters in .tsv file
			pid|0000000001022416913|tx|10224| is a cluster representative
			8 is the number of genes in the cluster
			15 is the phylostratum_gain from the perspective of the focal species
			19 is the phylostratum_lost from the perspective of the focal species
			

(CA-9) 
	Compute the gene families' gain/loss from the perspective of all focal species/taxID (excluded only bacterial and side-branch archaeal clusters).
	The script "nodesGTL.py" computes:
		(i) the total number of gene families (T)
		(ii) the gained number of gene families (G)
		(iii) the lost number of gene families (L)
		across the all nodes (including leaves) in the phylogeny.
	Prerequisite: 
		For each taxon (a leaf) in the phylogeny, the script "getGL.py" (see CA-8) should be called with the options -s to get the summary file with
		the counts of gene families gained/loss along the tree from the perspective of that taxon, i.e. that focal species.
		We use the Dollo's parsimony as described in (CA-8). 
	Note: 
		Singleton clusters are excluded from the analysis.

	Input: 
		(i) an input folder with subfolders containing gene families total/gain/loss (EDIT!!)
		(ii) the parents input folder, e.g. "Parents" (see (CA-5a))
		(iii) the names file (see DATA REQUIRED FOR THE CLUSTER ANALYSIS PIPELINE - (4); e.g. "names.txt" )

    Output: 
		A list of files in an output folder, where a single file corresponds to a single c-value (0.0 to 0.8) in the form 
		e.g.
		131567	Cellular_organisms	10942	10942	n.a.
		
		where 
			"131567" is the node id, 
			"Cellular_organisms" is the node name, 
			"10942" (the 3rd column) is the total number of gene families at the node
			"10942" (the 4th column) is the number of gene families gained at the node
			"n.a." (the 5th column) is the number of gene families lost at the node 
				(in this example "n.a."/"not available", since no gene families were at lost at this node)	

	> python3 nodesGTL.py -i INPUT_FOLDER -o OUTPUT_FOLDER -p PARENTS_FOLDER -n NAMES_FILE

	

----------- FUNCTIONAL ANALYSIS PIPELINE -----------
PREREQUISITES:
	- data computed using "Cluster analysis pipeline" (see above)

	- eggnog mapper annotations of the genes in your database
		install EggNog mapper (https://github.com/eggnogdb/eggnog-mapper)	
	
REQUIREMENTS: 
	- create your analysis folder
		> mkdir YOUR_PATH/Data/ClusterAnalysis	
	
	- create the output subfolders (for the COG and GO analysis; then the subfolder for each focal value)
		e.g. the subolder "Hsap" for H. sapiens
		> mkdir YOUR_PATH/Data/ClusterAnalysis/COG
		> mkdir YOUR_PATH/Data/ClusterAnalysis/GO
		> mkdir YOUR_PATH/Data/ClusterAnalysis/COG/Hsap
		> mkdir YOUR_PATH/Data/ClusterAnalysis/GO/Hsap	
	
	
(FA-1) Compute clusters with the functional annotations using getGL_genes_funct.py
	--> the output includes COG or GO annotations for clusters (and, if -g option is set, genes annotations for the genes in clusters are also printed) 
	--> the script is similar to "getGL.py" from the cluster analysis pipeline, 
		but in this case (besides the info on gene families gain/loss) for each gene and cluster COG or GO annotations are computed
			
	e.g. 
		H. sapiens/taxID = 9606 --> -f 9606
		c = 0.0 --> the input tsv file: -i YOUR_PATH/Data/results_0_0/db_clu_all.tsv
		for the COG analysis: -C
		for the GO analysis: -G														
	
	> mkdir mkdir YOUR_PATH/Data/ClusterAnalysis/COG/0_0
		--> create the subfolder "0_0" for c-value = 0.0
		
	> python3 getGL_genes_funct.py -f 9606 \
					-p YOUR_PATH/Data/Parents/ \
					-i YOUR_PATH/Data/results_0_0/db_clu_all.tsv \
					-l YOUR_PATH/Data/allLCA.txt \
					-C
					-o YOUR_PATH/Data/ClusterAnalysis/COG/0_0/Hsap_0_0_clusters_funct.txt \
					-g YOUR_PATH/Data/ClusterAnalysis/COG/0_0/Hsap_0_0_genes_funct.txt \
					-s YOUR_PATH/Data/ClusterAnalysis/COG/0_0/summary_GL_Hsap_0_0_funct.txt
	
	INPUT FILES:
		YOUR_PATH/Data/Parents/ (step (CA-5a))
		YOUR_PATH/Data/results_0_0/db_clu_all.tsv  (mmseqs cluster results; step (CA-3))
		YOUR_PATH/Data/allLCA.txt (step (CA-6))
		
	OUTPUT FILES:
		(i) YOUR_PATH/Data/ClusterAnalysis/COG/0_0/summary_GL_H_sap_0_0.txt
		(ii) YOUR_PATH/Data/ClusterAnalysis/COG/0_0/Hsap_0_0_clusters_funct.txt
		(iii) YOUR_PATH/Data/ClusterAnalysis/COG/0_0/Hsap_0_0_genes_funct.txt
	
	Output:
		(i) a summary output file for each focal species and each c (e.g. summary_GL_H_sap_0_0.txt ):
			e.g.
			
			ps	GF_Gained	GF_Lost
			1	8393	0
			2	193	0
			...
	
		(ii) Hsap_0_0_clusters_funct.txt contains a list of clusters with assigned COG or GO annotations (here: COG annotations)
		- in the example below:
			- the first row contains a cluster_ordinal_number (22), a cluster representative/geneId (pid|0000000000102284460|tx|10228|),
				the number of genes in the cluster (313), phylostratum_gain (10), phylostratum_loss (-) (from the perspective of the focal species)
			- the following rows are COG (or GO) annotations assigned to the cluster (if obtained using emapper)
			(in this example: the 'O' annotation is assigned to 286 genes in the cluster, and the 'U' annotation is assigned to 26 genes in the cluster;
			there are in total two (2) annotations assigned to the cluster)
	
			e.g.
			22	pid|0000000000102284460|tx|10228|	313	10	-
			O	286	2
			U	26	2
		
		(iii) Hsap_0_0_genes_funct.txt contains a list of clusters, where each cluster info is 
			accompanied by the list of its genes and their COG or GO annotations (if obtained using emapper) (in this example: COG annotations)
		
			e.g.
			3003	pid|0000000005745667408|tx|574566|	3	6	13
			pid|0000000005745667408|tx|574566|	M
			pid|0000000000860183965|tx|86018|	G
			pid|0000000094636211029|tx|946362|	M

	Note:
	(*) getGL_genes_funct.py has an option -S, which refers to the "strict function assignment".
		If the option -S is set, the clusters will be assigned only those functions that are assigned to at least half of the cluster's members.
		By default, this options is excluded.
	
	(**) Printing the detailed list of all genes (the -g option) and their functions is optional, 
		since, in the case of numerous GO annotations, it takes longer time to compute and the output is large.


(FA-2) Parse results of the clusters' functional analysis and output the results in the form of summary tables and total counts of
	COG or GO annotations to the output folder --> using "parseClusterFunct.py"
	Two output files are produced for gene families gain analysis and the other two files one for gene families loss analysis.
	
	> python3 parseClusterFunct.py -P 40 -S Hsap \
					-i PATH/ClusterAnalysis/COG/0_0/Hsap_0_0_clusters_funct.txt  \
					-c PATH/COGs.txt \
					-o PATH/ClusterAnalysis/COG/0_0
	
	OPTIONS: -P 40 (H. sapiens has 40 phylostrata in our analysis)
			-S Hsap (-S species_short_name; the analysis was adjusted for the four focal species H.sap, D.mel, A.thal, S.cer, 
						but can be extended to any other focal species; in this example: the focal species is H. sap.)
			
	INPUT FILES:
		PATH/ClusterAnalysis/COG/0_0/Hsap_0_0_clusters_funct.txt (step ((FA-1)))
		PATH/COGs.txt (see "Data required for functional analysis"(8))
		
	OUTPUT FILES:
		(i)  PATH/ClusterAnalysis/COG/0_0/phyl_multi_fun_table_Hsap_gain_cog.txt
		(ii) PATH/ClusterAnalysis/COG/0_0/phyl_multi_fun_table_Hsap_loss_cog.txt
		(iii)PATH/ClusterAnalysis/COG/0_0/phyl_multi_fun_total_counts_Hsap_gain_cog.txt 
		(iv) PATH/ClusterAnalysis/COG/0_0/phyl_multi_fun_total_counts_Hsap_loss_cog.txt 
	e.g.
		Summary analysis of COG cluster annotations, for H. sapiens (the number of phylostrata = 40), c=0.0,  
		
		Results in four output files:		
		(i) phyl_multi_fun_table_Hsap_gain_cog.txt --> summary table, i.e. the number of clusters with an assigned COG (or GO) category for a specific phylostratum (GF gain)
		(ii) phyl_multi_fun_table_Hsap_loss_cog.txt --> summary table, i.e. the number of clusters with an assigned COG (or GO) category for for a specific phylostratum (GF loss)
		
		-------		
		(i)-(ii) for each functional category and phylostratum: the number of clusters to which the functional category was assigned
		E.g. (for categories D and M and phylostrata 1 to 4)	
			category_id	name	1	2	3	4	...
			D	Cell cycle control, cell division, chromosome partitioning	1	1	4	2	...
			M	Cell wall/membrane/envelope biogenesis	8	7	11	9	...

		
		(iii) phyl_multi_fun_total_counts_Hsap_gain_cog.txt --> the total number of GF/clusters (GF gain), to which a functional category was assigned 
		(iv) phyl_multi_fun_total_counts_Hsap_loss_cog.txt --> the total number of GF/clusters (GF loss), to which a functional category was assigned 
		
		-------
		(iii)-(iv) for each functional category: the total number of gene families/clusters (gained/lost), to which a functional category was assigned 
		E.g. (for category D)
			category	total_count
			D	102
			...

	
(FA-3) Pyhlostratum analysis - analysis of function enrichment per phylostratum;
		--> using the hypergeometric test (see phylostratum_function_enrichment_analysis.py; )
	
	> time python3 phylostratum_function_enrichment_analysis.py \
				-i /PATH/ClusterAnalysis_2023/COG/Hsap/0_0 \
				-s H \
				-C \
				-c 0.8 \ 
				-p \
				-t YOUR_PATH/Data/COGs.txt \
				-P YOUR_PATH/Data/PS_names_short \
				-o OUTPUT_FOLDER
	
	Options: 
		-i input folder (the output files from (FA-2)),
		-s species (H, D, A, S for H.sap, D.mel, A.thal, S.cer, but it can be extended to any other focal species),
		-C|G
			-C analyze COG functions (either -C or -G is set),
			-G analyze GO functions (either -C or -G is set),
		-c c-value
		-t input folder for COG or GO functions 
			(if -C, then e.g. YOUR_PATH/Data/COGs.txt;
			if -G, then e.g. YOUR_PATH/Data/go-basic.obo
		-P input folder for phylostrata short names (YOUR_PATH/Data/PS_names_short, see (CA-7))
		[-p] plot results (optional: when the option is included, the functions analysis is plotted in the pdf format; 
					the default is false; it may be very slow in case of GO functions)		

	Note: the hypergeometric test is done for each phylostratum and each GO function/COG category across all c-values 0.0 to 0.8.
			We look for the overexpression of each GO function/COG category in that phylostratum.
	For each GO function or COG category, we compute the following parameters for the hypergeometric test:
		(i) quant: the number of successes in sample  (the total number of assignments of a given category for a given phylostratum)
		(ii) sample: the sample size (the total number of assigned categories for a given phylostratum)
		(iii) hit: the number of successes in population (the total number of assignments of a given category over all phylostrata)
		(iv) total: the population size (the total number of assigned categories over all phylostrata)
		
	over_rep = hypergeom.sf(quant - 1, total, hit, sample)
	under_rep = hypergeom.cdf(quant, total, hit, sample)
	p_value = min(over_rep, under_rep) * 2	
	We further adjust the computed p-value using a correction for multiple tests 
		(here: the Benjamini/Hochberg (non-negative) method) W
	p_adjusted = multi.multipletests(p_values, method="fdr_bh")[1]

----------- MCL ANALYSIS PIPELINE -----------
PREREQUISITES:
	- for MMseqs2 search analysis:
		install program MMseqs2 (https://github.com/soedinglab/MMseqs2) [1]
	
	- for MCL analysis:
		install MCL [6]

(MCL-1) 
	a) Run the mmseqs createdb [1] program to create the "db_all" from "db_all.faa" (if not already created at the step (CA-3))
	> mmseqs createdb db_all.faa db_all

	b) run the mmseqs search [1] on the db_all for each c = 0.0 to 0.8 (the results of the mmseqs search analysis are the input for MCL analysis): 
	> mmseqs search db_all db_all resultDB tmp -e 0.001 -c 0.8 --max-seqs 400 -s 4
		(-c 0.8 can be replaces with any other c-value 0.0 to 0.7)
	
	c) convert the results of the mmseqs search analysis to the input for MCL analysis: 
	> mmseqs convertalis db_all db_all resultDB result_0_8.tab --format-output "query,target,evalue,bits,alnlen,qlen,tlen"
	
(MCL-2)
	a) sorting the mmseqs search results (sorting by query taxId(numerically), subject taxId (numerically), query pid (numerically), subject pid (numerically): 
	e. g.
	--> sort -t '|' -k4,4n -k8,8n -k2,2n -k6,6n result_0_8.tab -o result_0_8_sorted.tab
	
	b) normalize the mmseqs search bit scores using the adapted OrthoFinder's algorithm [7] --> using the "normBitScore.py" script 
	e.g.
	> python3 normBitScore.py -i result_0_8_sorted.tab -o norm_0_8_sorted.tab
	
	The adapted OrthoFinder's algorithm [7]:
	# for each query_taxID # q
	#	for each subject_taxID # h
	# 		read all rows for q, h (since the data are sorted by query_taxID, subject_taxID)
	#		while reading compute for each row: L_q_h (L_q x L_h)
	#		count rows for q, h
	# 		sort by L_q_h
	#		split rows in bins of size = 10^3 according to L_q_h
	#		for each bin
	#			select the top 5% S'(i.e. B_q_h) 
	#			compute a, b using lin. regression (least squares) in log10(B_q_h) = a x log10(L_q_h) + b
	#			normalise each B_q_h --> B_q_h' = B_q_h / (10^b x L_q_h^a)
	#

(MCL-3) Compute taxId_ps_Hsap.txt --> for each taxon in the database find the LCA when to compared to H. sapiens as the focal species 
	Required input:
	- the focal species is set using the -f option; here the focal species is H. sapiens (taxId = 9606), so we use the option "-f 9606"
	- the file "allParents.txt" (see step (CA-5b)), which comprises the list of parents for each taxon in the database
	
	> python3 getAllPS.py -f 9606 -i allParents.txt -p ./Parents -o taxId_ps_Hsap.txt
		
(MCL-4)
	-- Analysis of the mmseqs search output using the neg. logarithm of the E-value as a distance ---
	a) transform the mmseqs search output to abc input from mcxload (done for each c-value; in the example below shown for c=0.8)
	> cut -f 1,2,3 result_0_8.tab > result_0_8.abc

	b) run mcxload [6] to get binary input for mcl (done for each c-value; in the example below shown for c=0.8)
	> mcxload -abc result_0_8.abc --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o result_0_8_bin.mci -write-tab result_0_8_bin.tab --write-binary
	
	c) run MCL [6] for I=2, I=1.5 for each c-value (in the example below shown for c=0.8)
	> mcl result_0_8_bin.mci -I 2 -use-tab result_0_8_bin.tab -o mcl_I_2_0_8.out	
	> mcl result_0_8_bin.mci -I 1.5 -use-tab result_0_8_bin.tab -o mcl_I_1_5_0_8.out
	
	d) parse MCL for I=2, I=1.5 output to get gene families gain/loss over all clusters for each c-value (in the example below shown for c=0.8)
	> python3 parseMCL_2_strict.py -i mcl_I_2_0_8.out.txt -p taxId_ps_Hsap.txt -P 40 \
							-o clu_mcl_I_2_0_8.txt -s summary_mcl_I_2_0_8.txt
	
	> python3 parseMCL_2_strict.py -i mcl_I_1_5_0_8.out.txt -p taxId_ps_Hsap.txt -P 40 \
							-o clu_mcl_I_1_5_0_8.txt -s summary_mcl_I_1_5_0_8.txt
	
(MCL-5)
	-- Analysis of the mmseqs search output using the normalized bit-score as a distance (the output of the step (MCL-2)) ---

	b) run mcxload [6] to get binary input for mcl (done for each c-value; in the example below shown for c=0.8)
	# Please note that the default course of action for mcxload is to use the best value found between a pair of labels. 
	# The bit scores were already normalized in (MCL-2), so I don't use the options --stream-neg-log10 -stream-tf 'ceil(200)' (see (MCL-4))
	# Instead of the option --stream-mirror, I use -ri max, since it uses less memory.

	> mcxload -abc norm_0_8_sorted.tab -ri max  -o norm_0_8_bin.mci -write-tab norm_0_8_bin.tab --write-binary
	
	c) run MCL [6] for I=2, I=1.5 for each c-value (in the example below shown for c=0.8)
	> mcl norm_0_8_bin.mci -I 2 -use-tab norm_0_8_bin.tab -o mcl_norm_I_2_0_8.out	
	> mcl norm_0_8_bin.mci -I 1.5 -use-tab norm_0_8_bin.tab -o mcl_norm_I_1_5_0_8.out

	d) parse MCL for I=2, I=1.5 output to get gene families gain/loss over all clusters for each c-value (in the example below shown for c=0.8)
	> python3 parseMCL_2_strict.py -i mcl_norm_I_2_0_8.out.txt -p taxId_ps_Hsap.txt -P 40 \
							-o clu_mcl_norm_I_2_0_8.txt -s summary_mcl_norm_I_2_0_8.txt
	
	> python3 parseMCL_2_strict.py -i mcl_norm_I_1_5_0_8.out.txt -p taxId_ps_Hsap.txt -P 40 \
							-o clu_mcl_norm_I_1_5_0_8.txt -s summary_mcl_norm_I_1_5_0_8.txt



	


----------- REFERENCES ----------- 

[1] Steinegger M and Soeding J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. 
	Nature Biotechnology, doi: 10.1038/nbt.3988 (2017).

[2] Mosè Manni, Matthew R Berkeley, Mathieu Seppey, Felipe A Simão, Evgeny M Zdobnov, 
	BUSCO Update: Novel and Streamlined Workflows along with Broader and Deeper Phylogenetic 
	Coverage for Scoring of Eukaryotic, Prokaryotic, and Viral Genomes. 
	Molecular Biology and Evolution, Volume 38, Issue 10, October 2021, Pages 4647–4654


[3] eggNOG-mapper v2: functional annotation, orthology assignments, and domain 
    prediction at the metagenomic scale. Carlos P. Cantalapiedra, 
    Ana Hernandez-Plaza, Ivica Letunic, Peer Bork, Jaime Huerta-Cepas. 2021.
    Molecular Biology and Evolution, msab293, https://doi.org/10.1093/molbev/msab293


[4] http://geneontology.org/docs/download-ontology/

[5] http://clovr.org/docs/clusters-of-orthologous-groups-cogs/

[6] van Dongen, Stijn, Graph clustering via a discrete uncoupling process, 
	Siam Journal on Matrix Analysis and Applications 30-1, p121-141, 2008. (https://doi.org/10.1137/040608635)
	
[7] Emms DM, Kelly S. OrthoFinder: solving fundamental biases in whole genome comparisons dramatically 
	improves orthogroup inference accuracy. Genome Biol. 2015 Aug 6;16(1):157. doi: 10.1186/s13059-015-0721-2



