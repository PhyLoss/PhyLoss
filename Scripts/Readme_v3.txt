# Author: Mirjana Domazet-Loso
# 2022 - 2023
# Contact: Mirjana.Domazet-Loso@fer.hr
#------------------------------------------
#
# Reference:
# 	Mirjana Domazet-Lošo, Tin Široki, Korina Šimičević, Tomislav Domazet-Lošo, 2023, Macroevolutionary dynamics of gene family gain and loss across multicellular eukaryotic lineages
# 
# ------------------------------------------
# Code contributors:
#	Mirjana Domazet-Lošo and Tin Široki
#	University of Zagreb Faculty of Electrical Engineering and Computing, Zagreb, Croatia

-- PREREQUISITES --
	- installed programs: 
		python 3

----------- DATA -----------
# create a folder "Data" in your working directory, e.g.
	> cd YOUR_PATH
	> mkdir Data
	> mkdir Data/DB	

***** Data required for the cluster analysis pipeline *****
(1) a database that contains protein sequences in fasta format, i.e. a list of fasta files (a single file for each species) --> save in YOUR_PATH/Data/DB
		- we used the following format for fasta files (required in subsequent scripts):
			- a separate file for each species named taxID.faa, where "taxID" is repaced with a true taxID for a species 
				e. g. for H. sapiens (taxID = 9606): 9606.faa
			- each .faa file contains protein sequences for the corresponding species
			- header lines in our .faa files have the following format, e.g.: 
			>pid|0000000000000012601|tx|1260|	KXA09781
			
			where	
				pid|0000000000000012601|tx|1260| represents a protein name (taxID/tx = 1260; proteinID/pid = 0000000000000012601)
				KXA09781 represents an original protein name		
		
		Note:
		- to get the above format for protein sequences (i.e. pid|0000000000000012601|tx|1260|), 
			you may use the script "filter_E.sh" that generates the .faa files from the original fasta files and 
			a file containg all genomes' taxID-s e.g. taxIDNames.txt (see (2))
			> bash filter_E.sh taxIDNames.txt path_to_original_fasta_folder YOUR_PATH/Data/DB > filter_E.log

(2) taxIDNames.txt --> save in YOUR_PATH/Data
		(see an example file in PhyLoss/ExampleData)
		- the list of all species in the database in the format: 
			taxID	species_name
		E.g.
			9606	Homo_sapiens

(3a) names.txt --> save in YOUR_PATH/Data
		(see an example file in PhyLoss/ExampleData)
		- the list of all nodes in the input phylogeny: from the species/leaves in the tree up to the root/cellular organism;
		- the required format: 
			nodeID	node_name
		E.g.
			9606	Homo_sapiens
			131567	Cellular_organisms

(3b) nodes.txt --> save in YOUR_PATH/Data
		(see an example file in PhyLoss/ExampleData)
		- the list of all nodes' pairs in the input phylogeny in the format:
			nodeID	parent_nodeID
		E.g. 
			9606	131567
		where		
			9606	Homo_sapiens
			131567	Cellular_organisms

***** Data required for the functional analysis pipeline *****
(4) the clusters computed using MMseqs2 and the input phylogeny data (1) - (3)
	The clusters are copmuted using scripts in "CLUSTER ANALYSIS PIPELINE".

(5) eggnog files --> save in YOUR_PATH/Data
		(see an example file in PhyLoss/ExampleData/chunk_sample.emapper.annotations)
		- the precomputed genes functional annotations (COG and GO annotations):
			- extract to the subfolder "PATH_TO_ANALYSIS_FOLDER/Data/eggnogResults"
			- for more information, see "FUNCTIONAL ANALYSIS PIPELINE", step (2)

(6) go-basic.obo --> save in YOUR_PATH/Data
		- the GO-basic functions file
		(see http://geneontology.org/docs/download-ontology/)

(7) COGs.txt --> save in YOUR_PATH/Data
		- the COG functional categories file
		(see http://clovr.org/docs/clusters-of-orthologous-groups-cogs/)

----------- CLUSTER ANALYSIS PIPELINE -----------

PREREQUISITES:
	- install programs: 
		python 3.6 or later (http://www.python.org/getit/)
	
	- for MMseqs2 cluster analysis:
		install program MMseqs2 (https://github.com/soedinglab/MMseqs2)

REQUIREMENTS: 
	- list of SCRIPTS files
		- getAllLCA.py
		- getAllPS.py
		- getGL.py
		- submit_cluster_c.sge
		- treeTaxId_all.py
		- upd_NodesWith1Child.py		

(1) concatenate all faa files in a single faa file "all_db.faa" (in the subolder "DB")	
	> cat YOUR_PATH/Data/DB/*.faa > YOUR_PATH/Data/DB/all_db.faa

(2) position to the folder "YOUR_PATH/Data"
	> cd YOUR_PATH/Data

(3) running the mmseqs cluster program (with options -e 0.001, -threads 12, --max-seqs 400, -c 0.0 to 0.8)
	> bash submit_cluster_c.sge &
	Alternatively: run the script using qsub
	
	Output: 
		- the list of subfolders named "YOUR_PATH/Data/results_0_0" to "YOUR_PATH/Data/results_0_8"
			each containing the mmseqs output tsv file "db_clu_all.tsv"
	
(4) If there are parental nodes with only one leaf, then nodes, names and taxa file should be updated using upd_NodesWith1Child.py:  
	> python3 upd_NodesWith1Child.py -i YOUR_PATH/Data/nodes.txt \
										-j YOUR_PATH/Data/names.txt \
										-t YOUR_PATH/Data/taxIDNames.txt \
										-o YOUR_PATH/Data/nodes_with1Child.txt \
										-n YOUR_PATH/Data/nodes_UPD.txt \ 
										-m YOUR_PATH/Data/names_UPD.txt \
										-s YOUR_PATH/Data/taxIDNames_UPD.txt > out.txt

(5) creating a list of parents for each focal species (in "Parents" subdirectory);
	- parents are parental nodes from the focal species up to the root (i.e. "Cellular organisms")
	> mkdir YOUR_PATH/Data/Parents
	> time python3 treeTaxId_all.py -i YOUR_PATH/Data/names_UPD.txt -j YOUR_PATH/Data/nodes_UPD.txt \
									-t YOUR_PATH/Data/taxIDNames_UPD.txt -o YOUR_PATH/Data/Parents/ > log.txt
	
		Output:
			- a list of files in the subfolder "Parents": 
				a file for each taxID, e.g. "Parents/7227parents.txt" and "allParents.txt";
				each file in the subfolder "Parents" contains a list od parental nodes starting from a taxID up to the node "Cellular_organisms"
	
	>> cat YOUR_PATH/Data/Parents/*parents.txt > YOUR_PATH/Data/allParents.txt


(6) computing LCA (Last Common Ancestor) matrix for all taxID (n x n; in pur case: 667 x 667); taxon's names are truncated to 10 char-s:
	> python3 getAllLCA.py -i YOUR_PATH/Data/allParents.txt -o YOUR_PATH/Data/allLCA.txt  

	E.g. output (allLCA.txt):
									880073     526226     526227    1924735	... 
								Caldithrix Gordonia_b Meiothermu Salpingoec 
			880073 (Caldithrix)         -1          2          2     131567 
			526226 (Gordonia_b)          2         -1          2     131567 
			526227 (Meiothermu)          2          2         -1     131567 
			...

(7) computing gene families' gain/lost for a focal species/taxID 
	(e.g. if taxID = 9606, then compute the gene families' gain/loss for H. sap.)
	
	E. g. for H. sapiens and c = 0.0 (i.e. the input subfolder "YOUR_PATH/Data/results_0_0")
	> python3 getGL.py -f 9606 -p YOUR_PATH/Data/Parents -i YOUR_PATH/Data/results_0_0/db_clu_all.tsv \
								-l YOUR_PATH/Data/allLCA.txt -o YOUR_PATH/Data/res_dbAll_H_sap_0_0.txt -s YOUR_PATH/Data/summary_GL_H_sap_0_0.txt

	Output: 
		(1) a summary output file for each focal species and each c (e.g. summary_GL_H_sap_0_0.txt ):
		-------------	
			ps	GF_Gained	GF_Lost
			1	8393	0
			2	193	0
		...
		-------------
		
		(2) a report file containing cluster analysis for clusters involved in gene families gain/loss for a focal species and c
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
			
	Running time: 
		15-20 seconds per focal species and c

----------- FUNCTIONAL ANALYSIS PIPELINE -----------
PREREQUISITES:
	- data computed using "Cluster analysis pipeline" (see above)

	- to retrieve your own eggnog annotations for the files in your database
		install EggNog mapper (https://github.com/eggnogdb/eggnog-mapper)	
		Note: 
			- we applied EggNog mapper to our database 
				- emapper version: emapper-1.0.3-5-g6972f60 emapper DB: 4.5.1
	
REQUIREMENTS: 
	- create your analysis folder
		> mkdir YOUR_PATH/Data/ClusterAnalysis	
	
	- create the output subfolders (for the COG and GO analysis; for each c-value a subfolder "c_0_0" to "c_0_8")
		> mkdir YOUR_PATH/Data/ClusterAnalysis/COG
		> mkdir YOUR_PATH/Data/ClusterAnalysis/GO
		> mkdir YOUR_PATH/Data/ClusterAnalysis/COG/c_0_0 c_0_1 c_0_2 c_0_3 c_0_4 c_0_5 c_0_6 c_0_7 c_0_8
		> mkdir YOUR_PATH/Data/ClusterAnalysis/GO/c_0_0 c_0_1 c_0_2 c_0_3 c_0_4 c_0_5 c_0_6 c_0_7 c_0_8
	
(1) computing clusters for the functional analysis using getGL_genes.py
	Note:
	the common clusters' data required for both the COG and GO functional analysis, will be stored in the subfolder "YOUR_PATH/Data/ClusterAnalysis/COG" 
	
	Use getGL_genes.py to get list of clusters for a focal species and a c-value (or getGL_genes.sh for all c-values);
	e.g. for H. sapiens (taxID = 9606 --> the option -f 9606; 
							c = 0.0 --> the input tsv file: -i YOUR_PATH/Data/results_0_0/db_clu_all.tsv)
	
	> python3 getGL_genes.py -f 9606 \
					-p YOUR_PATH/Data/Parents/ \
					-i YOUR_PATH/Data/results_0_0/db_clu_all.tsv \
					-l YOUR_PATH/Data/allLCA.txt \
					-g YOUR_PATH/Data/ClusterAnalysis/COG/0_0/res_dbAll_Hsap_0_0_genes.txt \
					-s YOUR_PATH/Data/ClusterAnalysis/COG/0_0/summary_GL_Hsap_0_0.txt > out.txt
	Running time: 15-20 seconds
	
	Use getGL_genes.sh to get list of clusters and their genes for a focal species and all c-values (0.0 to 0.8);
	e.g. for H. sapiens (taxID = 9606 --> the option -f 9606; 
							c = 0.0 --> the input tsv file: -i YOUR_PATH/Data/results_0_0/db_clu_all.tsv)
	> bash getGL_genes.sh 9606 \
			YOUR_PATH/Data \
			YOUR_PATH/Data/Parents \ 
			YOUR_PATH/Data/allLCA.txt \
			YOUR_PATH/Data/ClusterAnalysis/COG \
			Hsap > geteGL_genes_Hsap.out
	Running time: 9 x 15-20 seconds ~ 3 mins
	
	Output:
		- each cluster contains a list of its genes
		- in the example below:
			- the first row contains a cluster_ordinal_number (1), a cluster representative/geneId (pid|0000000000000100201|tx|10020|),
				the number of genes in the cluster (3), phylostratum_gain (28), phylostratum_loss (30) (from the perspective of the focal species)
			- the following rows are genes in the cluster
	
		E.g.
		1       pid|0000000000000100201|tx|10020|       3       28      30
		pid|0000000000000100201|tx|10020|
		pid|0000000000422542821|tx|42254|
		pid|0000000005133713766|tx|51337|
	

(2) - extract eggnog functional annotations to a folder, e. g. YOUR_PATH/Data/eggnogResults
	E.g. (see an example file "chunk_sample.emapper.annotations" in PhyLoss/ExampleData)
	
	- extracted files (emapper.annotations) contain list of genes withe their GO and COG annotations
	
(3) Computing COG and GO annotations for clusters
	--- COG --- 
	> bash getGL_genes_funct_cog.sh 9606 YOUR_PATH/Data H_sapiens_9606 > getGL_genes_funct_H_sapiens_9606.out
	(calling getGL_genes_funct.py; 
	the script includes option STRICT, if set to "True" (default: "False") --> stricter assignment of functions to clusters; at least 50% of the cluster members have to be assigned a function, for a function to be assigned to the whole cluster)
	
	Running time: 
		total (c=0.0 to 0.8): 12-13 minutes	
		(for each c: ~1.5 min)
	
	--- GO --- 
	> bash getGL_genes_funct_go.sh 9606 YOUR_PATH/Data H_sapiens_9606 > getGL_genes_funct_H_sapiens_9606_go.out
	(calling getGL_genes_funct.py; 
	the script includes option STRICT, if set to "True" (default: "False") --> stricter assignment of functions to clusters; at least 50% of the cluster members have to be assigned a function, for a function to be assigned to the whole cluster)

	Output:
		- each cluster contains a list of genes in the cluster and, if possible, for each gene a list of COG or GO functions (if found using the eggnog mapper)
		- in the example below: 
			- the first row contains a cluster_ordinal_number (2), a cluster representative/geneId (pid|0000000000010020321|tx|10020|),
				the number of genes in the cluster (4), phylostratum_gain (20), phylostratum_loss (35) (from the focal species perspective)
			- the following rows correspond to the genes in the cluster, and, where possible, a list of assigned COG o GO functions for each gene
			(in this example: the gene "pid|0000000000010020321|tx|10020|" is assigned the COG function "T")
	
		E.g.
		2       pid|0000000000010020321|tx|10020|       4       20      35
		pid|0000000000010020321|tx|10020|       T
		pid|0000000000955532907|tx|9555|        T
		pid|0000000000850215858|tx|8502|
		pid|0000000000786810283|tx|7868|

(4) COG anaylsis using parseClusterFunct.py
	Output: 
		two files for gene families gain analysis (in the form (i)-(ii)) and two output files for gene families lost analysis (in the form (i)-(ii))
	
		(i) the count of gene families/clusters (gained/lost) assigned a specific functional category
		E.g. (for category D)
			category	total_count
			D	102
			...
		
		(ii) the count of clusters with assigned functional category (gained/lost) per phylostratum
		E.g. (for categories D and M and phylostrata 1 to 4)	
			category_id	name	1	2	3	4	...
			D	Cell cycle control, cell division, chromosome partitioning	1	1	4	2	...
			M	Cell wall/membrane/envelope biogenesis	8	7	11	9	...

	
	E.g. for H. sapiens: 
	> bash parseClusterFunct.sh 39 Hsap PATH/Data/ClusterAnalysis/COG PATH/Data/COGs.txt
		
	Note: the script parseClusterFunct.sh calls parseClusterFunct.py for c-values 0.0 to 0.8, 
			but you can call the script parseClusterFunct.py directly for a specific c-value data set, e.g. for c=0.0:
			
			> parseClusterFunct.py -P 39 -S Hsap \
				-i YOUR_PATH/Data/ClusterAnalysis/COG/0_0/res_dbAll_H_sapiens_9606_0_0_clusters_funct.txt 
				-c YOUR_PATH/Data/COGs.txt 
				-o YOUR_PATH/Data/ClusterAnalysis/COG/0_0
	
		Options: -P the number of phylostrata, 
				-S the short name of the focal species, 
				-i the input data file from the step (3)
				-c path to COG file, 
				-o output folder
	Running time (for all c-values): 9 x 0.3 seconds ~ 3 seconds
	
(5) pyhlostratum analysis - analysis of function enrichment per phylostratum (including hypergeometric test)
	> time python3 phylostratum_function_enrichment_analysis_v2_new.py \
				-i YOUR_PATH/Data \
				-s H -C -c 0.8 -p
	Options: -i input folder,
				-s species,
				-C|G
					-C analyze COG functions (either -C or -G is set),
					-G analyze GO functions (either -C or -G is set),
				-c c-value
				[-p] plot results (optional: when the option is included, the functions analysis is plotted in the pdf format; 
									the default is false)
