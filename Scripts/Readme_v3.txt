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

***** Data required for the cluster analysis pipeline *****
(1) a database that contains proteins in fasta format (a fasta file for each species in the database)
		- we used the following format:
			- each genome is in a separate file with an extension "faa" and is named taxID.faa; 
				e. g. for H. sapiens (taxID = 9606): 9606.faa
			- each file contains protein sequences for a species
			- header lines in files have the following form; e.g.: 
			>pid|0000000000000012601|tx|1260|	KXA09781
			
			where	
				pid|0000000000000012601|tx|1260| represents a protein name (taxID = 1260; proteinID = 0000000000000012601)
				KXA09781 represents an original protein name
		
		
		- to get the above format for protein sequences (i.e. pid|0000000000000012601|tx|1260|), 
			you may use the script "filter_E.sh" that generates the .faa output from the original fasta files
			- the required input: a path to a folder containing original fasta files and 
				a file containing a list of all taxa/genomes with their taxID-s and names (taxIDNames.txt; see (2))
			> bash filter_E.sh taxIDNames.txt path_to_original_fasta_files path_to_edited_fasta_files > filter_E.log

(2) taxIDNames.txt --> save in YOUR_PATH/Data
		(see an example file in PhyLoss/ExampleData)
		- the list of all species in the database: taxID + species name
		E.g.
			9606	Homo_sapiens

(3) names.txt --> save in YOUR_PATH/Data
		(see an example file in PhyLoss/ExampleData)
		- the list of all nodes in the phylogeny - from the species (leaves in the tree) up to the cellular organism: nodeID + node name
		E.g.
			9606	Homo_sapiens
			131567	Cellular_organisms

(3) nodes.txt --> save in YOUR_PATH/Data
		(see an example file in PhyLoss/ExampleData)
		- the list of all nodes' pairs in the phylogeny: nodeID + parent_nodeID
		E.g. 
			9606	131567
		where		
			9606	Homo_sapiens
			131567	Cellular_organisms

***** Data required for the functional analysis pipeline *****
(4) the clusters computed using MMseqs2 and the input phylogeny data (1) - (3)
	The clusters are copmuted using scripts in "CLUSTER ANALYSIS PIPELINE".

(5) eggnog files --> save in YOUR_PATH/Data
		(see an example file in PhyLoss/ExampleData/chunk_example0.emapper.annotations
		E.g.
		-------------
		# emapper version: emapper-1.0.3-5-g6972f60 emapper DB: 4.5.1
		# time: Thu Jul 23 17:46:52 2020
		#query_name	seed_eggNOG_ortholog	seed_ortholog_evalue	seed_ortholog_score	predicted_gene_name	GO_terms	KEGG_KOs	BiGG_reactions	Annotation_tax_scope	OGs	bestOG|evalue|score	COG cat	eggNOG annot
		pid|0000000000000008821|tx|882|	882.DVU0001	5.8e-253	878.2	DNAA	GO:0000166,GO:0001882,GO:0001883,GO:0003674,GO:0003676,GO:0003677,GO:0003688,GO:0003824,GO:0005488,GO:0005524,GO:0005575,GO:0005618,GO:0005622,GO:0005623,GO:0005737,GO:0005886,GO:0006139,GO:0006152,GO:0006163,GO:0006164,GO:0006172,GO:0006195,GO:0006200,GO:0006461,GO:0006725,GO:0006753,GO:0006793,GO:0006796,GO:0006807,GO:0008150,GO:0008152,GO:0009056,GO:0009058,GO:0009116,GO:0009117,GO:0009119,GO:0009123,GO:0009124,GO:0009125,GO:0009126,GO:0009127,GO:0009128,GO:0009132,GO:0009133,GO:0009135,GO:0009136,GO:0009141,GO:0009143,GO:0009144,GO:0009146,GO:0009150,GO:0009152,GO:0009154,GO:0009156,GO:0009158,GO:0009161,GO:0009163,GO:0009164,GO:0009165,GO:0009166,GO:0009167,GO:0009168,GO:0009169,GO:0009179,GO:0009180,GO:0009185,GO:0009188,GO:0009199,GO:0009203,GO:0009205,GO:0009207,GO:0009259,GO:0009260,GO:0009261,GO:0009987,GO:0016020,GO:0016043,GO:0016311,GO:0016462,GO:0016787,GO:0016817,GO:0016818,GO:0016887,GO:0017076,GO:0017111,GO:0018130,GO:0019438,GO:0019439,GO:0019637,GO:0019693,GO:0022607,GO:0030312,GO:0030554,GO:0032549,GO:0032550,GO:0032553,GO:0032555,GO:0032559,GO:0034641,GO:0034654,GO:0034655,GO:0035639,GO:0036094,GO:0042278,GO:0042451,GO:0042454,GO:0042455,GO:0043167,GO:0043168,GO:0043565,GO:0043933,GO:0044085,GO:0044237,GO:0044238,GO:0044248,GO:0044249,GO:0044270,GO:0044271,GO:0044281,GO:0044424,GO:0044464,GO:0044710,GO:0046031,GO:0046034,GO:0046128,GO:0046129,GO:0046130,GO:0046390,GO:0046434,GO:0046483,GO:0046700,GO:0051259,GO:0051260,GO:0055086,GO:0065003,GO:0070271,GO:0071704,GO:0071822,GO:0071840,GO:0071944,GO:0072521,GO:0072522,GO:0072523,GO:0090407,GO:0097159,GO:1901135,GO:1901136,GO:1901137,GO:1901265,GO:1901292,GO:1901293,GO:1901360,GO:1901361,GO:1901362,GO:1901363,GO:1901564,GO:1901565,GO:1901566,GO:1901575,GO:1901576,GO:1901657,GO:1901658,GO:1901659	K02313		bactNOG[38]	05CI4@bactNOG,0GB44@delNOG,16QDA@proNOG,COG0593@NOG	NA|NA|NA	L	it binds specifically double-stranded DNA at a 9 bp consensus (dnaA box) 5'-TTATC CA A CA A-3'. DnaA binds to ATP and to acidic phospholipids (By similarity)
		-------------
		
		- the precomputed genes functional annotations (COG and GO annotations):
			- extract to the subfolder "PATH_TO_ANALYSIS_FOLDER/Data/eggnogResults"
			- for more information, see "FUNCTIONAL ANALYSIS PIPELINE", step (2)

(6) go-basic.obo --> save in YOUR_PATH/Data
		(see an example file in PhyLoss/ExampleData)
		- the GO-basic functions file

(7) COGs.txt --> save in YOUR_PATH/Data
		(see an example file in PhyLoss/ExampleData)
		- the COG functions file

----------- CLUSTER ANALYSIS PIPELINE -----------

PREREQUISITES:
	- install programs: 
		python 3.6 or later (http://www.python.org/getit/)
	
	- for MMseqs2 cluster analysis:
		install program MMseqs2 (https://github.com/soedinglab/MMseqs2)
	

REQUIREMENTS: 
	- create an analysis folder
	- copy DATA files to the analysis folder
	- copy SCRIPTS files copy to the analysis folder
	
(1) create the subolders "DB", "All_faa" in the analysis folder 
	>> mkdir DB
	>> mkdir All_faa

(2) copy database files to the subfolder "All_faa"

(3) concatenate all faa files in a single faa file "all_db.faa" (in the subolder "DB")	
	>> cat All_faa/*.faa > DB/all_db.faa

(4) running mmseqs cluster (with options -e 0.001, -threads 12, --max-seqs 400, -c 0.0 to 0.8)
	>> bash submit_cluster_c.sge &
	Alternatively: run the script using qsub
	
	Output: 
		- the list of subfolders named "PATH_TO_TSV_FOLDER/results_0_0" to "PATH_TO_TSV_FOLDER/results_0_8", 
			each containing the mmseqs output tsv file "db_clu_all.tsv"
	
(5a) If there are parental nodes with only one leaf, then nodes, names and taxa file should be updated using upd_NodesWith1Child.py:  
	>> python3 upd_NodesWith1Child.py -i /storage/home/mdomazet/ProteinsT/Data/nodes_All_2022.txt \
										-j YOUR_PATH/Data/names.txt \
										-t YOUR_PATH/Data/taxIDNames.txt \
										-o YOUR_PATH/Data/nodes_with1Child.txt \
										-n YOUR_PATH/Data/nodes.txt \ 
										-m YOUR_PATH/Data/names_UPD.txt \
										-s YOUR_PATH/Data/taxIDNames_UPD.txt > out.txt

(5b) creating list of parents for each focal species (in "Parents" subdirectory);
	- parents are parental nodes from the focal species up to the root (i.e. "Cellular organisms")
	>> mkdir YOUR_PATH/Data/Parents
	>> time python3 treeTaxId_all.py -i YOUR_PATH/Data/names_UPD.txt -j YOUR_PATH/Data/nodes_UPD.txt \
									-t YOUR_PATH/Data/taxIDNames_UPD.txt -o YOUR_PATH/Data/Parents/ > log.txt
	
		Output:
			- a list of files in the subfolder "Parents" - a file for each taxID, e.g. "Parents/7227parents.txt" and "allParents.txt"
			- each file in the subfolder "Parents" contains a list od parental nodes starting from a taxID up to the node "Cellular_organisms"
	
	>> cat YOUR_PATH/Data/Parents/*parents.txt > YOUR_PATH/Data/allParents.txt


(6) computing LCA (Last Common Ancestor) matrix for all taxID (n x n; in pur case: 667 x 667); taxon's names are truncated to 10 char-s:
	>> python3 getAllLCA.py -i YOUR_PATH/Data/allParents.txt -o YOUR_PATH/Data/allLCA.txt  

	Output (allLCA.txt):
									880073     526226     526227    1924735	... 
								Caldithrix Gordonia_b Meiothermu Salpingoec 
			880073 (Caldithrix)         -1          2          2     131567 
			526226 (Gordonia_b)          2         -1          2     131567 
			526227 (Meiothermu)          2          2         -1     131567 
			...

(7) computing gene gain/lost for a focal species (given taxID, e.g. 9606 H. sap. and the number of its phylostrata: 39), and a c-value
	E. g. for H. sapiens and c = 0.0 (i.e. the subfolder "PATH_TO_TSV_FOLDER/results_0_0")
	>> python3 getGL.py -f 9606 -p YOUR_PATH/Data/Parents -i PATH_TO_TSV_FOLDER/results_0_0/db_clu_all.tsv 
								-l YOUR_PATH/Data/allLCA.txt -o YOUR_PATH/Data/res_dbAll_H_sap_0_0.txt -s YOUR_PATH/Data/summary_GL_H_sap_0_0.txt > out.txt

	Output: 
		(1) summary output file for each focal species and each c (e.g. summary_GL_H_sap_0_0.txt ):
		-------------	
			ps	GF_Gained	GF_Lost
			1	8393	0
			2	193	0
		...
		-------------
		
		(2) report file containing cluster analysis for clusters involved in the GF gain/loss for a focal species and c
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

	- to retrieve your own eggnog annotations for the files in the database
		install EggNog mapper (https://github.com/eggnogdb/eggnog-mapper)	
		Note: 
			- we applied EggNog mapper to our database 
				- emapper version: emapper-1.0.3-5-g6972f60 emapper DB: 4.5.1
	
REQUIREMENTS: 
	- position in the analysis folder
		>> cd PATH_TO_ANALYSIS_FOLDER
	
	
(1) computing the clusters for functional analysis using getGL_genes.py (output folder: ClusterAnalysis)
		- in the analysis folder create the subolder "ClusterAnalysis", and its subfolders "c_0_0" to "c_0_8"
			>> mkdir YOUR_PATH/ClusterAnalysis/COG
			>> mkdir YOUR_PATH/ClusterAnalysis/GO
			>> mkdir YOUR_PATH/ClusterAnalysis/COG/c_0_0 c_0_1 c_0_2 c_0_3 c_0_4 c_0_5 c_0_6 c_0_7 c_0_8
			>> mkdir YOUR_PATH/ClusterAnalysis/GO/c_0_0 c_0_1 c_0_2 c_0_3 c_0_4 c_0_5 c_0_6 c_0_7 c_0_8
	
	Use getGL_genes.py to get list of clusters for a focal species and a c-value (or getGL_genes.sh for all c-values);
	e.g. for H. sapiens (taxID = 9606, the option -f; c=0.0, the input tsv file: -i PATH_TO_TSV_FOLDER/results_0_0/db_clu_all.tsv)
	>> time python3 getGL_genes.py -f 9606 \
					-p YOUR_PATH/Data/Parents/ \
					-i PATH_TO_TSV_FOLDER/results_0_0/db_clu_all.tsv \
					-l YOUR_PATH/Data/allLCA.txt \
					-g YOUR_PATH/Data/ClusterAnalysis/0_0/res_dbAll_Hsap_0_0_genes.txt \
					-s YOUR_PATH/Data/ClusterAnalysis/0_0/summary_GL_Hsap_0_0.txt > out.txt
	Running time: 15-20 seconds
	
	Use getGL_genes.sh to get list of clusters and their genes for a focal species and all c-values (0.0 to 0.8);
	e.g. for H. sapiens (taxID = 9606, the option -f; c=0.0, the input tsv file: -i PATH_TO_TSV_FOLDER/results_0_0/db_clu_all.tsv)
	>> bash getGL_genes.sh 9606 \
			YOUR_PATH/Data/tsv_new \
			YOUR_PATH/Data/Parents \ 
			YOUR_PATH/Data/allLCA.txt \
			YOUR_PATH/Data/ClusterAnalysis \
			Hsap > geteGL_genes_Hsap.out
	Running time: 9 x 15-20 seconds ~ 3 mins
	
	Output:
		- each cluster contains a list of genes in the cluster
		- in the example below:
			- the first row contains a cluster_ordinal_number (1), a cluster representative/geneId (pid|0000000000000100201|tx|10020|)
				+ the number of genes in the cluster (3), phylostratum_gain (28), phylostratum_loss (30)
			- the following rows are genes in the cluster
	
		E.g.
		1       pid|0000000000000100201|tx|10020|       3       28      30
		pid|0000000000000100201|tx|10020|
		pid|0000000000422542821|tx|42254|
		pid|0000000005133713766|tx|51337|
	

(2) - extract eggnog functional annotations to a folder, e. g. YOUR_PATH/Data/eggnogResults
	E.g. (see an example file in PhyLoss/ExampleData or DATA-->(5))
	- extracted files (emapper.annotations) 
		contain list of genes withe their GO and COG annotations
	
	
(3) Computing COG and GO annotations for clusters
	--- COG --- 
	> bash getGL_genes_funct_cog.sh 9606 YOUR_PATH/Data H_sapiens_9606 > getGL_genes_funct_H_sapiens_9606.out
	(calling getGL_genes_funct.py; 
	the script includes option STRICT, which is initially set to "True"; which means that functions are addigned to clusters more strictly)
	
	Running time: 
		total (c=0.0 to 0.8): 12-13 minutes	
		(for each c: ~1.5 min)
	
	--- GO --- 
	> bash getGL_genes_funct_go.sh 9606 YOUR_PATH/Data H_sapiens_9606 > getGL_genes_funct_H_sapiens_9606_go.out
	(calling getGL_genes_funct.py; 
	the script includes option STRICT, if set to "True" --> stricter assignment of functions to clusters)

	Output:
		- each cluster contains a list of genes in the cluster and, if possible, each gene is assigned a list of COR or GO functions
		- in the example below: 
			- the first row contains a cluster_ordinal_number (2), a cluster representative/geneId (pid|0000000000010020321|tx|10020|)
				+ the number of genes in the cluster (4), phylostratum_gain (20), phylostratum_loss (35)
			- the following rows are genes in the cluster, and, where possible, the assigned list of functions for each gene
			(the gene "pid|0000000000010020321|tx|10020|" is assigned the COG function "T")
	
		E.g.
		2       pid|0000000000010020321|tx|10020|       4       20      35
		pid|0000000000010020321|tx|10020|       T
		pid|0000000000955532907|tx|9555|        T
		pid|0000000000850215858|tx|8502|
		pid|0000000000786810283|tx|7868|

(4) COG anaylsis using parseClusterFunct.py
	E.g. for H. sapiens (-P the number of phylostrata, -S the short name, -i the input data file from the step (3)
			-c path to COG file, -o output folder)
	
	> bash parseClusterFunct.sh 39 Hsap PATH/Data/ClusterAnalysis/COG PATH/Data/COGs.txt
	(calling
		parseClusterFunct.py -P 39 -S Hsap \
				-i YOUR_PATH/Data/ClusterAnalysis/COG/0_0/res_dbAll_H_sapiens_9606_0_0_clusters_funct.txt 
				-c YOUR_PATH/Data/COGs.txt 
				-o YOUR_PATH/Data/ClusterAnalysis/COG/0_0 > cog_out.txt
	
	Running time: 9 x 0.3 seconds ~ 3 seconds
	
(5) pyhlostratum analysis - analysis of function enrichment per phylostratum (including hypergeometric test)
	> time python3 phylostratum_function_enrichment_analysis_v2_new.py \
				-i YOUR_PATH/Data \
				-s H -C -c 0.8 -p
	(Options: -i input folder; -s species; -C COG functions; -G GO functions; -c c-value; -p plot results)
	




	
