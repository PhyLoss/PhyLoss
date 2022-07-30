# Author: Mirjana Domazet-Loso
# July, 2022
# Contact: Mirjana.Domazet-Loso@fer.hr
#------------------------------------------
#
# Reference:
# 	Mirjana Domazet-Lošo, Tin Široki, Tomislav Domazet-Lošo, 2022, Macroevolutionary dynamics of gene family gain and loss across multicellular eukaryotic lineages
# 
# ------------------------------------------
# Code contributors:
#	Mirjana Domazet-Lošo and Tin Široki
#	University of Zagreb Faculty of Electrical Engineering and Computing, Zagreb, Croatia


----------- DATA -----------

***** Data required for the cluster analysis pipeline *****
(1) db.zip
		- the database contains 667 files in fasta format (a fasta file for each species in the database)
		- each file in the database has extension "ff" and is named taxID.ff; e. g. for H. sapiens (taxID = 9606): 9606.ff
		- each file contains protein sequences for a species
		- header lines in files have the following form; e.g.: 
			>pgi|0000000000000012601|ti|1260|pi|0|	KXA09781
		where	
			pgi|0000000000000012601|ti|1260|pi|0| represents a protein name (taxID = 1260; proteinID = 0000000000000012601)
			KXA09781 represents an original protein name

(2) taxIDNames_All_2022.txt
		- the list of all 667 species in the database: taxID + species name
		E.g.
			9606	Homo_sapiens

(3) names_All_2022.txt 
		- the list of all nodes in the phylogeny - from the species (leaves in the tree) up to the cellular organism: nodeID + node name
		E.g.
			9606	Homo_sapiens
			131567	Cellular_organisms

(3) nodes_All_2022.txt
		- the list of all nodes' pairs in the phylogeny: nodeID + parent_nodeID
		E.g. 
			1148	2
		where		
			1148 --> Synechocystis_sp_pcc_6803
			2 --> Bacteria

***** Data required for the functional analysis pipeline *****
(4) clusters.tar.gz
		- the clusters computed using MMseqs2 and the phylogeny data  
			- each extracted file corresponds to a focal species (H. sapiens, D. melanogaster, S.cerevisiae, or A. thaliana) and a c-value (c=0.0 to 0.8)
			- the files are extracted to the folder "Clusters_2022" and its subfolders "c_0_0" to "c_0_8" (a subfolder for each c-value)
			(e. g. the subfolder "PATH_TO_ANALYSIS_FOLDER/Data/Clusters_2022/c_0_0" contains the following files: Athal_2022_clusters.txt, Dmel_2022_clusters.txt, Hsap_2022_clusters.txt, and Scer_2022_clusters.txt)
		- for the more information on sample data, see "FUNCTIONAL ANALYSIS PIPELINE", step (1)
		
		NOTE:
		Alterantively, you can compute yur own clusters using scripts in "CLUSTER ANALYSIS PIPELINE"

(5) eggnog.tar.gz
	- the precomputed genes functional annotations (COG and GO annotations):
		- extract to the subfolder "PATH_TO_ANALYSIS_FOLDER/Data/eggnogResults"
		- for the more information on sample data, see "FUNCTIONAL ANALYSIS PIPELINE", step (2)


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
		- computeParents.sh
		- getAllLCA.awk
		- getAllPS.py
		- getGL_Athal.awk
		- getGL_Dmel.awk
		- getGL_HS.awk
		- getGL_Scer.awk
		- numCluster_Athal.sge
		- numCluster_Dmel.sge
		- numCluster_HS.sge
		- numClusters_Athal.awk
		- numCluster_Scer.sge
		- numClusters_Dmel.awk
		- numClusters_HS.awk
		- numClusters.py
		- numClusters_Scer.awk
		- submit_cluster_c.sge
		- treeTaxId.awk		

(1) create the subolders "DB", "All_ff" in the analysis folder 
	>> mkdir DB
	>> mkdir All_ff

(2) unzip db.zip to the subfolder "All_ff"

(3) concatenate all ff files in a single ff file "all_db.ff" (in the subolder "DB")	
	>> cat All_ff/*.ff > DB/all_db.ff

(4) running mmseqs search (with options -e 0.001, -threads 12, --max-seqs 400, -c 0.0 to 0.8)
	>> bash submit_cluster_c.sge &
	Alternatively: run the script using qsub
	
	Output: 
		- the list of subfolders named "results_0_0" to "results_0_8", 
			each containing the mmseqs output tsv file "db_clu_all.tsv"
	
(5) creating list of parents for each focal species (in "Parents" subdirectory) 
	>> bash computeParents.sh nodes_All_2022.txt Parents > computeParents.log
	>> cat Parents/*parents.txt > allParents.txt 
	
	NOTE: 
		computeParents.sh calls treeTaxId.awk
	Output:
		- the list of files in the subfolder "Parents" - a file for each taxID, e.g. "Parents/7227parents.txt" and "allParents.txt"
		- each file in the subfolder "Parents" contains a list od parental nodes starting from a taxID up to the node "Cellular_organisms"

(6) computing LCA (Last Common Ancestor) matrix for all taxID (667 x 667); taxon's names are truncated to 10 char-s:
	>> awk -f getAllLCA.awk allParents.txt > allLCA.txt

	Output (allLCA.txt):
									880073     526226     526227    1924735	... 
								Caldithrix Gordonia_b Meiothermu Salpingoec 
			880073 (Caldithrix)         -1          2          2     131567 
			526226 (Gordonia_b)          2         -1          2     131567 
			526227 (Meiothermu)          2          2         -1     131567 
			...

(7) for each focal species (H. sapiens, D. melanogaster, S.cerevisiae, A. thaliana): 
		from the perspective of the focal species compute phylostrata for each species (taxID) in the database
	
	>> python3 getAllPS.py -f 9606 -i allParents.txt -p ./Parents -o taxId_ps_Hsap.txt
	>> python3 getAllPS.py -f 7227 -i allParents.txt -p ./Parents -o taxId_ps_Dmel.txt
	>> python3 getAllPS.py -f 559292 -i allParents.txt -p ./Parents -o taxId_ps_Scer.txt
	>> python3 getAllPS.py -f 3702 -i allParents.txt -p ./Parents -o taxId_ps_Athal.txt

(9) For each c value (0.0 to 0.8):
		call numClusters.py to process mmseqs output (tsv file) for the c-value 
			and construct separate auxiliary output for each focal species
	
	Example for c = 0.8:
		- the output folder (here "./results_0_8/") should have subfolders for each focal species
		- these subfolders should be named "Hsap", "Dmel", "Scer", "Athal" and can be created using mkdir:
		>> mkdir -p ./results_0_8/Hsap Dmel Scer Athal
		>> python3 numClusters.py -i names_All_2022.txt -j . -t ./results_0_8/db_clu_all.tsv -o ./results_0_8/
		
	Output: 
		an output file for each focal species (e. g. for S. cerevisiae: "Scer/dbAllPlus_clu_0_8_Scer_out.txt")
		-------------	
			(1) Cluster representative: pgi|0000000001002011297|ti|10020|pi|0|      Number of members in the cluster: 2
					Cluster members:
					(1) pgi|0000000001002011297|ti|10020|pi|0|
					Dipodomys_ordii 9
					(2) pgi|0000000000010090805|ti|10090|pi|0|
					Mus_musculus    9
			...
		-------------
			
			where
				the protein pgi|0000000000010090805|ti|10090|pi|0| belongs to Mus_musculus (taxID = 10090), 
					which is in the phylostratum 9 from the perspective of the focal species S. cerevisiae
	
	Running time: 
		several minutes for each numClusters.py call
	
(10) form allLCA_n_667.txt
	>> echo "667" >> allLCA_n_667.txt
	>> cat allLCA.txt >> allLCA_n_667.txt

(11) For each focal species: (H. sapiens, D. melanogaster, S. cerevisiae, A. thaliana) (steps 11a - 11d)
		For each c-value 0.0 to 0.8: 
			compute gene families' gain/loss (based on mmseqs results stored in subfolders "./results_0_0/", .., "./results_0_8/")
	
	Note: 
		- each sge script uses an accompanying awk scipt (getGL_HS.awk, getGL_Dmel.awk, getGL_Scer.awk, getGL_Athal.awk)
		- in each script replace "PATH_TO_ANALYSIS_FOLDER" with the path to the analysis folder
	>> bash numCluster_HS.sge & 
	>> bash numCluster_Athal.sge &
	>> bash numCluster_Dmel.sge &
	>> bash numCluster_Scer.sge &

	Output files (for each focal species and each c-value): gene cluster gain/loss, and an auxiliary results file containing phylostratum for each gene in each cluster
		e. g. "./results_0_0/Athal/gain_loss_0_0_Athal.txt" (for A. thaliana and c = 0.0) 
			and "./results_0_0/Athal/res_dbAllPlus_clu_0_0_Athal_out_EXT.txt"
		
	Running time of each sge script: 
		for each c-value: 5-7 minutes
		total for all c-values: around 60 minutes

	
----------- FUNCTIONAL ANALYSIS PIPELINE -----------
PREREQUISITES:
	- install programs: 
		python 3
	- data computed using "Cluster analysis pipeline" (Data items 4 - 7) (see above)

	- to retrieve your own eggnog annotations for the files in the database
		install EggNog mapper (https://github.com/eggnogdb/eggnog-mapper)	
		Note: 
			- we applied EggNog mapper to our database 
				- emapper version: emapper-1.0.3-5-g6972f60 emapper DB: 4.5.1
	
REQUIREMENTS: 
	- position in the analysis folder
	- in the analysis folder (replace PATH_TO_ANALYSIS_FOLDER with the folder name) create the folders for the functional analysis reults:
		>> mkdir PATH_TO_ANALYSIS_FOLDER/results_2022
		>> mkdir PATH_TO_ANALYSIS_FOLDER/results_2022/cluster_hyper_cog cluster_hyper_go cluster_hyper_cog_summary_v2 cluster_hyper_go_summary_v2 phylostratum_hyper_cog_v2 phylostratum_hyper_go_v2	
		
	- copy SCRIPTS files copy to the analysis folder
		- clu_funct_analysis_v2.sh
		- clu_funct_summary_v2.sh
		- cluster_function_enrichment_analysis_mdl_v4.py
		- cluster_function_enrichment_summary_cog_v2.py
		- cluster_function_enrichment_summary_go_v4.py
		- parseClusters_for_FA_NEW.py
		- parseClusters_NEW.sh
		- phylostratum_function_enrichment_analysis_v2.py
		- phylostratum_function_enrichment_analysis_v3.py
		- phylostr_overall_go.sh

(1) a) Precomputed clusters ready for functional analysis are stored in clusters.tar.gz
	- to extract the archive file to the folder "Clusters_2022"
		>> tar -xzvf clusters.tar.gz -C Clusters_2022

	OR, ALTERNATIVELY: 
	b) If you want to compute the clusters using parseClusters_for_FA_NEW.py (output folder: Clusters_2022)
		- in the analysis folder create the subolder "Clusters_2022", and its subfolders "c_0_0" to "c_0_8"
			>> mkdir ./Clusters_2022
			>> mkdir ./Clusters_2022/c_0_0 c_0_1 c_0_2 c_0_3 c_0_4 c_0_5 c_0_6 c_0_7 c_0_8
	Input: path_to_mmseqs_results_root_folder, a focal species short name, output folder, the number of phylostrata for the focal species
		>> nohup bash parseClusters_NEW.sh path_to_mmseqs_results_root_folder Hsap ./Clusters_2022 39 > parseClusters_Hsap.out &
		>> nohup bash parseClusters_NEW.sh path_to_mmseqs_results_root_folder Dmel ./Clusters_2022 31 > parseClusters_Dmel.out &
		>> nohup bash parseClusters_NEW.sh path_to_mmseqs_results_root_folder Athal ./Clusters_2022 24 > parseClusters_Athal.out &
		>> nohup bash parseClusters_NEW.sh path_to_mmseqs_results_root_folder Scer ./Clusters_2022 28 > parseClusters_Scer.out &
	
	Output sample (e.g. from file "Clusters_2022/c_0_0/Athal_2022_clusters.txt"):
	-------------
		(19)	pgi|0000000010525853528|ti|1052585|pi|0|	1052585	11	14	3
		pgi|0000000010525853528|ti|1052585|pi|0|
		pgi|0000000001294265217|ti|1294265|pi|0|
		pgi|0000000012942653678|ti|1294265|pi|0|
	-------------
	where
		- each cluster contains a list of genes in the cluster
		- the first row contains a geneID (a cluster representative) (pgi|0000000010525853528|ti|1052585|pi|0|), 
			its taxID (1052585), phylostratum_gain (11), phylostratum_loss (14), the number of genes in the cluster (3)

(2) To use precomputed genes functional annotations:
	- extract eggnog.tar.gz to the subfolder "eggnogResults"
		>> mkdir PATH_TO_ANALYSIS_FOLDER/Data/eggnogResults
		>> tar -xzvf eggnog.tar.gz -C eggnogResults
	
	- extracted files (chunk_0.emapper.annotations to chunk_150.emapper.annotations) 
		contain list of genes withe their GO and COG annotations
	
	Note:
		- the data was generated using emapper-1.0.3-5-g6972f60 (emapper DB: 4.5.1) on our database
		- emapper was run with option -m diamond

	
(3) Computing COG and GO annotations for clusters
	IMPORTANT NOTE:
	# In script "cluster_function_enrichment_analysis_v4.py" (called from clu_funct_analysis_v2.sh) change "PATH_TO_THE_ANALYSIS_FOLDER" to a correct path!
    # In script "cluster_function_enrichment_analysis_v4.py" CATEGORIES_IN_PATH = "PATH_TO_THE_ANALYSIS_FOLDER/Data/eggnogResults"
		set to a correct path where eggNog results files (functional annotations) are stored

	Input - an eggNog results file (a sample of our eggNog results file "chunk_0.emapper.annotations")
	-------------
	# emapper version: emapper-1.0.3-5-g6972f60 emapper DB: 4.5.1
	# time: Thu Jul 23 17:46:52 2020
	#query_name	seed_eggNOG_ortholog	seed_ortholog_evalue	seed_ortholog_score	predicted_gene_name	GO_terms	KEGG_KOs	BiGG_reactions	Annotation_tax_scope	OGs	bestOG|evalue|score	COG cat	eggNOG annot
	pgi|0000000000000008821|ti|882|pi|0|	882.DVU0001	5.8e-253	878.2	DNAA	GO:0000166,GO:0001882,GO:0001883,GO:0003674,GO:0003676,GO:0003677,GO:0003688,GO:0003824,GO:0005488,GO:0005524,GO:0005575,GO:0005618,GO:0005622,GO:0005623,GO:0005737,GO:0005886,GO:0006139,GO:0006152,GO:0006163,GO:0006164,GO:0006172,GO:0006195,GO:0006200,GO:0006461,GO:0006725,GO:0006753,GO:0006793,GO:0006796,GO:0006807,GO:0008150,GO:0008152,GO:0009056,GO:0009058,GO:0009116,GO:0009117,GO:0009119,GO:0009123,GO:0009124,GO:0009125,GO:0009126,GO:0009127,GO:0009128,GO:0009132,GO:0009133,GO:0009135,GO:0009136,GO:0009141,GO:0009143,GO:0009144,GO:0009146,GO:0009150,GO:0009152,GO:0009154,GO:0009156,GO:0009158,GO:0009161,GO:0009163,GO:0009164,GO:0009165,GO:0009166,GO:0009167,GO:0009168,GO:0009169,GO:0009179,GO:0009180,GO:0009185,GO:0009188,GO:0009199,GO:0009203,GO:0009205,GO:0009207,GO:0009259,GO:0009260,GO:0009261,GO:0009987,GO:0016020,GO:0016043,GO:0016311,GO:0016462,GO:0016787,GO:0016817,GO:0016818,GO:0016887,GO:0017076,GO:0017111,GO:0018130,GO:0019438,GO:0019439,GO:0019637,GO:0019693,GO:0022607,GO:0030312,GO:0030554,GO:0032549,GO:0032550,GO:0032553,GO:0032555,GO:0032559,GO:0034641,GO:0034654,GO:0034655,GO:0035639,GO:0036094,GO:0042278,GO:0042451,GO:0042454,GO:0042455,GO:0043167,GO:0043168,GO:0043565,GO:0043933,GO:0044085,GO:0044237,GO:0044238,GO:0044248,GO:0044249,GO:0044270,GO:0044271,GO:0044281,GO:0044424,GO:0044464,GO:0044710,GO:0046031,GO:0046034,GO:0046128,GO:0046129,GO:0046130,GO:0046390,GO:0046434,GO:0046483,GO:0046700,GO:0051259,GO:0051260,GO:0055086,GO:0065003,GO:0070271,GO:0071704,GO:0071822,GO:0071840,GO:0071944,GO:0072521,GO:0072522,GO:0072523,GO:0090407,GO:0097159,GO:1901135,GO:1901136,GO:1901137,GO:1901265,GO:1901292,GO:1901293,GO:1901360,GO:1901361,GO:1901362,GO:1901363,GO:1901564,GO:1901565,GO:1901566,GO:1901575,GO:1901576,GO:1901657,GO:1901658,GO:1901659	K02313		bactNOG[38]	05CI4@bactNOG,0GB44@delNOG,16QDA@proNOG,COG0593@NOG	NA|NA|NA	L	it binds specifically double-stranded DNA at a 9 bp consensus (dnaA box) 5'-TTATC CA A CA A-3'. DnaA binds to ATP and to acidic phospholipids (By similarity)
	pgi|0000000000000008822|ti|882|pi|0|	882.DVU0002	3.3e-220	769.2	DNAN	GO:0003674,GO:0003824,GO:0003887,GO:0005575,GO:0005576,GO:0005618,GO:0005622,GO:0005623,GO:0005737,GO:0005829,GO:0006139,GO:0006259,GO:0006260,GO:0006261,GO:0006271,GO:0006725,GO:0006807,GO:0006950,GO:0006974,GO:0008150,GO:0008152,GO:0009058,GO:0009059,GO:0009360,GO:0009987,GO:0016740,GO:0016772,GO:0016779,GO:0022616,GO:0030312,GO:0032991,GO:0033554,GO:0034061,GO:0034641,GO:0034645,GO:0042575,GO:0043170,GO:0043234,GO:0044237,GO:0044238,GO:0044249,GO:0044260,GO:0044424,GO:0044444,GO:0044464,GO:0044699,GO:0044763,GO:0046483,GO:0050896,GO:0051716,GO:0071704,GO:0071944,GO:0090304,GO:1901360,GO:1901576	K02338		bactNOG[38]	05CZ8@bactNOG,0GATM@delNOG,0HDN0@dproNOG,16QFW@proNOG,COG0592@NOG	NA|NA|NA	L	DNA polymerase III is a complex, multichain enzyme responsible for most of the replicative synthesis in bacteria. This DNA polymerase also exhibits 3' to 5' exonuclease activity. The beta chain is required for initiation of replication once it is clamped onto DNA, it slides freely (bidirectional and ATP- independent) along duplex DNA (By similarity)
	pgi|0000000000000008823|ti|882|pi|0|	883.DvMF_2083	0.0	1336.6	GYRB	GO:0000166,GO:0001882,GO:0001883,GO:0003674,GO:0003824,GO:0003916,GO:0003918,GO:0005488,GO:0005524,GO:0005575,GO:0005622,GO:0005623,GO:0005737,GO:0006139,GO:0006152,GO:0006163,GO:0006195,GO:0006200,GO:0006259,GO:0006260,GO:0006261,GO:0006265,GO:0006351,GO:0006725,GO:0006753,GO:0006793,GO:0006796,GO:0006807,GO:0007059,GO:0008094,GO:0008150,GO:0008152,GO:0009056,GO:0009058,GO:0009059,GO:0009116,GO:0009117,GO:0009119,GO:0009123,GO:0009125,GO:0009126,GO:0009128,GO:0009141,GO:0009143,GO:0009144,GO:0009146,GO:0009150,GO:0009154,GO:0009158,GO:0009161,GO:0009164,GO:0009166,GO:0009167,GO:0009169,GO:0009199,GO:0009203,GO:0009205,GO:0009207,GO:0009259,GO:0009261,GO:0009295,GO:0009330,GO:0009987,GO:0010467,GO:0016070,GO:0016462,GO:0016787,GO:0016817,GO:0016818,GO:0016853,GO:0016887,GO:0017076,GO:0017111,GO:0018130,GO:0019438,GO:0019439,GO:0019637,GO:0019693,GO:0030554,GO:0032549,GO:0032550,GO:0032553,GO:0032555,GO:0032559,GO:0032774,GO:0032991,GO:0034641,GO:0034645,GO:0034654,GO:0034655,GO:0035639,GO:0036094,GO:0042221,GO:0042278,GO:0042454,GO:0042493,GO:0042623,GO:0043167,GO:0043168,GO:0043170,GO:0043234,GO:0044237,GO:0044238,GO:0044248,GO:0044249,GO:0044260,GO:0044270,GO:0044271,GO:0044281,GO:0044424,GO:0044464,GO:0044699,GO:0044710,GO:0044763,GO:0046034,GO:0046128,GO:0046130,GO:0046434,GO:0046483,GO:0046700,GO:0050896,GO:0055086,GO:0061505,GO:0071103,GO:0071704,GO:0072521,GO:0072523,GO:0090304,GO:0097159,GO:1901135,GO:1901136,GO:1901265,GO:1901292,GO:1901360,GO:1901361,GO:1901362,GO:1901363,GO:1901564,GO:1901565,GO:1901575,GO:1901576,GO:1901657,GO:1901658	K02470		bactNOG[38]	05C7D@bactNOG,0GB49@delNOG,0HC60@dproNOG,16Q2I@proNOG,COG0187@NOG	NA|NA|NA	L	DNA gyrase negatively supercoils closed circular double- stranded DNA in an ATP-dependent manner and also catalyzes the interconversion of other topological isomers of double-stranded DNA rings, including catenanes and knotted rings (By similarity)
	-------------

	# COG (SET OPTION -C in clu_funct_analysis_v2.sh: 
		cluster_function_enrichment_analysis_v4.py -C ... )
	# GO (SET OPTION -G IN clu_funct_analysis_v2.sh: 
		cluster_function_enrichment_analysis_v4.py -G ... )
	
		>> nohup bash clu_funct_analysis_v2.sh H > clu_funct_analysis_Hsap.out &		
		>> nohup bash clu_funct_analysis_v2.sh D > clu_funct_analysis_Dmel.out &
		>> nohup bash clu_funct_analysis_v2.sh A > clu_funct_analysis_Athal.out &
		>> nohup bash clu_funct_analysis_v2.sh S > clu_funct_analysis_Scer.out &

		E.g. H. sapiens, GO functions' analysis (clu_funct_analysis_v2.sh H):
			Running time: 9 x ~20 mins = ~180 mins
			Running-time for c=0.8: ~20 mins
			The resulting annotations file: 190 MB

		Example output - COG:
			PATH_TO_THE_ANALYSIS_FOLDER/cluster_hyper_cog/c_0_0 and A. thaliana:
			Athal_2022_cog.txt
			-------------
				pgi|0000000019284864591|ti|1928486|pi|0|
				K	1	3	1661	18327	2.33
				J	1	3	948	18327	3.20
				L	1	3	1275	18327	2.74
			-------------
			
			Athal_2022_cog_annotations.txt
			-------------			
				pgi|0000000001924736336|ti|1924736|pi|0|
				0000000010478227485#104782###T
				0000000001047827713#104782###T
				0000000010478242631#104782###T
				0000000003475156726#347515###I
			-------------
			
		Example output - GO:
			PATH_TO_THE_ANALYSIS_FOLDER/cluster_hyper_go/c_0_0 and A. thaliana:
			
			Athal_2022_go.txt
			-------------
				pgi|0000000001924736336|ti|1924736|pi|0|
				GO:0035637      1       215     315     974136  3.86
				GO:0016740      1       215     1234    974136  1.88
				GO:0016202      1       215     77      974136  5.90
				GO:0016773      1       215     350     974136  3.70
				GO:0051128      1       215     850     974136  2.42
			-------------

(4) Analysis of function enrichment per cluster
	IMPORTANT NOTE:
	# In script "cluster_function_enrichment_summary_cog_v2.py" change "PATH_TO_THE_ANALYSIS_FOLDER" to the correct path.
	# In script "cluster_function_enrichment_summary_go_v4.py" change "PATH_TO_THE_ANALYSIS_FOLDER" to the correct path.

	a) processing of function enrichment --per cluster-- for COG terms (computes enrichment for all focal species at once)
		>> time python3 cluster_function_enrichment_summary_cog_v2.py > cog_v2.out &
		
		Running time: ~3 mins
	
	b) processing of function enrichment --per cluster-- for GO terms (computes enrichment for a single f. species; -p --> plotting)
		>> time python3 cluster_function_enrichment_summary_go_v4.py -s H -c 0.8 -p > clu_funct_summary_HS_2022_c_0_8_go.out &
		>> time python3 cluster_function_enrichment_summary_go_v4.py -s A -c 0.8 -p > clu_funct_summary_Athal_2022_c_0_8_go.out &
		>> time python3 cluster_function_enrichment_summary_go_v4.py -s D -c 0.8 -p > clu_funct_summary_Dmel_2022_c_0_8_go.out &
		>> time python3 cluster_function_enrichment_summary_go_v4.py -s S -c 0.8 -p > clu_funct_summary_Scer_2022_c_0_8_go.out &
	
		E.g. H. sapiens and c=0.8 + plotting (cluster_function_enrichment_summary_go_v4.py -s H -c 0.8 -p)
			Running time: 6 mins
	
		c) calling cluster_function_enrichment_summary_go_v4.py for c=0.0 to c=0.8; no plotting (since plotting is too slow)
		>> nohup bash clu_funct_summary_v2.sh H > clu_funct_summary_Hsap.out &
		>> nohup bash clu_funct_summary_v2.sh S > clu_funct_summary_Scer.out &
		>> nohup bash clu_funct_summary_v2.sh D > clu_funct_summary_Dmel.out &
		>> nohup bash clu_funct_summary_v2.sh A > clu_funct_summary_Athal.out &


(5) Analysis of function enrichment per phylostratum
	# In script "phylostratum_function_enrichment_analysis_v2.py" change "PATH_TO_THE_ANALYSIS_FOLDER" to the correct path.
	# In script "phylostratum_function_enrichment_analysis_v3.py" change "PATH_TO_THE_ANALYSIS_FOLDER" to the correct path.

	a) enrichment analysis of COG functions gained and lost --per phylostratum-- with a hypergeometric test
		Set PARAMETER_C = 0.8, COG = True
		>> time python3 phylostratum_function_enrichment_analysis_v3.py > phylo_cog_0_8.out &  
		>> nohup time python3 phylostratum_function_enrichment_analysis_v2.py -s H -C -c 0.8 -p > clu_funct_analysis_HS_2022_c_0_8_cog.out &
	
		Running time for all species: ~3 mins
	
	b) enrichment analysis of GO functions gained and lost --per phylostratum-- with a hypergeometric test
		>> nohup bash phylostr_funct_summary_v2.sh A > phylostr_funct_summary_Athal.out &
			For each c-value =.0 to 0.8:
				--> calling phylostratum_function_enrichment_analysis_v2.py -s A -G -c 0.$c
		
		Do the same for all focal species (H, D, A, S)





