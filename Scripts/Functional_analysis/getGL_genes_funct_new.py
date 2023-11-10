__version__ = "1.0"
__author__ = "Mirjana Domazet-Lošo"
# February 6, 2023

"""
Usage:
	>> time python3 getGL_genes_funct_new.py -f 9606 \
					-p YOUR_PATH/Data/Parents \
					-i YOUR_PATH/Data/results_0_0/db_clu_all.tsv \
					-l YOUR_PATH/Data/allLCA.txt \
					-e YOUR_PATH/Data/eggnogResultsNew \
					-C \
					-o YOUR_PATH/Data/ClusterAnalysis/COG/0_0/res_dbAll_Hsap_0_0_clusters_funct.txt \
					-g YOUR_PATH/Data/ClusterAnalysis/COG/0_0/res_dbAll_Hsap_0_0_genes_funct.txt \
					-s YOUR_PATH/Data/ClusterAnalysis/COG/0_0/summary_GL_Hsap_0_0_funct.txt > out_genes_funct.txt
"""
# Count gene gain/lost for a focal species (given taxID, e.g. 9606 H. sap.) and annotate clusters and the genes in each cluster with COG/GO functions using emapper results
# # Rules for computing gene families' gain/loss  (Dollo's parsimony applied):
# 	if a cluster contains genes of the focal species
# 	(1a) if a cluster contains only members (genes) of the focal species, then I'll suppose that they originated at its minPS (== maxPS) (PS --> phylostratum), 
# 	(2) if a cluster also contains other species' genes, then the gene family was gained in minPS and still exists (not lost) 
#	
# 	(1b) if a cluster contains only members (genes) of one (non-focal) species, 
# 		then I'll suppose that they originated at this species PS (or its lineage) --> CAN BE IGNORED FROM THE POINT OF FOCAL SPECIES
#	
# 	(3) multiple species, and belonging to different PS; from the focal species' point of view: gene family gained at minPS, gene family lost at (maxPS+1)
#	
# 	(4) multiple species, and belonging to the same PS --> CAN BE IGNORED FROM THE POINT OF VIEW OF A FOCAL SPECIES
#	
# 	In conclusion:
# 	Only cases (1a), (2) and (3) are USED FOR COMPUTING GENE FAMILIES GAINED/LOST FROM THE FOCAL SPECIES PERSPECTIVE
	# 	Note: each cluster included in the gene families' gain/loss should have at least two members (singleton-clusters are not included in the analysis).

#import os
from os import walk
import getopt
import sys
from sys import getsizeof

DEBUG = False

def parseArgv(argv):
	print("Program arguments:")

	options = "hf:p:i:l:e:GCo:s:g:S" # options that require arguments should be followed by :
	long_options = ["help", "focal_species_id", "parents_directory", "input_tsv_file", "lca_file", "emapper_directory", "GO_function", "COG_function", "output_file", "summary_file", "genes_file", "strict_function_assignment"]
	try:
		opts, args = getopt.getopt(argv, options, long_options)
	except getopt.GetoptError:
		print("Error. Usage: getGL_genes_funct.py -f <focal_species_id> -p <parents_directory> -i <input_tsv_file> -l <lca_file> -e <emapper_directory> -G[|-C] -o <output_file> -s <summary_file> [-g <genes_file>] [-S]")
		sys.exit(2)
	
	FS_ID = ""
	PARENTS_DIR = ""
	INPUT_TSV_FILE = ""
	LCA_FILE = ""
	OUTPUT_FILE = ""
	SUMMARY_FILE = ""
	GENES_FILE = ""
	EMAPPER_DIR = ""
	STRICT = False
	GO = True
	
	for opt, arg in opts:
		print(opt + "\t" + arg)
		if opt in ("-h", "--help"):
			print("Usage: getGL_genes_funct.py -f <focal_species_id> -p <parents_directory> -i <input_tsv_file> -l <lca_file> -e <emapper_directory> -G[|-C] -o <output_file> -s <summary_file> [-g <genes_file>] [-S]")
			sys.exit()
		elif opt in ("-f", "--focal_species_id"):
			FS_ID = arg
		elif opt in ("-p", "--parents_directory"):
			PARENTS_DIR = arg
		elif opt in ("-i", "--input_tsv_file"):
			INPUT_TSV_FILE = arg
		elif opt in ("-l", "--lca_file"):
			LCA_FILE = arg
		elif opt in ("-e", "--emapper_directory"):
			EMAPPER_DIR = arg
		elif opt in ("-G", "--GO_function"):
			GO = True
		elif opt in ("-C", "--COG_function"):
			GO = False
		elif opt in ("-o", "--output_file"):
			OUTPUT_FILE = arg
		elif opt in ("-s", "--summary_file"):
			SUMMARY_FILE = arg
		elif opt in ("-g", "--genes_file"):
			GENES_FILE = arg
		elif opt in ("-S", "--strict_function_assignment"):
			STRICT = True
	# end for
	return FS_ID, PARENTS_DIR, INPUT_TSV_FILE, LCA_FILE, EMAPPER_DIR, GO, OUTPUT_FILE, SUMMARY_FILE, GENES_FILE, STRICT
# end parseArgv


# compute dictPS (input: parents of the focal species, e.g. Parents/9606.txt)
# for each focal species' phylostrata (1 to maxPS):
#	store its nodeId
def getPS(parentFile, dictPS):
	cnt = 0
	temp = [] # contains ps in reverse order e.g. from the focal species H. sap perspective: H. sap. will be in 0-th position, instead of 38-th (since it is at ps = 39)
	for line in parentFile: # eg. 9606	Homo_sapiens --> dictPS[1] = "9606"
		arr = line.rstrip().split()
		temp.append(arr[0])
		cnt += 1
	MAX_PS = cnt

	for i in range(MAX_PS):
		dictPS[temp[i]] = MAX_PS - i 

	return MAX_PS
# end getPS

# compute dictLCA (input: LCA file, focal species taxId)
#							880073	 526226	 526227	1924735 ...
#						Caldithrix Gordonia_b Meiothermu Salpingoec ...
#	880073 (Caldithrix)		 -1		  2		  2	 131567 ...
#	526226 (Gordonia_b)		  2		 -1		  2	 131567 ...
# for each species: store its LCA with the focal species
def getLCA(lcaFile, taxId, dictLCA):
	
	# read taxId from the first row
	line = next(lcaFile)
	listTaxId = line.strip().split()
	
	# skip the second row
	next(lcaFile)
	
	for line in lcaFile:
		arr = line.rstrip().split()
		if arr[0] == str(taxId):
			break
	
	for i in range(2, len(arr)):
		dictLCA[listTaxId[i - 2]] = arr[i]	
# end getLCA

# read genes' functions from emapper files and store them to a dictionary
def getGeneFunct(emapperDir, dictGeneFunct):
	
	# get all emapper files
	file_paths = []
	print("Reading files from folder: " + emapperDir)
	sys.stdout.flush()
	for (dirpath, dirnames, filenames) in walk(emapperDir):
		for f in filenames:
			if f.endswith(".annotations"):
				file_paths.append(dirpath + "/" + f)
	
	try:
		file_path = ""
		geneId = ""
		print(file_paths)
		file_paths_sorted = sorted(file_paths)
		print(file_paths_sorted)
		for file_path in file_paths_sorted:
			print("sizeof(dictGeneFunct) = " + str(getsizeof(dictGeneFunct)) + " bytes.")
			with open(file_path, "r") as f:
				if DEBUG: 
					print("Inputing categories from file " + file_path)
					sys.stdout.flush()
				#sys.stdout.flush()
				for line in f:
					if line.startswith("pid"):
						line_split = line.split("\t")
	
						#OLD version of emapper
						#if GO: # category string
						#	cat_str = line_split[5].strip()
						#else:
						#	cat_str = line_split[11].strip()
						
						# NEW version of emapper
						if GO: # category string
							cat_str = line_split[9].strip()
						else:
							cat_str = line_split[6].strip()
	
						if cat_str == "" or cat_str == "-" or (not GO and (cat_str == "R" or cat_str == "S")):
							continue
	
						if GO: # split list of categories
							cats = cat_str.split(",")
						else:
							cats = list(cat_str) # e.g. "EH" --> ['E', 'H']
						#id_elements = line_split[0].split("|")
						#geneId = id_elements[1] + "#" + id_elements[3]
						# or just take the whole geneId, i.e. line_split[0]
						geneId = line_split[0]
	
						if geneId not in dictGeneFunct:
							dictGeneFunct[geneId] = set()
	
						for cat in cats:
							cat = cat.strip()
							cat = cat.strip("\t")
							# if category/function unknown or COG category R or S
							# the followig two lines are obsolete, since the conditions are checked above (l.180)
							if len(cat) == 0: continue
							if not GO and (cat == "R" or cat == "S" or cat=="-"): continue
							# if GO function, remove "GO:"; append later during print
							if GO: cat = cat.strip("GO:")
							dictGeneFunct[geneId].add(cat)
				#end for
	except:
		print("Error occurred while reading emapper file " + file_path + "\t" + str(geneId))
		sys.stdout.flush()
			
# end getGeneFunct

# compute taxID from a geneId of the form "pid|0000000003120175935|tx|312017|"
def getTaxId(geneId):
	arr = geneId.split("|")
	return arr[3]
# end getTaxId

# (1) return the string containing a function, the number of genes to whom this function was assigned and the number of functions in a cluster that pass the threshold;
# (2) update dictFunctAllClusters
def getStrFunction(dictFunct, STRICT, numMembers, dictFunctAllClusters, GO):
	strFunct = ""
	if GO: prefixFunct = "GO:"
	else: prefixFunct = ""
	
	if STRICT: threshold = numMembers / 2
	else: threshold = 1
	if DEBUG: print("Threshold: " + str(threshold))
	
	if len(dictFunct) > 0: # if there are some functions assigned to a cluster
		thresholdNumFunctions = 0 # number of functions that are assigned to [at least 1 gene | at least half of the genes | some other criterion] genes in a cluster
		if STRICT:
			for funct in dictFunct: # for functions annotated for at least half of members
				if dictFunct[funct] >= threshold: # this criteria must be removed if all functions are included
					thresholdNumFunctions += 1 # number of functions that are assigned to at least half of the genes in a cluster
					
		else:
			thresholdNumFunctions = len(dictFunct) # all functions; no threshold
		# end if
		for funct in dictFunct: # for functions annotated for at least half of members
			if dictFunct[funct] >= threshold: 
				# print quant, sample_count, cat_2_total_count[cat_id], total_count
				strFunct += prefixFunct + funct + "\t" + str(dictFunct[funct]) + "\t" + str(thresholdNumFunctions) + "\n"
				if funct in dictFunctAllClusters:
					dictFunctAllClusters[funct] += 1
				else:
					dictFunctAllClusters[funct] = 1
				# write a separate file with summary data!!
	return strFunct
	
# end getStrFunction

# list of genes with assigned functions
def getStrGeneFunction(dictGeneFunct, GO, geneId):
	strFunctGene = ""
	#if geneId in dictGeneFunct: # if there are annotated functions for this gene
	for funct in dictGeneFunct[geneId]: # dictGeneFunct[geneId] is a list of functions
		if strFunctGene != "": strFunctGene += ","
		if GO: strFunctGene += "GO:"
		strFunctGene += funct
	# end for
	return strFunctGene 
# end getStrGeneFunction

# parse tsv file
# (1) print summary data:
#   for each pyhlostratum: number of gene families gained/lost)
# (2) print cluster data:
#   for each cluster: clusterId representative/geneId number-of-genes phylostratumGain pyhlostratumLoss
def parseTsv(tsvFile, distPS, dictLCA, outputFile, printClusters, summaryFile, genesFile, printGenes, MAX_PS, dictGeneFunct, STRICT, GO):

	prevReprId = ""
	clusterId = 0
	reprId = ""
	text = "" # just cluster representatives
	text2 = "" #  cluster representatives + list of their genes
	lineGenesIds = "" # list of genes in a cluster
	
	totalGFGained = [0 for i in range(MAX_PS)] # set total GF gained/lost counters to 0
	totalGFLost = [0 for i in range(MAX_PS)]
	dictFunctAllClusters = {} # count the number of clusters to which a function was assigned
	cntFunctClusters = 0 # the overall number of functions assigned to any of the clusters in the dataset

	# Count gene gain/lost for a focal species (given taxID, e.g. 9606 H. sap.) 
	# and annotate clusters and the genes in each cluster with COG/GO functions using emapper results
	# # Rules for computing gene families' gain/loss  (Dollo's parsimony applied):
	# 	if a cluster contains genes of the focal species
	# 	(1a) if a cluster contains only members (genes) of the focal species, then I'll suppose that they originated at its minPS (== maxPS) (PS --> phylostratum), 
	# 	(2) if a cluster also contains other species' genes, then the gene family was gained in minPS and still exists (not lost) 
	#	
	# 	(1b) if a cluster contains only members (genes) of one (non-focal) species, 
	# 		then I'll suppose that they originated at this species PS (or its lineage) --> CAN BE IGNORED FROM THE POINT OF FOCAL SPECIES
	#	
	# 	(3) multiple species, and belonging to different PS; from the focal species' point of view: gene family gained at minPS, gene family lost at (maxPS+1)
	#	
	# 	(4) multiple species, and belonging to the same PS --> CAN BE IGNORED FROM THE POINT OF VIEW OF A FOCAL SPECIES
	#	
	# 	In conclusion:
	# 	Only cases (1a), (2) and (3) are USED FOR COMPUTING GENE FAMILIES GAINED/LOST FROM THE FOCAL SPECIES PERSPECTIVE
	# 	Note: each cluster included in the gene families' gain/loss should have at least two members (singleton-clusters are not included in the analysis).
	
	try:
		for line in tsvFile: # representativeID taxID
			array = line.split()
			if array[0] != prevReprId: 
				prevReprId = reprId # old representative
			reprId = array[0]
			geneId = array[1]
			if reprId == geneId: # new cluster (the first memeber of the cluster is its representative)
				if DEBUG: print(line.rstrip() + "\n")
				#sys.stdout.flush()

				# if this was not the first cluster, process the previous cluster
				if clusterId > 0: 
					# printing clusterId representative/geneId number-of-genes phylostratumGain pyhlostratumLoss
					if maxPS == MAX_PS: psLost = "-"
					else: psLost = str(maxPS + 1)
					psGained = str(minPS)
					
					# print only those clusters that include the focal species and have >= 2 members or whose minPS != maxPS (GF was gained/lost along the way to the focal species)
					if (psLost == "-" and numMembers >= 2) or minPS != maxPS:
						if printClusters:
							if DEBUG: print(str(clusterId) + "\t" + prevReprId + "\t" + str(numMembers) + "\t" + psGained + "\t" + psLost + "\n")
							#print("STRICT: " + str(STRICT))
							if DEBUG: print(dictFunct)
							
							strFunct = getStrFunction(dictFunct, STRICT, numMembers, dictFunctAllClusters, GO)
							if DEBUG: print(dictFunctAllClusters)
							sys.stdout.flush()

							text += str(clusterId) + "\t" + prevReprId + "\t" + str(numMembers) + "\t" + psGained + "\t" + psLost + "\n" + strFunct + "\n"
						# end if printClusters
						if printGenes: # always print all functions for a gene (it doesn't matter if a function hasn't been assigned to a whole cluster)
							#print("lineGenesIds: " + lineGenesIds)
							text2 += str(clusterId) + "\t" + prevReprId + "\t" + str(numMembers) + "\t" + psGained + "\t" + psLost + "\n" + lineGenesIds + "\n"
						totalGFGained[minPS - 1] += 1 # -1 since indices start at 0
						if psLost != "-":
							totalGFLost[maxPS] += 1 # -1 since indices start at 0, but +1, since the cluster is lost maxPS + 1
					
					
					if clusterId % 10000 == 0: # print in batches
						if printClusters:
							outputFile.write(text)
						if printGenes: 
							genesFile.write(text2)
						text = ""
						text2 = ""
						print("Finished reading " + str(clusterId) + " clusters.")
						sys.stdout.flush()

				# set values for the processing of the new cluster
				minPS = MAX_PS
				maxPS = 0
				numMembers = 0
				clusterId += 1
				lineGenesIds = ""
				dictFunct = {} # count the number a function is assigned to a gene in the current cluster
				
			# current cluster
			numMembers += 1
			member_taxId = getTaxId(geneId)
			# update counts of annotated functions in the current cluster
			if geneId in dictGeneFunct: 
				for funct in dictGeneFunct[geneId]: # dictGeneFunct[geneId] is a list of functions
					if funct not in dictFunct:
						dictFunct[funct] = 1
					else: dictFunct[funct] += 1

			if printGenes:
				# add gene's functions
				if geneId in dictGeneFunct:
					strFunctGene = "\t" + getStrGeneFunction(dictGeneFunct, GO, geneId)
				else: strFunctGene = ""
				lineGenesIds += geneId + strFunctGene + "\n"
				
			# find phylostratum for member_taxId
			member_ps = dictLCA[member_taxId]
			if member_ps == "-1": ps = MAX_PS
			else: ps = int(dictPS[member_ps])
			if ps < minPS: minPS = ps
			if ps > maxPS: maxPS = ps
		# end for		
		tsvFile.close()
		
		# process the last cluster
		if maxPS == MAX_PS: psLost = "-"
		else: psLost = str(maxPS + 1)
		psGained = str(minPS)
		#if psLost == "-" or minPS != maxPS:
		# NEW; Nov 28, 2022 - print only those clusters that include the focal species and have at least 2 members or whose minPS != maxPS (GF was gained/lost along the way to the focal species)
		if (psLost == "-" and numMembers >= 2) or minPS != maxPS:
			if printClusters: 
				strFunct = getStrFunction(dictFunct, STRICT, numMembers, dictFunctAllClusters, GO)
				text += str(clusterId) + "\t" + prevReprId + "\t" + str(numMembers) + "\t" + psGained + "\t" + psLost + "\n" + strFunct
			if printGenes: 
				text2 += str(clusterId) + "\t" + prevReprId + "\t" + str(numMembers) + "\t" + psGained + "\t" + psLost + "\n" + lineGenesIds
			totalGFGained[minPS - 1] += 1 # -1 since indices start at 0
			if psLost != "-":
				totalGFLost[maxPS] += 1 # -1 since indices start at 0, but +1, since the cluster is lost maxPS + 1
		
		# write the last batch (if non-empty text)
		if text != "" and printClusters: 
			outputFile.write(text)
		if text2 != "" and printGenes: 
			genesFile.write(text2)		
		
		if printClusters: 
			outputFile.close()
		if printGenes:
			genesFile.close()
		
		# print totals - gene families gained/lost at each phylostratum from the perspective of the focal species
		summaryFile.write("ps" + "\t" + "GF_Gained" + "\t" + "GF_Lost" + "\n")
		for i in range(MAX_PS):
			summaryFile.write(str(i + 1) + "\t" + str(totalGFGained[i]) + "\t" + str(totalGFLost[i]) + "\n")
		summaryFile.write("\n\n")
		
		# for each category/function: print the total number of clusters to whom this function was assigned
		summaryFile.write("Category/Function" + "\t" + "Number of clusters" + "\n")
		lineFunct = ""
		cntLine = 0
		for funct in dictFunctAllClusters:
			if GO: lineFunct += "GO:"
			lineFunct += funct + "\t" + str(dictFunctAllClusters[funct]) + "\n"
			cntFunctClusters += dictFunctAllClusters[funct]
			cntLine += 1
			if cntLine % 10000 == 0: # print in batches
				summaryFile.write(lineFunct)
				lineFunct = ""
				cntLine = 0
		# print last batch
		summaryFile.write(lineFunct)
		# print total number of functions over all clusters
		summaryFile.write("\n\nTotal number of functions over all clusters:" + str(cntFunctClusters) + "\n")
		
		summaryFile.close()
	except:
		print("Error occurred while reading tsv file.")
		sys.stdout.flush()

# end parseTsv


# __main__
if __name__ == "__main__":
	print("__main__: Starting ... ")
	FS_ID, PARENTS_DIR, INPUT_TSV_FILE, LCA_FILE, EMAPPER_DIR, GO, OUTPUT_FILE, SUMMARY_FILE, GENES_FILE, STRICT = parseArgv(sys.argv[1:])
	#print("TSV input file:" + INPUT_TSV_FILE)
	#print("LCA input file:" + LCA_FILE)
	#print("Emapper input folder:" + EMAPPER_DIR)
	print("GO functions:" + str(GO))
	#print("Output file:" + OUTPUT_FILE)
	#print("Summary file:" + SUMMARY_FILE)
	#print("Genes file:" + GENES_FILE)
	sys.stdout.flush()
	
	# read phylostrata from the perspective of the focal species (FS_ID)
	parentFile = open(PARENTS_DIR + "/" + str(FS_ID) + "parents.txt", "r")
	dictPS = {}
	MAX_PS = getPS(parentFile, dictPS)
	print(dictPS)
	parentFile.close()
	
	# read LCA file and get dictLCA for the focal species
	lcaFile = open(LCA_FILE, "r")
	dictLCA = {}
	getLCA(lcaFile, FS_ID, dictLCA)
	#print(dictLCA)
	lcaFile.close()
	print("Finished reading " + LCA_FILE)
	sys.stdout.flush()
	
	# read genes' functions from emapper files and store them to a dictionary
	dictGeneFunct = {}
	getGeneFunct(EMAPPER_DIR, dictGeneFunct)
	print("Finished reading " + EMAPPER_DIR)
	#if DEBUG:
	#	cnt = 0
	#	for g in dictGeneFunct:
	#		print(g + "\t" + str(dictGeneFunct[g]))
	#		cnt += 1
	#		if cnt >= 100: break
	#	sys.stdout.flush()
	
	# parseTsv(tsvFile, listPS, outputFile, summaryFile, MAX_PS):
	tsvFile = open(INPUT_TSV_FILE, "r")
	summaryFile = open(SUMMARY_FILE, "w")
	
	outputFile = None
	printClusters = False
	if OUTPUT_FILE != "":
		outputFile = open(OUTPUT_FILE, "w")
		printClusters = True
	
	genesFile = None
	printGenes = False
	if GENES_FILE != "":
		genesFile = open(GENES_FILE, "w")
		printGenes = True
	#STRICT = True # strict criterion for assigning a category/function to a cluster

	print("Parsing tsv file.")
	sys.stdout.flush()
	parseTsv(tsvFile, dictPS, dictLCA, outputFile, printClusters, summaryFile, genesFile, printGenes, MAX_PS, dictGeneFunct, STRICT, GO)
	print("Done.")
	# tsvFile, outputFile, summaryFile are closed in parseTsv


