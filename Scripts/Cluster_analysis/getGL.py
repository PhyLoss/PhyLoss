__version__ = "1.0"
__author__ = "Mirjana Domazet-Lošo"
# October 27, 2022
# Previous versions: 2017-2019

"""
Usage: time python3 getGL.py -f 9606 -p YOUR_PATH/Data/Parents/ -i YOUR_PATH/Data/tsv_new/results_0_0/db_clu_all.tsv \
					-l YOUR_PATH/Data/allLCA.txt -o YOUR_PATH/Data/res_dbAll_H_sap_0_0.txt -s YOUR_PATH/Data/summary_GL_H_sap_0_0.txt > out.txt
"""
# Count gene gain/lost for a focal species (given taxID, e.g. 9606 H. sap. and the number of its phylostrata: 39)
# Rules:
	# if (foundT) { # a cluster contains genes of the focal species
		# (1a) if a cluster contains only members (genes) of the focal species, then I'll suppose that they originated at its minPS (== maxPS), and are not lost (clusterRepresentative is from the focal species)
		# (2) if a cluster also contains other species' genes, then the gene family was gained in minPS and still exists (not lost) (clusterRepresentative doesn't have to be from the focal species)
	
	# (1b) if a cluster contains only members (genes) of one (non-focal) species, 
	# then I'll suppose that they originated at the species PS (or its lineage) --> CAN BE IGNORED FROM THE POINT OF FOCAL SPECIES
	
	# (3) multiple species, and belonging to different PS; from t's point of view: gene family gained at minPS, gene family lost at (maxPS+1)

	# (4) multiple species, and belonging to the same PS --> CAN BE IGNORED FROM THE POINT OF FOCAL SPECIES
		# (4a) if LCA of these species equals the cluster's minPS (== maxPS), than from t's point of view: gene family gained at minPS, gene family lost at (maxPS+1)
		# (4b) if LCA of these species does not equal the cluster's minPS (== maxPS), than UNRESOLVED (i.e. for example their LCA with the focal species is "Cell. organisms", but their LCA is "Viridiplantae")
	
	# In conclusion:
	# Only cases (foundT) (i.e. 1a and 2) and (3) are USED FOR COMPUTING GENE FAMILIES GAINED/LOST FROM THE FOCAL SPECIES PERSPECTIVE

import os
import getopt
import sys

def parseArgv(argv):
	print("Program arguments:")

	options = "hf:p:i:l:o:s:" # options that require arguments should be followed by :
	long_options = ["help", "focal_species_id", "parents_directory", "input_tsv_file", "lca_file", "output_file", "summary_file"]
	try:
		opts, args = getopt.getopt(argv, options, long_options)
	except getopt.GetoptError:
		print("Error. Usage: getGL.py -f <focal_species_id> -p <parents_directory> -i <input_tsv_file> -l <lca_file> -o <output_file> -s <summary_file>")
		sys.exit(2)
	
	FS_ID = ""
	PARENTS_DIR = ""
	INPUT_TSV_FILE = ""
	LCA_FILE = ""
	OUTPUT_FILE = ""
	SUMMARY_FILE = ""
	
	for opt, arg in opts:
		print(opt + "\t" + arg)
		if opt in ("-h", "--help"):
			print("Usage: getGL.py -f <focal_species_id> -p <parents_directory>  -i <input_tsv_file> -l <lca_file> -o <output_file> -s <summary_file>")
			sys.exit()
		elif opt in ("-f", "--focal_species_id"):
			FS_ID = arg
		elif opt in ("-p", "--parents_directory"):
			PARENTS_DIR = arg
		elif opt in ("-i", "--input_tsv_file"):
			INPUT_TSV_FILE = arg
		elif opt in ("-l", "--lca_file"):
			LCA_FILE = arg
		elif opt in ("-o", "--output_file"):
			OUTPUT_FILE = arg
		elif opt in ("-s", "--summary_file"):
			SUMMARY_FILE = arg
	# end for
	return FS_ID, PARENTS_DIR, INPUT_TSV_FILE, LCA_FILE, OUTPUT_FILE, SUMMARY_FILE
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

	return MAX_PS;
# edn getPS

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

# compute taxID from a geneId of the form "pid|0000000003120175935|tx|312017|"
def getTaxId(geneId):
	arr = geneId.split("|")
	return arr[3]
# end getTaxId


# parse tsv file
# (1) print summary data:
#   for each pyhlostratum: number of gene families gained/lost)
# (2) print cluster data:
#   for each cluster: clusterId representative/geneId number-of-genes phylostratumGain pyhlostratumLoss
def parseTsv(tsvFile, distPS, dictLCA, outputFile, summaryFile, MAX_PS):

	prevReprId = ""
	clusterId = 0
	reprId = ""
	text = ""
	totalGFGained = [0 for i in range(MAX_PS)] # set total GF gained/lost counters to 0
	totalGFLost = [0 for i in range(MAX_PS)]

	for line in tsvFile: # representativeID taxID
		array = line.split()
		if array[0] != prevReprId: 
			prevReprId = reprId # old representative
		reprId = array[0]
		geneId = array[1]
		if reprId == geneId: # new cluster (the first memeber of the cluster is its representative)
			
			# if this was not the first cluster, process the previous cluster
			if clusterId > 0: 
				# printing clusterId representative/geneId number-of-genes phylostratumGain pyhlostratumLoss
				if maxPS == MAX_PS: psLost = "-"
				else: psLost = str(maxPS + 1)
				psGained = str(minPS)
				# print only those clusters that include the focal species or whose minPS != maxPS (GF was gained/lost along the way to the focal species)
				# these clusters are also those included in the total counts of GF gained/lost at a specific phylostratum
				if psLost == "-" or minPS != maxPS:
					text += str(clusterId) + "\t" + prevReprId + "\t" + str(numMembers) + "\t" + psGained + "\t" + psLost + "\n"
					totalGFGained[minPS - 1] += 1 # -1 since indices start at 0
					if psLost != "-":
						totalGFLost[maxPS] += 1 # -1 since indices start at 0, but +1, since the cluster is lost maxPS + 1
				if clusterId % 10000 == 0: # print in batches
					outputFile.write(text)
					text = ""

			# set values for the processing of the new cluster
			minPS = MAX_PS
			maxPS = 0
			numMembers = 0
			clusterId += 1
			
		# current cluster
		numMembers += 1
		member_taxId = getTaxId(geneId)
		# find phylostratum for member_taxId
		member_ps = dictLCA[member_taxId]
		if member_ps == "-1": ps = MAX_PS
		else: ps = int(dictPS[member_ps])
		if ps < minPS: minPS = ps
		if ps > maxPS: maxPS = ps
	# end for		
	
	# process the last cluster
	if maxPS == MAX_PS: psLost = "-"
	else: psLost = str(maxPS + 1)
	psGained = str(minPS)
	if psLost == "-" or minPS != maxPS:
		text += str(clusterId) + "\t" + reprId + "\t" + str(numMembers) + "\t" + psGained + "\t" + psLost + "\n"
		totalGFGained[minPS - 1] += 1 # -1 since indices start at 0
		if psLost != "-":
			totalGFLost[maxPS] += 1 # -1 since indices start at 0, but +1, since the cluster is lost maxPS + 1
	
	# write the last batch (if non-empty text)
	if text != "": 
		outputFile.write(text) 
	tsvFile.close()
	outputFile.close()
	
	# print totals - gene families gained/lost at each phylostratum from the perspective of the focal species
	summaryFile.write("ps" + "\t" + "GF_Gained" + "\t" + "GF_Lost" + "\n")
	for i in range(MAX_PS):
		summaryFile.write(str(i + 1) + "\t" + str(totalGFGained[i]) + "\t" + str(totalGFLost[i]) + "\n")
	
	summaryFile.close()
# end parseTsv


# __main__
if __name__ == "__main__":
	print("__main__: Starting ... ")
	FS_ID, PARENTS_DIR, INPUT_TSV_FILE, LCA_FILE, OUTPUT_FILE, SUMMARY_FILE = parseArgv(sys.argv[1:])
	
	# read phylostrata from the perspective of the focal species (FS_ID)
	parentFile = open(PARENTS_DIR + str(FS_ID) + "parents.txt", "r")
	dictPS = {}
	MAX_PS = getPS(parentFile, dictPS)
	#print(dictPS)
	parentFile.close()
	
	# read LCA file and get dictLCA for the focal species
	lcaFile = open(LCA_FILE, "r")
	dictLCA = {}
	getLCA(lcaFile, FS_ID, dictLCA)
	#print(dictLCA)
	lcaFile.close()

	# parseTsv(tsvFile, listPS, outputFile, summaryFile, MAX_PS):
	tsvFile = open(INPUT_TSV_FILE, "r")
	outputFile = open(OUTPUT_FILE, "w")
	summaryFile = open(SUMMARY_FILE, "w")
	parseTsv(tsvFile, dictPS, dictLCA, outputFile, summaryFile, MAX_PS)
	# tsvFile, outputFile, summaryFile are closed in parseTsv
