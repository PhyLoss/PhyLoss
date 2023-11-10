__version__ = "1.0"
__author__ = "Mirjana Domazet-Lošo"
# February 6, 2023

"""
Usage:
	>> time python3 parseClusterFunct.py -P 40 -S Hsap \
										-i PATH/ClusterAnalysis/COG/0_0/Hsap_0_0_clusters_funct.txt  \
										-c PATH/COGs.txt \
										-o PATH/ClusterAnalysis/COG/0_0
"""
# Parse results of the clusters' functional analysis and output the results in the form of tables and total counts to the output directory (-o)
# A separate table is produced for gain and the other one for loss; the total counts are also separately produced

#import os
from os import walk
import getopt
import sys
from sys import getsizeof

DEBUG = False

def parseArgv(argv):
	print("Program arguments:")

	options = "hP:S:i:c:o:" # options that require arguments should be followed by :
	long_options = ["help", "number-of-phylostrata", "input_species", "input_file", "categories_file", "output_directory"]
	try:
		opts, args = getopt.getopt(argv, options, long_options)
	except getopt.GetoptError:
		print("Error. Usage: parseClusterFunct.py -P <number-of-phylostrata> -S <input_species> -i <input_file> -c <categories_file> -o <output_directory>")
		sys.exit(2)
	
	INPUT_FILE = ""
	CATEGORIES_FILE = ""
	OUTPUT_DIR = ""
	MAX_PS = 0
	INPUT_SPECIES = ""
	
	for opt, arg in opts:
		print(opt + "\t" + arg)
		if opt in ("-h", "--help"):
			print("Usage: parseClusterFunct.py -P <number-of-phylostrata> -S <input_species> -i <input_file> -c <categories_file> -o <output_directory>")
			sys.exit()
		elif opt in ("-i", "--input_file"):
			INPUT_FILE = arg
		elif opt in ("-o", "--output_directory"):
			OUTPUT_DIR = arg
		elif opt in ("-c", "--categories_file"):
			CATEGORIES_FILE = arg
		elif opt in ("-P", "--number-of-phylostrata"):
			MAX_PS = arg
		elif opt in ("-S", "--input_species"):
			INPUT_SPECIES = arg
	# end for
	return INPUT_FILE, OUTPUT_DIR, CATEGORIES_FILE, MAX_PS, INPUT_SPECIES
# end parseArgv

# get COGs (id + name)
def getCOGs(categoriesFile):
	dictCatNames = {}
	if categoriesFile:
		with open(categoriesFile, "r") as f:
			for line in f:
				# e. g. [D] Cell cycle control, cell division, chromosome partitioning
				if "[R]" in line or "[S]" in line:
					continue # skip: [R] General function prediction only; [S] Function unknown
				elif "[" in line:
					id = line[1]
					dictCatNames[id] = line[4:].rstrip() # from the 4th position to the eol
	print("Finished reading categories.", flush = True)

	# explicitly add id = "no_annotation", name = "Cluster without assigned function"
	dictCatNames["no_annotation"] = "Cluster without assigned function"	
	return dictCatNames
# end getCOGs

# get GOs (id + name)
def getGOs(categoriesFile):
	dictCatNames = {}
	
	# input e.g.
	# [Term]
	# id: GO:0000001
	# name: mitochondrion inheritance
	# namespace: biological_process
	# ...

	if categoriesFile:
		with open(categoriesFile, "r") as f:
			id = ""
			for line in f:
				if line.startswith("id: "):
					# get id without "GO:", e.g. "id: GO:0000001" --> "0000001"
					id = line.strip("id: ").strip() # strip(characters) --> remove a set of characters at leading/trailing positions
					if "GO:" in id:
						id = id.strip("GO: ").strip()
					else:
						id = "ignore"
				if line.startswith("name: ") and id != "ignore":
					name = line[6:].strip()
					dictCatNames[id] = name
	dictCatNames["no_annotation"] = "Cluster without assigned function"
	#print(dictCatNames)
	return dictCatNames
# end getGOs

# auxiliary function called from parseClusters
def updDictGainLoss(dictCatPS, dictCatTotal, c, ps):
	if c not in dictCatPS:
		dictCatPS[c] = {}
	
	if ps not in dictCatPS[c]:
		dictCatPS[c][ps] = 0
	dictCatPS[c][ps] += 1

	#if not (gain_loss == "loss" and ps == "-"):
	#	if c not in dictCatTotal: dictCatTotal[c] = 0
	#	dictCatTotal[c] += 1
	if c not in dictCatTotal: 
		dictCatTotal[c] = 0
	dictCatTotal[c] += 1

# end updDictGainLoss

# parse input file, i.e. clusters
def parseClusters(inputFile, dictCatPSGain, dictCatPSLoss, dictCatTotalGain, dictCatTotalLoss, cog_go):

	# Input. e.g.
	#2187	pid|0000000008153320723|tx|81533|	31	10	15
	#
	#2281	pid|0000000000945713148|tx|945713|	2	1	7
	#T	1	3
	#O	1	3
	#G	1	3
	prev = ""
	try:
		for line in inputFile: # representativeID taxID
			if line in ['\n', '\r\n']:
				if prev == "Cluster line": # if the prev. cluster was without annotations (i.e. this cuurent line is below the cluster data)
					updDictGainLoss(dictCatPSGain, dictCatTotalGain, "no_annotation", psGain)
					if psLoss != "-": # count only those categories/functions that were lost
						updDictGainLoss(dictCatPSLoss, dictCatTotalLoss, "no_annotation", psLoss)
				if prev !="" and DEBUG: 
					print(dictCatPSGain)
					print(dictCatPSLoss)
					print(dictCatTotalGain)
					print(dictCatTotalLoss)
				prev = ""
			else: # non-empty line
				if DEBUG: print(line.rstrip())

				arr = line.split()
				if "pid" in line: # cluster data: id reprId numMembers psGain psLoss
					psGain = arr[3]
					psLoss = arr[4] # possibly "-"
					prev = "Cluster line"
				else: # category/function line
					if cog_go == "cog": c = arr[0] 
					else: c = arr[0][3:].strip() # input without "GO:"

					updDictGainLoss(dictCatPSGain, dictCatTotalGain, c, psGain)
					if psLoss != "-":
						updDictGainLoss(dictCatPSLoss, dictCatTotalLoss, c, psLoss)
					prev = "Function line"

		if prev !="" and DEBUG: 
			print(dictCatPSGain)
			print(dictCatPSLoss)
			print(dictCatTotalGain)
			print(dictCatTotalLoss)
	except:
		print("Error occurred while reading input file.", flush = True)
# end parseClusters


def printCategoriesPS(outputFile, dictCatNames, dictCatPS, MAX_PS, cog_go):
	with open(outputFile, "w") as out:
		ps_list = ""
		for ps in range (1, int(MAX_PS)):
			ps_list += str(ps) + "\t"
		ps_list += str(MAX_PS)
		
		if DEBUG: print("dictCatPS\n" + str(dictCatPS))
		if cog_go == "go": prefix_c = "GO:"
		else: prefix_c = ""

		out.write("category_id\tname\t" + ps_list + "\n")
		for c in dictCatNames:
			cnt_list = ""
			#if c not in dictCatPS:
			#	for ps in range (1, int(MAX_PS)):
			#		cnt_list += "0\t"
			#	cnt_list += "0"
			#	
			#else:
			if c in dictCatPS: # print only functions that were assigned
				for ps in range (1, int(MAX_PS)):
					if str(ps) in dictCatPS[c]:
						cnt_list += str(dictCatPS[c][str(ps)]) + "\t"
					else:
						cnt_list += "0\t"
				# end for
				if MAX_PS in dictCatPS[c]:
					cnt_list += str(dictCatPS[c][MAX_PS])
				else:
					cnt_list += "0"

				if c == "no_annotation": prefix_c = "" # delete GO prefix in this case
				out.write(prefix_c + c + "\t" + dictCatNames[c] + "\t" + cnt_list + "\n")
			
			#out.write(c + "\t" + dictCatNames[c] + "\t" + cnt_list + "\n")
# end printCategoriesPS

def printCategoriesTotal(outputFile, dictCatNames, dictCatTotal, cog_go):
	with open(outputFile, "w") as out:
		out.write("category\ttotal_count\n")

		if cog_go == "go": prefix_c = "GO:"
		else: prefix_c = ""		

		for c in dictCatNames:
			if c in dictCatTotal:
				if c == "no_annotation": prefix_c = "" # delete GO prefix in this case
				out.write(prefix_c + c + "\t" + str(dictCatTotal[c]) + "\n")
			#else:
			#	out.write(prefix_c + c + "\t0\n")
			
# end printCategoriesTotal


# __main__
if __name__ == "__main__":
	print("__main__: Starting ... ")
	INPUT_FILE, OUTPUT_DIR, CATEGORIES_FILE, MAX_PS, INPUT_SPECIES = parseArgv(sys.argv[1:])
	dictCatNames = {}
	#dictCatNames = getCOGs(CATEGORIES_FILE)
	#print(dictCatNames, flush = True)
	
	if "COG" in CATEGORIES_FILE: 
		COG_GO = "cog" 
		dictCatNames = getCOGs(CATEGORIES_FILE) # COGs.txt
	else:
		COG_GO = "go" 	
		dictCatNames = getGOs(CATEGORIES_FILE) # go-basic.obo
	#print(dictCatNames, flush = True)
	
	# parse clusters to get number of functions gained|lost for each phylostratum
	inputFile = open(INPUT_FILE)
	dictCatPSGain = {}
	dictCatPSLoss = {}
	dictCatTotalGain = {}
	dictCatTotalLoss = {}
	parseClusters(inputFile, dictCatPSGain, dictCatPSLoss, dictCatTotalGain, dictCatTotalLoss, COG_GO)
	inputFile.close()
	
	# print results - ps and categories
	outputFile = OUTPUT_DIR + "/" + "phyl_multi_fun_table_" + INPUT_SPECIES + "_gain_" + COG_GO + ".txt"
	printCategoriesPS(outputFile, dictCatNames, dictCatPSGain, MAX_PS, COG_GO)

	# categories lost at ps == "-" won't be printed
	outputFile = OUTPUT_DIR + "/" + "phyl_multi_fun_table_" + INPUT_SPECIES + "_loss_" + COG_GO + ".txt"
	printCategoriesPS(outputFile, dictCatNames, dictCatPSLoss, MAX_PS, COG_GO)

	# print results - total categories
	# categories lost at ps == "-" won't be printed
	outputFile = OUTPUT_DIR + "/" + "phyl_multi_fun_total_counts_" + INPUT_SPECIES + "_gain_" + COG_GO + ".txt"
	printCategoriesTotal(outputFile, dictCatNames, dictCatTotalGain, COG_GO)
	
	outputFile = OUTPUT_DIR + "/" + "phyl_multi_fun_total_counts_" + INPUT_SPECIES + "_loss_" + COG_GO + ".txt"
	printCategoriesTotal(outputFile, dictCatNames, dictCatTotalLoss, COG_GO)
	print("Done.")














