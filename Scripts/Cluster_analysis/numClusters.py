__version__ = "1.0"
__author__ = "Mirjana Domazet-Lošo"
# April 15, 2022

"""
Usage: python3 numClusters.py -i names_All_2022.txt -j . -t PATH_TO_ANALYSIS_FOLDER/results_0_8/db_clu_all.tsv -o PATH_TO_ANALYSIS_FOLDER/results_0_8/
--> -i option: input file with all species' taxId-s + names
--> -j option: input directory where ps files for each focal species are stored (ps files: taxId_ps_Athal.txt, taxId_ps_Dmel.txt, ...)
--> -c option: input file with the results of mmseqs clustering (.tsv file)
--> -d option: directory for output files of the form: dbAllPlus_clu_0_8_out.txt (one file for each focal species)
"""

import os
import getopt
import sys

# ps files
PS_FILES = {"H" : "taxId_ps_Hsap.txt", "A" : "taxId_ps_Athal.txt", "D" : "taxId_ps_Dmel.txt", "S" : "taxId_ps_Scer.txt"}
FS = {"H" : "Hsap", "A" : "Athal", "D" : "Dmel", "S" : "Scer"}
FS_TAXID = {"H" : 9606, "A" : 3702, "D" : 7227, "S" : 559292}

#PS_FILES = {"H" : "taxId_ps_Hsap.txt"}
#FS = {"H" : "Hsap"}
#FS_TAXID = {"H" : 9606}

def parseArgv(argv):
	print("Program arguments:")

	options = "hi:i:j:t:o:" # options that require arguments should be followed by :
	long_options = ["help", "names-file", "ps-directory", "tsv-file", "output-directory"]
	try:
		opts, args = getopt.getopt(argv, options, long_options)
	except getopt.GetoptError:
		print("Error. Usage: numClusters.py -i <names-file> -j <ps-directory> -t <tsv-file> -o <output-directory>")
		sys.exit(2)
	
	NAMES_FILE = ""
	PS_DIR = ""
	TSV_FILE = ""
	OUTPUT_DIR = ""
	
	for opt, arg in opts:
		print(opt + "\t" + arg)
		if opt in ("-h", "--help"):
			print("Usage: numClusters.py -i <names-file> -j <ps-directory> -t <tsv-file> -o <output-directory>")
			sys.exit()
		elif opt in ("-i", "--names-file"):
			NAMES_FILE = arg
		elif opt in ("-j", "--ps-directory"):
			PS_DIR = arg
		elif opt in ("-t", "--tsv-file"):
			TSV_FILE = arg
		elif opt in ("-o", "--output-directory"):
			OUTPUT_DIR = arg
	# end for
	return NAMES_FILE, PS_DIR, TSV_FILE, OUTPUT_DIR
# end parseArgv


# read ps and names file
def getInput(inputFile, dictPS):
	for line in inputFile: # eg. 10020	6
		array = line.split()
		dictPS[array[0]] = array[1]
# end getInput


# copute taxId from geneId
def getTaxId(geneId):
	#if geneId == pgi|0000000003120175935|ti|312017|pi|0| then get taxid 312017
	array = geneId.split("|")
	return array[3]	
# end getTaxId

def getTsv(tsvFile, outputFiles, dictNames, dictAll):

	prevReprId = ""
	members = []
	clusterId = 0
	for line in tsvFile: # representativeID taxID
		array = line.split()
		reprId = array[0]
		geneId = array[1]
		
		if reprId == geneId: # new cluster (the first memeber of the cluster is its representative)
			if len(members) > 0: # this was not the first cluster, so the data of the previous cluster should be processed
				# printing
				taxId = getTaxId(members[0])
				for fs in FS:
					outputFile = outputFiles[fs]
					# choose only those clusters that contain at least 2 members or contain only focal species (in this case, the representative must be also member of the focal sp.)
					if len(members) > 1 or int(taxId) == int(FS_TAXID[fs]):
						outputFile.write("(" + str(clusterId) + ") Cluster representative: " + members[0] + "\tNumber of members in the cluster: " + str(len(members)) + "\n")					
						# print members
						outputFile.write("\tCluster members: \n")
						for j in range(len(members)):
							memberTaxid = getTaxId(members[j])
							outputFile.write("\t(" + str(j+1) + ") " + members[j] + "\n\t" + dictNames[memberTaxid] + "\t" + dictAll[fs][memberTaxid] + "\n")
						# end for
						outputFile.write("\n\n")				
			# end if
			members.clear() # delete all members
			clusterId += 1
			members.append(geneId) # add first member (it is representative itself)
			
		else: # current cluster
			members.append(geneId)

# end getTsv


# __main__
if __name__ == "__main__":
	print("__main__: Starting ... ")
	NAMES_FILE, PS_DIR, TSV_FILE, OUTPUT_DIR = parseArgv(sys.argv[1:])
	
	dictAll = {}	
	for key in PS_FILES:
		inputFileName = PS_DIR + "/" + PS_FILES[key]
		inputFile = open(inputFileName, "rt")
		dictAll[key] = {}
		getInput(inputFile, dictAll[key])
		print("Read " + inputFileName + ".")
		inputFile.close()
		#break
	
	dictNames = {}
	namesFile = open(NAMES_FILE, "rt")
	getInput(namesFile, dictNames) # 526226	Gordonia_bronchialis_dsm_43247
	namesFile.close()
		
	# process tsv file and construct output for each focal species
	tsvFile = open(TSV_FILE, "rt")
	outputFiles = {}
	for key in FS:
		i = OUTPUT_DIR.find("0_")
		c = OUTPUT_DIR[i:(i+3)] # e.g "0_8"
		outputFiles[key] = open(OUTPUT_DIR + FS[key] + "/" + "dbAllPlus_clu_" + c + "_" + FS[key] + "_out.txt", "wt") # dbAllPlus_clu_0_8_HS_out.txt
	getTsv(tsvFile, outputFiles, dictNames, dictAll)
	tsvFile.close()
	
	for key in FS:
		outputFiles[key].close()
	print("Done.")
