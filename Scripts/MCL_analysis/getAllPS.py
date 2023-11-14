__version__ = "1.0"
__author__ = "Mirjana Domazet-Lošo"
# 2022-2023

"""
Usage: python3 getAllPS.py -f 9606 -i allParents.txt -p ./Parents -o taxId_ps_Hsap.txt
"""


import os
import getopt
import sys

CELL_ORGANISMS = 131567

def parseArgv(argv):
    print("Program arguments:")

    options = "hf:i:p:o:" # options that require arguments should be followed by :
    long_options = ["help", "focal-species", "input-file", "parents-directory", "output-file"]
    try:
        opts, args = getopt.getopt(argv, options, long_options)
    except getopt.GetoptError:
        print("Error. Usage: getAllPS.py -f <focal-species> -i <input-file> -p <parents-directory> -o <output-file>")
        sys.exit(2)
    
    FOCAL_SPECIES = ""
    INPUT_FILE = ""
    PARENTS_DIR = ""
    OUTPUT_FILE = ""
    
    for opt, arg in opts:
        print(opt + "\t" + arg)
        if opt in ("-h", "--help"):
            print("Usage: getAllPS.py -f <focal-species> -i <input-file> -p <parents-directory> -o <output-file>")
            sys.exit()
        elif opt in ("-f", "--focal-species"):
            FOCAL_SPECIES = arg
        elif opt in ("-i", "--input-file"):
            INPUT_FILE = arg
        elif opt in ("-p", "--parents-directory"):
            PARENTS_DIR = arg
        elif opt in ("-o", "--output-file"):
            OUTPUT_FILE = arg
    # end for
    return FOCAL_SPECIES, INPUT_FILE, PARENTS_DIR, OUTPUT_FILE
# end parseArgv


def printOutput(inputFile, outputFile, dictPS, FOCAL_SPECIES):
	dict1 = {}    
	
	nextTaxon = True    	
	for line in inputFile: # e.g. 10020	Dipodomys_ordii
		array = line.split()
		nodeId = array[0]	 
		if nextTaxon:
			taxId = nodeId
			nextTaxon = False
			#print("taxId:" + taxId)
			if int(taxId) == int(FOCAL_SPECIES):
				dict1[taxId] = dictPS[taxId]
		else:
			if taxId not in dict1: 
				if nodeId in dictPS: # get the lowest common ancestor (maximal possible ps) of the focal species and the current taxon (including the focal species itself)
					dict1[taxId] = dictPS[nodeId]
			if int(nodeId) == CELL_ORGANISMS:
				nextTaxon = True # get ready for the next taxon
				#print("Next")
	# end for		
	for key in dict1:
		outputFile.write(key + "\t" + str(dict1[key]) + "\n")			
	# end for	
# end printOutput

# read focal species ps
def getFocalSpeciesPS(psFile, dictPS):
	cnt = 0
	listPS = []
	for line in psFile: # e.g. 9606	Homo_sapiens ... up to ... 131567	Cellular_organisms
		array = line.split()
		listPS.append(array[0]) # add 9606
		cnt += 1
	#print("cnt = " + str(cnt))

	for i in range(0, cnt):
		ps = cnt - i
		dictPS[listPS[i]] = ps # e.g. dictPS[9606] = 39
	#print(dictPS)
# end getFocalSpeciesPS


# __main__
if __name__ == "__main__":
	print("__main__: Starting ... ")
	FOCAL_SPECIES, INPUT_FILE, PARENTS_DIR, OUTPUT_FILE = parseArgv(sys.argv[1:])
	psFile = open(PARENTS_DIR + "/" + FOCAL_SPECIES + "parents.txt", "rt")
	# read focal species ps
	dictPS = {}
	getFocalSpeciesPS(psFile, dictPS)
	psFile.close()
	
	inputFile = open(INPUT_FILE, "rt")
	outputFile = open(OUTPUT_FILE, "wt")
	printOutput(inputFile, outputFile, dictPS, FOCAL_SPECIES)
	print("Done.")
	inputFile.close()
	outputFile.close()


