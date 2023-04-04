__version__ = "1.0"
__author__ = "Mirjana Domazet-Loso"

#
# get_NodesWith1Child.py - find parental nodes with only one child/taxon
# November, 2022
#							
#	--> Run: python3 upd_NodesWith1Child.py -i /storage/home/mdomazet/ProteinsT/Data/nodes_All_2022.txt \
#													-j PATH/Data/names_All_2022.txt \
#													-t PATH/Data/taxIDNames_All_2022.txt \
#													-o PATH/Data/nodes_with1Child.txt \
#													-n PATH/Data/nodes_All_2022_UPD.txt \ 
#													-m PATH/Data/names_All_2022_UPD.txt \
#													-s PATH/Data/taxIDNames_All_2022_UPD.txt > out.txt
#
import os, sys, getopt

def parseArgv(argv):
    print("Program arguments:")

    options = "hi:j:t:o:n:m:s:" # options that require arguments should be followed by :
    long_options = ["help", "in_nodes_file", "in_names_file", "taxa_file", "out_file", "nodes_upd_file", "names_upd_file", "taxa_upd_file"]
    try:
        opts, args = getopt.getopt(argv, options, long_options)
    except getopt.GetoptError:
        print("Error. Usage: upd_NodesWith1Child.py -i <in_nodes_file> -j <in_names_file> -t <taxa_file> -o <out_file> -n <nodes_upd_file> -m <names_upd_file> -s <taxa_upd_file>")
        sys.exit(2)
    
    NODES_FILE = ""
    TAXA_FILE = ""
    NAMES_FILE = "" 
    OUTPUT_FILE = "" 
    NODES_UPD_FILE = "" 
    NAMES_UPD_FILE = "" 
    TAXA_UPD_FILE = "" 

    for opt, arg in opts:
        print(opt + "\t" + arg)
        if opt in ("-h", "--help"):
            print("Usage: upd_NodesWith1Child.py -i <in_nodes_file> -j <in_names_file> -t <taxa_file> -o <out_file> -n <nodes_upd_file> -m <names_upd_file> -s <taxa_upd_file>")
            sys.exit()
        elif opt in ("-i", "--in_nodes_file"):
            NODES_FILE = arg
        elif opt in ("-j", "--in_names_file"):
            NAMES_FILE = arg
        elif opt in ("-t", "--taxa_file"):
            TAXA_FILE = arg
        elif opt in ("-o", "--out_file"):
            OUTPUT_FILE = arg
        elif opt in ("-n", "--nodes_upd_file"):
            NODES_UPD_FILE = arg
        elif opt in ("-m", "--names_upd_file"):
            NAMES_UPD_FILE = arg
        elif opt in ("-s", "--taxa_upd_file"):
            TAXA_UPD_FILE = arg
    # end for
    return NODES_FILE, NAMES_FILE, TAXA_FILE, OUTPUT_FILE, NODES_UPD_FILE, NAMES_UPD_FILE, TAXA_UPD_FILE
# end parseArgv

# read data from file to dict
def readFile(inputFile, dict1):
	inFile = open(inputFile, "r")
	for line in inFile:
		arr = line.strip().split("\t")
		dict1[arr[0]] = arr[1]
	inFile.close()
# end readFile

# write data from dict to file
def writeFile(outputFile, dict1):
	outFile = open(outputFile, "w")
	text = ""
	for id in dict1:
		text += id + "\t" + dict1[id] + "\n"
	outFile.write(text)
	outFile.close()
	#print(dict1)
# end readFile


# main
if __name__ == "__main__":
	NODES_FILE, NAMES_FILE, TAXA_FILE, OUTPUT_FILE, NODES_UPD_FILE, NAMES_UPD_FILE, TAXA_UPD_FILE = parseArgv(sys.argv[1:])

	print("Starting ... ")
	print("\nOpening nodes file:" + NODES_FILE)
	print("\nOpening names file:" + NAMES_FILE)
	print("\nOpening taxa file:" + TAXA_FILE)
	print("\nOpening output file:" + OUTPUT_FILE)	
	print("\nOpening nodes-updated file:" + NODES_UPD_FILE)	
	print("\nOpening names-updated file:" + NAMES_UPD_FILE)	
	print("\nOpening taxa-updated file:" + TAXA_UPD_FILE)	

	# read nodes
	dictNodes = {} # store all nodes and their parents
	readFile(NODES_FILE, dictNodes);   

	# read taxa: e.g. 10020	Dipodomys_ordiis
	dictTaxa = {}
	readFile(TAXA_FILE, dictTaxa);
	#print("dictTaxa: " + str(dictTaxa))

	# read all nodes' names; e.g. 8292	Amphibia
	dictNames = {}
	readFile(NAMES_FILE, dictNames);
	#print("dictNames: " + str(dictNames))	

	dictParentNodes = {} # store parental nodes and count the number of their children
	for nodeId in dictNodes:
		parentId = dictNodes[nodeId]
		if parentId in dictParentNodes:
			dictParentNodes[parentId] += 1
		else:
			dictParentNodes[parentId] = 1
			
	dictNodesTaxa = {} # store only taxa parental nodes who have one child
	for taxId in dictTaxa:
		taxParentId = dictNodes[taxId]
		if dictParentNodes[taxParentId] == 1:
			dictNodesTaxa[taxParentId] = taxId
	#print(dictNodesTaxa)
	
	# print nodeIds for the nodes which have only one child
	outFile = open(OUTPUT_FILE, "w")
	text = ""
	for nodeId in dictNodesTaxa:
		taxId = dictNodesTaxa[nodeId]
		text += nodeId + "\t" +  taxId + "\n"
		
		# update nodes and taxa names
		dictNodes[taxId] = dictNodes[nodeId] # UPD: taxId + grandParent_nodeId
		dictNodes.pop(nodeId) # DEL: parent_nodeId + grandParent_nodeId

		dictTaxa[taxId] += " (" + dictNames[nodeId] + ")" # UPD taxon's name --> needed for names-upd file
		dictNames[taxId] += " (" + dictNames[nodeId] + ")" # UPD taxon's name --> needed for taxa-upd file
		dictNames.pop(nodeId) # DEL parental node from names file
	# end for
	outFile.write(text) # print nodes with one child
	outFile.close()
	
	# update nodes
	writeFile(NODES_UPD_FILE, dictNodes)

	# update names
	writeFile(NAMES_UPD_FILE, dictNames)

	# update taxa
	writeFile(TAXA_UPD_FILE, dictTaxa)
		
	print("Done.")
# end if
