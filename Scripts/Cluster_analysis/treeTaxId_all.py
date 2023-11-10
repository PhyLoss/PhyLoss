__version__ = "1.0"
__author__ = "Mirjana Domazet-Loso"

#
# treeTaxId_all.py - get the phylogeny tree for a taxId (print all parents up to the root)
#	--> Run: time python3 treeTaxId_all.py -i PATH_TO_FOLDER/names_All_2022.txt -j PATH_TO_FOLDER/nodes_All_2022.txt -t PATH_TO_FOLDER/taxIDNames_All_2022.txt -o PATH_TO_FOLDER/Parents/ > out.txt     
#
import os, sys, getopt

def parseArgv(argv):
    print("Program arguments:")

    options = "hi:j:t:o:" # options that require arguments should be followed by :
    long_options = ["help", "in_file_names", "in_file_nodes", "in_file_taxId", "out_file_parents"]
    try:
        opts, args = getopt.getopt(argv, options, long_options)
    except getopt.GetoptError:
        print("Error. Usage: treeTaxId_all.py -i <in_file_names> -j <in_file_nodes> -t <in_file_taxId> -o <out_file_parents>")
        sys.exit(2)
    
    INPUT_FILE_NAMES = "" 
    INPUT_FILE_NODES = ""
    INPUT_FILE_TAXID = ""
    OUTPUT_FOLDER = ""
	

    for opt, arg in opts:
        print(opt + "\t" + arg)
        if opt in ("-h", "--help"):
            print("Usage: treeTaxId_all.py -i <in_file_names> -j <in_file_nodes> -t <in_file_taxId> -o <out_file_parents>")
            sys.exit()
        elif opt in ("-i", "--in_file_names"):
            INPUT_FILE_NAMES = arg
        elif opt in ("-j", "--in_file_nodes"):
            INPUT_FILE_NODES = arg
        elif opt in ("-t", "--in_file_taxId"):
            INPUT_FILE_TAXID = arg
        elif opt in ("-o", "--out_file_parents"):
            OUTPUT_FOLDER = arg
    # end for
    return INPUT_FILE_NAMES, INPUT_FILE_NODES, INPUT_FILE_TAXID, OUTPUT_FOLDER
# end parseArgv

# main
if __name__ == "__main__":
    INPUT_FILE_NAMES, INPUT_FILE_NODES, INPUT_FILE_TAXID, OUTPUT_FOLDER = parseArgv(sys.argv[1:])
    
    # read nodes
    dictNodes = {}
    inFile = open(INPUT_FILE_NODES, "r")
    for line in inFile: # e.g. 526226    2
        arr = line.strip().split("\t")
        dictNodes[arr[0]] = arr[1]
    #print(dictNodes)
    inFile.close()

    # read names
    dictNames = {}
    inFile = open(INPUT_FILE_NAMES, "r")
    for line in inFile: # e.g. 526226    Gordonia_bronchialis_dsm_43247
        arr = line.strip().split("\t")
        dictNames[arr[0]] = arr[1]
    #print(dictNames)
    inFile.close()
    dictNames[str(131567)] = "Cellular_organisms" # manually add
	
    # read taxId
    dictTaxId = {}
    inFile = open(INPUT_FILE_TAXID, "r")
    for line in inFile: # e.g. 526226    Gordonia_bronchialis_dsm_43247
        arr = line.strip().split("\t")
        dictTaxId[arr[0]] = arr[1]
    #print(dictTaxId)
    inFile.close()

    for taxId in dictTaxId: # for each taxId: create a new file with all parents up to the root
        OUTPUT_FILE = OUTPUT_FOLDER + taxId + "parents.txt"
        outFile = open(OUTPUT_FILE, "w")
        t = taxId
        text = ""
        while True: # get all parents
            if t in dictNodes: # child: t
                # find taxId + name
                text += str(t) + "\t" + dictNames[t] + "\n"
                t = dictNodes[t]; # new child        
            if t == "131567": # stop when the root (Cellular_organisms) in encountered
                text += str(t) + "\t" + dictNames[t] + "\n"
                break
        outFile.write(text)
        outFile.close()

