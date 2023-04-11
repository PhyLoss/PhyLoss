__version__ = "1.0"
__author__ = "Mirjana Domazet-Loso"

"""
    getAllLCA.py - get LCA of all pairs of taxa (taxId1, taxId2)
    Usage: time python3 getAllLCA.py -i PATH_TO_FOLDER/Data/allParents.txt -o PATH_TO_FOLDER/Data/allLCA.txt    
"""


import os, sys, getopt

def parseArgv(argv):
    print("Program arguments:")

    options = "hi:o:" # options that require arguments should be followed by :
    long_options = ["help", "in_file", "out_file"]
    try:
        opts, args = getopt.getopt(argv, options, long_options)
    except getopt.GetoptError:
        print("Error. Usage: getAllLCA.py -i <in_file> -o <out_file>")
        sys.exit(2)
    
    INPUT_FILE = "" 
    OUTPUT_FILE = ""
	

    for opt, arg in opts:
        print(opt + "\t" + arg)
        if opt in ("-h", "--help"):
            print("Usage: getAllLCA.py -i <in_file> -o <out_file>")
            sys.exit()
        elif opt in ("-i", "--in_file"):
            INPUT_FILE = arg
        elif opt in ("-o", "--out_file"):
            OUTPUT_FILE = arg
    # end for
    return INPUT_FILE, OUTPUT_FILE
# end parseArgv

# main
if __name__ == "__main__":
    INPUT_FILE, OUTPUT_FILE = parseArgv(sys.argv[1:])
    
    # read all taxa from input file
    # for each taxon read its parents up to "Cellular_organisms"
    inFile = open(INPUT_FILE, "r")
    dictTaxa = {}
    dictParents = {} # dictionary of lists
    newTaxon = True
    for line in inFile: # e.g. 10020	Dipodomys_ordii
        arr = line.strip().split("\t")
        if newTaxon:            
            taxId = arr[0]
            dictTaxa[taxId] = arr[1]
            newTaxon = False
            dictParents[taxId] = []
        else: # parental nodes
            dictParents[taxId].append(arr[0])
            if arr[0] == "131567": # Cellular_organisms
                newTaxon = True # get ready for the next taxon

    inFile.close()

    # compute LCA for each pair of taxa
    LCA = {}
    for t in dictTaxa: # for each taxon find its LCA with all other taxa t2 != t
        for t2 in dictTaxa:
            if t == t2:
                LCA[(t, t2)]= -1
            else:
                found = False; # flag for found LCA
                for i in range(len(dictParents[t])): # check all t's parents to find LCA with t2's parents
                    for j in range(len(dictParents[t2])):
                        if dictParents[t][i] == dictParents[t2][j]:
                            found = True
                            break
                    if found: 
                        LCA[(t, t2)] = dictParents[t][i] # or: dictParents[t2][j]
                        break            
    

    # print all LCA
    outFile = open(OUTPUT_FILE, "w")
    
    text = 23 * ' ' # 23 spaces at the beginning of the header row    
    for t in dictTaxa: # taxID
        text += ' ' + t.rjust(10, " ") # print taxId on 11 positions
    outFile.write(text + "\n")
    
    text = 23 * ' ' # 23 spaces at the beginning of the taxa names' row    
    for t in dictTaxa: # taxname
        text += ' ' + dictTaxa[t][:10].rjust(10, " ")
    outFile.write(text + "\n")

    # print rows starting with taxa names and then followed by LCA values
    for t2 in dictTaxa: 
        text = t2.rjust(10, " ") + " (" + dictTaxa[t2][:10].rjust(10, " ") + ") "
        for t in dictTaxa:
            text += str(LCA[(t2, t)]).rjust(10, " ") + " "
        outFile.write(text + "\n")
    # end for
    outFile.close()
