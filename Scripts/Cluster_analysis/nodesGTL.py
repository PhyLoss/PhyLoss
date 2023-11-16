"""
nodesGTL.py - compute GF total/gain/loss for all nodes (including leaves) in an input directory
    --> Input: input folder with subfolders containing GF total/gain/loss (+ parents folder, + names file)
    --> Output: a file for each c-value in an output folder
    --> Run: time python3 nodesGTL.py -i INPUT_FOLDER [-i INPUT_FOLDER2 ...] \
                                -o OUTPUT_FOLDER -p PARENTS_FOLDER -n NAMES_FILE
"""

__version__ = "1.0"
__author__ = "Mirjana Domazet-Lošo"
# November, 2023

import getopt
import sys
import os
from glob import glob

def parseArgv(argv):
    print("Program arguments:")

    options = "hi:o:p:n:" # options that require arguments should be followed by :
    long_options = ["help", "input-folder", "output-folder", "parents-folder", "names-file"]
    try:
        opts, args = getopt.getopt(argv, options, long_options)
    except getopt.GetoptError:
        print("Error. Usage: nodesGTL.py -i <input-folder> -o <output-folder> -p <parents-folder> -n <names-file>")
        sys.exit(2)
    
    INPUT_FOLDERS = [] # possibly multiple input folders
    OUTPUT_FOLDER = ""
    PARENTS_FOLDER = ""
    NAMES_FILE = ""

    for opt, arg in opts:
        print(opt + "\t" + arg)
        if opt in ("-h", "--help"):
            print("Usage: nodesGTL.py -i <input-folder> -o <output-folder> -p <parents-folder> -n <names-file>")
            sys.exit()
        elif opt in ("-i", "--input-folder"):
            INPUT_FOLDERS.append(arg)
        elif opt in ("-o", "--output-folder"):
            OUTPUT_FOLDER = arg
        elif opt in ("-p", "--parents-folder"):
            PARENTS_FOLDER = arg
        elif opt in ("-n", "--names-file"):
            NAMES_FILE = arg
    # end for
    print("\n")
    return INPUT_FOLDERS, OUTPUT_FOLDER, PARENTS_FOLDER, NAMES_FILE
# end parseArgv


# get names for each node (leaves also)
def getNames(NAMES_FILE):
    dictNames = {}
    with open(NAMES_FILE, "r") as f:
        for line in f: # e.g. 1924735	Salpingoeca_punica
            arr = line.strip().split()
            dictNames[arr[0]] = arr[1]
    return dictNames
# end getNames

# get all parents for each leaf/taxId
def getParents(PARENTS_FOLDER):
    dictParents = {}
    
    files = os.listdir(PARENTS_FOLDER)
    #print(files)
    for inFile in files:
        # e.g. 882parents.txt --> 882
        taxId = inFile[: inFile.find("parents.txt")]
        dictParents[taxId] = {}
        with open(PARENTS_FOLDER + "/" + inFile, "r") as f: # get the number of lines in file
            count = sum(1 for _ in f)
            
        with open(PARENTS_FOLDER + "/" + inFile, "r") as f:
            # 882	Desulfovibrio_vulgaris_str_hildenborough
            # 2	Bacteria ...
            i = 0;
            for line in f:
                i += 1
                dictParents[taxId][str(count - i + 1)] = line.strip().split()[0] # store nodeId  
    # end for inFile
    return dictParents
# end getParents

if __name__ == "__main__":
    print("__main__: Starting ... ")
    INPUT_FOLDERS, OUTPUT_FOLDER, PARENTS_FOLDER, NAMES_FILE = parseArgv(sys.argv[1:])
    
    dictParents = getParents(PARENTS_FOLDER)
    dictNames = getNames(NAMES_FILE)

    dictNodes = {}
    for c in range(0, 9):
        cFolder = "0_" + str(c)
        dictNodes[str(c)] = {}

    for path in INPUT_FOLDERS:
        print(path)
        for c in range(0, 9):
            cFolder = "0_" + str(c)
            files = os.listdir(path + "/" + cFolder)
            for inFile in files: # from c=0.8 to 0.0
                #summary_GL_Actinidia_chinensis_3625_0_0.txt --> 3625
                arrFile = inFile.split("_")
                taxId = arrFile[len(arrFile) - 3]
                with open(path + "/" + cFolder + "/" + inFile, "r") as f:
                    # summary_GL_Actinidia_chinensis_3625_0_0.txt
                    #ps	GF_Gained	GF_Lost
                    #1	8496	0
                    #2	204	0
                    total = totalG = totalL = 0
                    for line in f:
                        if "ps" in line:
                            continue
                        arr = line.strip().split()
                        ps = arr[0]
                        totalG += int(arr[1])
                        totalL += int(arr[2])
                        total = totalG - totalL
                        nodeId = dictParents[taxId][ps]
                        
                        if nodeId not in dictNodes[str(c)]:
                            # e.g. Thermococcus_kodakarensis_kod1,500,50,100
                            dictNodes[str(c)][nodeId] = [total, arr[1], arr[2]]
                #break
            # end for inFile
    
    
    #list of nodes that should not be printed, since they are not shown in the phylogeny
    listExclude = ["1841598", "1163730", "28889", "273075", "1737403"]
            
    # print results to output folder; 
    # for each c-value print a separate file
    for c in range(0, 9):
        outFileName = OUTPUT_FOLDER + "/" + "TGL_0_" + str(c) + ".txt"
        outFile = open(outFileName, "w")
        text = ""
        count = 0
        for nodeId in dictNodes[str(c)]:
            count += 1
            l = dictNodes[str(c)][nodeId]
            if nodeId in listExclude: 
                continue
            text += nodeId + "\t" + dictNames[nodeId] + "\t" \
                                    + str(l[0]) + "\t" + str(l[1]) + "\t"+ str(l[2]) + "\n"
            if count % 100 == 0:
                outFile.write(text)
                text = ""
        # end for
        # print the last batch
        outFile.write(text)
        outFile.close()
    # end for