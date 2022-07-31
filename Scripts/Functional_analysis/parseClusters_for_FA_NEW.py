__version__ = "1.0"
__author__ = "Mirjana Domazet-Loso"

"""
parseClusters_for_FA_NEW.py - parse mmseqs cluster output (e.g. for c = 0.8) to get a shorten cluster representation used in functional analysis
	--> Run: time python3 parseClusters_for_FA_NEW.py -i PATH_TO_ANALYSIS_FOLDER/results_0_8/Hsap/dbAllPlus_clu_0_8_Hsap_out.txt -o ./Clusters_2022/c_0_8/Hsap_2022_clusters.txt -p 39 > parseClu_for_FA_NEW.out    
"""


import sys, getopt
PS_MIN = int(1)

def parseArgv(argv):
    print("Program arguments:")

    options = "hi:o:p:" # options that require arguments should be followed by :
    long_options = ["help", "ifile", "ofile", "psmax"]
    try:
        opts, args = getopt.getopt(argv, options, long_options)
    except getopt.GetoptError:
        print("Error. Usage: parseClusters_for_FA_NEW.py -i <inputfile> -o <outputfile> -p <psmax>")
        sys.exit(2)
    
    INPUT_FILE = "" 
    OUTPUT_FILE = "" 

    for opt, arg in opts:
        print(opt + "\t" + arg)
        if opt in ("-h", "--help"):
            print("Usage: parseClusters_for_FA_NEW.py -i <inputfile> -o <outputfile> -p <psmax>")
            sys.exit()
        elif opt in ("-i", "--ifile"):
            INPUT_FILE = arg
        elif opt in ("-o", "--ofile"):
            OUTPUT_FILE = arg
        elif opt in ("-p", "--psmax"):
            PS_MAX = int(arg)
    # end for
    return INPUT_FILE, OUTPUT_FILE, PS_MAX
# end parseArgv

def constructOutput(inFile, outFile, PS_MAX):
    geneList = {}
    firstCluster = True
    
    for line in inFile: 
        if line.strip(): # check only non-empty lines         
            if line.find("Cluster representative") != -1:
                if firstCluster == True:
                    firstCluster = False
                else:
                    # resolve the previous cluster, i.e. print cluster information to the output
                    if psMin != psMax or psMax == PS_MAX: # cluster with only species from the same phylostratum (and no focal species)
                        if psMax == PS_MAX: psLoss = '-' 
                        else: psLoss = psMax + 1
                        firstRow += str(psMin) + "\t" + str(psLoss) + "\t" + numGenes
                        outFile.write(firstRow + "\n")
                        for gene in geneList: # will this keep the same order as in the input file?
                            #outFile.write(gene + "\t" + geneList[gene] + "\n") # geneID, PS
                            outFile.write(gene + "\n") # only geneID
                    
                geneList.clear()
                cntGenes = 0
                psMin = PS_MAX
                psMax = PS_MIN

                clusterReprLine = line.split();
                firstRow = "\n" + clusterReprLine[0] + "\t" + clusterReprLine[3] + "\t";
                temp = clusterReprLine[3].split('|')
                firstRow += temp[3] + "\t"  # add taxid (temp[3])
                numGenes = clusterReprLine[10]
                #print(firstRow)
                for line in inFile: # read genes' lines (geneId + ps)
                    if line.find("pid") != -1: 
                        geneArr = line.split()
                        line2 = inFile.readline()
                        geneArr2 = line2.split()
                        geneList[geneArr[1]] = geneArr2[1]
                        #print(geneList)
                        if int(geneArr2[1]) < psMin: psMin = int(geneArr2[1])
                        if int(geneArr2[1]) > psMax: psMax = int(geneArr2[1])
                        cntGenes += 1
                        if cntGenes == int(numGenes):
                            break
                # end for
    # end for

    # resolve the last cluster
    if psMin != psMax or psMax == PS_MAX: # cluster with only species from the same phylostratum (and no focal species)    
        if psMax == PS_MAX: psLoss = '-' 
        else: psLoss = psMax + 1
        firstRow += str(psMin) + "\t" + str(psLoss) + "\t" + numGenes
        outFile.write(firstRow + "\n")
        for gene in geneList: # will this keep the same order as in the input file?
            outFile.write(gene + "\n") # only geneID
# end constructOutput


# main
if __name__ == "__main__":
    INPUT_FILE, OUTPUT_FILE, PS_MAX = parseArgv(sys.argv[1:])

    print("Starting ... ")
    print("\nOpening input file:" + INPUT_FILE)
    inFile = open(INPUT_FILE, 'r')
    print("\nOpening output file:" + OUTPUT_FILE)
    outFile = open(OUTPUT_FILE, 'w')
    constructOutput(inFile, outFile, PS_MAX)

    inFile.close()
    outFile.close()
    print("Done.")
# end if
