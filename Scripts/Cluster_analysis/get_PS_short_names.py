__version__ = "1.0"
__author__ = "Mirjana Domazet-Loso"
#2021-2023
#
# get_PS_short_names.py 	- read parents files from an input folder; 
#							- for each taxon: 
#							construct list of PS names starting with (1) cellular org. (use existing abbrev. for PS 1 to 4)
#							
#	--> Run: time python3 get_PS_short_names.py -i PATH_TO_INPUT_FILE -p PATH_TO_PARENTS_FOLDER -o PATH_TO_FOLDER_PS_NAMES > out.txt     
#
import os, sys, getopt

def parseArgv(argv):
    print("Program arguments:")

    options = "hi:p:o:" # options that require arguments should be followed by :
    long_options = ["help", "in_file", "parents_folder", "out_folder"]
    try:
        opts, args = getopt.getopt(argv, options, long_options)
    except getopt.GetoptError:
        print("Error. Usage: get_PS_short_names.py -i <in_file> -p <parents_folder> -o <out_folder>")
        sys.exit(2)
    
    INPUT_FILE = ""
    PARENTS_FOLDER = ""
    OUTPUT_FOLDER = "" 

    for opt, arg in opts:
        print(opt + "\t" + arg)
        if opt in ("-h", "--help"):
            print("Usage: get_PS_short_names.py -i <in_file> -p <parents_folder>-o <out_folder>")
            sys.exit()
        elif opt in ("-i", "--in_file"):
            INPUT_FILE = arg
        elif opt in ("-p", "--parents_folder"):
            PARENTS_FOLDER = arg
        elif opt in ("-o", "--out_folder"):
            OUTPUT_FOLDER = arg
    # end for
    return INPUT_FILE, PARENTS_FOLDER, OUTPUT_FOLDER
# end parseArgv

# main
if __name__ == "__main__":
    INPUT_FILE, PARENTS_FOLDER, OUTPUT_FOLDER = parseArgv(sys.argv[1:])

    print("Starting ... ")
    print("\nOpening input file:" + INPUT_FILE)
    print("\nOpening parents folder:" + PARENTS_FOLDER)
    print("\nOpening output folder:" + OUTPUT_FOLDER)
    
    inFile = open(INPUT_FILE, "r")
    dictTaxa = {}
    for line in inFile: # e.g. 10020    Dipodomys_ordii    D_ordii_10020
        arr = line.strip().split("\t")
        dictTaxa[arr[0]] = arr[1].split(" ")[0] + "_" + arr[0].strip() # short taxon's name
        #dictTaxa[arr[0]] = arr[2] # taxon's name + taxId
    inFile.close()

    #3000001    Archaea1    DPANN/Euryarchaeota/TACK/Asgard_archaea/Eukaryota
    #3000002    Archaea2    Euryarchaeota/TACK/Asgard_archaea/Eukaryota
    #3000003    Archaea3    TACK/Asgard_archaea/Eukaryota
    #3000004    Asgard/Eukaryota    Asgard_archaea/Eukaryota
    dictAbbrev = {}
    dictAbbrev["3000001"] = "3000001" + "\t" + "Archaea1" + "\t" + "DPANN/Euryarchaeota/TACK/Asgard_archaea/Eukaryota"
    dictAbbrev["3000002"] = "3000002" + "\t" + "Archaea2" + "\t" + "Euryarchaeota/TACK/Asgard_archaea/Eukaryota"
    dictAbbrev["3000003"] = "3000003" + "\t" + "Archaea3" + "\t" + "TACK/Asgard_archaea/Eukaryota"
    dictAbbrev["3000004"] = "3000004" + "\t" + "Asgard/Eukaryota" + "\t" + "Asgard_archaea/Eukaryota"
    

    # new, 2023
    dictAbbrev["3000216"] = "3000216" + "\t" + "Phaff./Saccharomycetaceae" + "\t" + "Phaffomycetaceae/Saccharomycetaceae"
    dictAbbrev["3000217"] = "3000217" + "\t" + "Saccharomyco./mycetaceae" + "\t" + "Saccharomycodaceae/mycetaceae"
    dictAbbrev["3000223"] = "3000217" + "\t" + "Kaz./Nau./Saccharomyces" + "\t" + "Kazachstania/Naumovozyma/Saccharomyces"


# 5000125	Bibiono./Cyclorrhapha	Bibionomorpha/Cyclorrhapha
# 5000100	Choano./Metazoa	Choanoflagellida/Metazoa
    dictAbbrev["5000125"] = "5000125" + "\t" + "Bibiono./Cyclorrhapha" + "\t" + "Bibionomorpha/Cyclorrhapha"
    dictAbbrev["5000100"] = "5000100" + "\t" + "Choano./Metazoa" + "\t" + "Choanoflagellida/Metazoa"


    dictAbbrev["7000116"] = "7000116" + "\t" + "Sap./Malv./Brassicales" + "\t" + "Sapindales/Malvales/Brassicales"

    dictAbbrev["6000007"] = "6000007" + "\t" + "D_melanogaster_group" + "\t" + "Drosophila_melanogaster_group"
    dictAbbrev["6000009"] = "6000009" + "\t" + "D_melanogaster_subgroup" + "\t" + "Drosophila_melanogaster_subgroup"
    dictAbbrev["6000011"] = "6000011" + "\t" + "sim./sech./mel. group" + "\t" + "Drosophila_similans/sechellia/melanogaster_group"


    filenames = os.listdir(PARENTS_FOLDER)
    for f in filenames:
        if f.endswith("parents.txt"):		
            taxId = f.split("parents.txt")[0]
            #print(f + "\t" + taxId)            
            
            if taxId in dictTaxa:
                print(f + "\t" + taxId)            
                pFile = open(PARENTS_FOLDER + f, "r")

                outFileName = dictTaxa[taxId] + ".txt"
                outFile = open(OUTPUT_FOLDER + outFileName, "w")
                outFile.write("ID	short_name	name\n") # print the last batch
                text = ""
                cntLine = 0
                for line in pFile: # read lines (e.g. 131567    Cellular_organisms    Cellular_organisms)
                    # since I want the reversed order, each new line is added to the beginning of the text
                    arr = line.strip().split("\t")
                    if arr[0] in dictAbbrev: # abbreviations
                        newLine = dictAbbrev[arr[0]] + "\n"
                    else:
                        newLine = arr[0] + "\t" + arr[1] + "\t" + arr[1] + "\n" # id, short_name, name
                    text = newLine + text
                # end for
                outFile.write(text) # print the last batch
                pFile.close()    
                outFile.close()    
                exit
            # end if taxId ...
    print("Done.")
# end if
