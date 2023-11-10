__version__ = "1.0"
__author__ = "Mirjana Domazet-Loso"

#
# convertAll.py - convert all protein sequence files to our standard faa format using clean.py
#       
# 
#	--> Run: time python3 convertAll.py -t taxId.txt -i PATH_TO_INPUT_FOLDER [-r PATH_TO_REPORT_FOLDER] -o PATH_TO_OUTPUT_FOLDER

import os, sys, getopt
from os import listdir
from os.path import isfile, join

def parseArgv(argv):
	print("Program arguments:")

	options = "hi:t:r:o:" # options that require arguments should be followed by :
	long_options = ["help", "input_folder", "taxId_file", "report_folder", "output_folder"]
	try:
		opts, args = getopt.getopt(argv, options, long_options)
	except getopt.GetoptError:
		print("Error. Usage: converAll.py -i <input_folder> -t <taxId_file> [-r <report_folder>] -o <output_folder>")
		sys.exit(2)
	
	INPUT_FOLDER = ""
	REPORT_FOLDER = ""
	TAXID_FILE = ""
	OUTPUT_FOLDER = ""	

	for opt, arg in opts:
		print(opt + "\t" + arg)
		if opt in ("-h", "--help"):
			print("Usage: converAll.py -i <input_folder> -t <taxId_file> [-r <report_folder>] -o <output_folder>")
			sys.exit()
		elif opt in ("-i", "--input_folder"):
			INPUT_FOLDER = arg
		elif opt in ("-r", "--report_folder"):
			REPORT_FOLDER = arg
		elif opt in ("-t", "--taxId_file"):
			TAXID_FILE = arg
		elif opt in ("-o", "--output_folder"):
			OUTPUT_FOLDER = arg
	# end for
	return INPUT_FOLDER, REPORT_FOLDER, TAXID_FILE, OUTPUT_FOLDER
# end parseArgv

# main
if __name__ == "__main__":
	INPUT_FOLDER, REPORT_FOLDER, TAXID_FILE, OUTPUT_FOLDER = parseArgv(sys.argv[1:])
	
	taxIdFile = open(TAXID_FILE, "r")
	dictTaxa = {}
	# read taxa id + names
	for line in taxIdFile:
		arr = line.strip().split() # e.g. 9606	Homo_sapiens
		dictTaxa[arr[1]] = arr[0]
	# end for
	taxIdFile.close()

	# scan all files in INPUT_FOLDER
	files = [f for f in listdir(INPUT_FOLDER) if isfile(join(INPUT_FOLDER, f))]		
	for i in range(len(files)):
		# get first part of the file name (up to the first ".")
		start = 0
		end = files[i].find(".")
		name = files[i][start:end]
		taxId = dictTaxa[name]
		# input file
		inputFile = join(INPUT_FOLDER, files[i])
		
		# report file
		if REPORT_FOLDER != "":
			fileReport = join(REPORT_FOLDER, taxId + "_report.txt")
			commandReport = " -r " + str(fileReport)
		else:
			commandReport = ""

		# call clean.py with arguments
		#time python3 clean.py -t 9606 -i PATH_TO_INPUT_FOLDER/Homo_sapiens.fa -r PATH_TO_FOLDER_REPORT/report.txt -o PATH_TO_FOLDER_EDITED_FASTA
		command = "python clean.py -t " + taxId + " -i " + str(inputFile) + commandReport + " -o " + OUTPUT_FOLDER
		os.system(command)
		if (i + 1) % 100 == 0:
			print("Analyzed " + str(i+1) + " files.\n", flush=True)
		#sys.exit()
	# end for
	print("Analyzed " + str(i+1) + " files.\n", flush=True)
	print("Done.")
