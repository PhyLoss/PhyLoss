__version__ = "1.0"
__author__ = "Korina Šimičević"

# changes made by MDL, 2023/2024  
# --> originally writing updated DB was too slow; 
# --> modified so that each db file is open and read/written only once

import os, sys, getopt

def parseArgv(argv):
	print("Program arguments:")

	options = "hi:d:" # options that require arguments should be followed by :
	long_options = ["help", "db_folder", "input_file"]
	try:
		opts, args = getopt.getopt(argv, options, long_options)
	except getopt.GetoptError:
		print("Error. Usage: filter_db.py -i <input_file> -d <db_folder>")
		sys.exit(2)
	
	DB_FOLDER = ""
	INPUT_FILE = ""	

	for opt, arg in opts:
		print(opt + "\t" + arg)
		if opt in ("-h", "--help"):
			print("Usage: filter_db.py -i <input_file> -d <db_folder>")
			sys.exit()
		elif opt in ("-d", "--db_folder"):
			DB_FOLDER = arg
		elif opt in ("-i", "--input_file"):
			INPUT_FILE = arg
	# end for
	return DB_FOLDER, INPUT_FILE
# end parseArgv

# auxiliary functions - start
def getTaxId(geneId):
	arr = geneId.split("|")
	return arr[3]

def get_pid(geneId):
	arr = geneId.split("|")
	return int(arr[1])

def extract_protein_from_fasta(line: str) -> str:
	return line[1:].split("\t")[0]

# auxiliary functions - end


def filter_db(db_path: str, filter_path: str) -> None:
	filter = [p.strip() for p in open(filter_path, 'r').readlines()]

	dictTaxa = {} # dictionary of taxa with self missing hits
	for protein in filter:
		taxId = getTaxId(protein)
		if taxId not in dictTaxa: # add protein
			dictTaxa[taxId] = {}
		# add protein to taxId's dictionary -- these proteins do not have self-hits
		geneId = get_pid(protein)
		dictTaxa[taxId][geneId] = 1

	# MDL, March 2024
	#for genomes that have genes without self-hits
	for taxId in dictTaxa:
		protein_file = db_path + taxId + ".faa"
		lines = open(protein_file, 'r').readlines()
		
		allLines = ""
		cnt = 0
		
		#with open(protein_file, 'w') as db_out:
		db_out = open(protein_file, 'w')
		print(protein_file)
		for line in lines:
			if line.startswith(">"):
				current_protein = extract_protein_from_fasta(line)

				if get_pid(current_protein) not in dictTaxa[taxId]:
					write = True
					cnt += 1
					if cnt % 1000 == 0: # print in batches 
						db_out.write(allLines)
						db_out.flush()
						allLines = ""
				else:
					write = False
			if write:
				allLines += line			
		# end for
		db_out.write(allLines)
		db_out.close()		
# end filter_db
		
if __name__ == '__main__':
	DB_FOLDER, INPUT_FILE = parseArgv(sys.argv[1:])
	filter_db(DB_FOLDER, INPUT_FILE)

