__version__ = "1.0"
__author__ = "Korina Šimičević"

#modifications MDL; 2023/2024
import os, sys, getopt

def parseArgv(argv):
	print("Program arguments:")

	options = "hd:o:" # options that require arguments should be followed by :
	long_options = ["help", "db_folder", "output_file"]
	try:
		opts, args = getopt.getopt(argv, options, long_options)
	except getopt.GetoptError:
		print("Error. Usage: find_all_genes.py -d <db_folder> -o <output_file>")
		sys.exit(2)
	
	DB_FOLDER = ""
	OUTPUT_FILE = ""	

	for opt, arg in opts:
		print(opt + "\t" + arg)
		if opt in ("-h", "--help"):
			print("Usage: find_all_genes.py -d <db_folder> -o <output_file>")
			sys.exit()
		elif opt in ("-d", "--db_folder"):
			DB_FOLDER = arg
		elif opt in ("-o", "--output_file"):
			OUTPUT_FILE = arg
	# end for
	return DB_FOLDER, OUTPUT_FILE
# end parseArgv

# main
if __name__ == "__main__":
	DB_FOLDER, OUTPUT_FILE = parseArgv(sys.argv[1:])
	with open(OUTPUT_FILE, 'w') as genes_out:
		for file in os.listdir(DB_FOLDER):
			with open(DB_FOLDER + "/" + file, 'r') as taxon_file:
				lines = taxon_file.readlines()
				for line in lines:
					if line.startswith(">"):
						parts = line.split()
						genes_out.write(f"{parts[0][1:]}\n")
