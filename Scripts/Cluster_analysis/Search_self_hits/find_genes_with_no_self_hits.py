__version__ = "1.0"
__author__ = "Korina Šimičević"

# changes made by MDL, 2023/2024  
import os, sys, getopt
import pandas as pd

def parseArgv(argv):
	print("Program arguments:")

	options = "hi:g:o:" # options that require arguments should be followed by :
	long_options = ["help", "input_file", "genes_file", "output_file"]
	try:
		opts, args = getopt.getopt(argv, options, long_options)
	except getopt.GetoptError:
		print("Error. Usage: find_genes_with_no_self_hits.py -i <input_file> -g <genes_file> -o <output_file>")
		sys.exit(2)
	
	INPUT_FILE = ""
	GENES_FILE = ""
	OUTPUT_FILE = ""	

	for opt, arg in opts:
		print(opt + "\t" + arg)
		if opt in ("-h", "--help"):
			print("Usage: find_genes_with_no_self_hits.py -i <input_file> -g <genes_file> -o <output_file>")
			sys.exit()
		elif opt in ("-i", "--input_file"):
			INPUT_FILE = arg
		elif opt in ("-g", "--genes_file"):
			GENES_FILE = arg
		elif opt in ("-o", "--output_file"):
			OUTPUT_FILE = arg
	# end for
	return INPUT_FILE, GENES_FILE, OUTPUT_FILE
# end parseArgv

# main
if __name__ == "__main__":
	INPUT_FILE, GENES_FILE, OUTPUT_FILE = parseArgv(sys.argv[1:])
	
	# genes with self-hits
	hits_df = pd.read_csv(INPUT_FILE, names=['query', 'target', 'eval'], sep='\t')
	hits_df['query'] = hits_df['query'].map(str.strip)
	
	# all genes
	genes_df = pd.read_csv(GENES_FILE, names=['gene'], sep='\t')
	genes_df['gene'] = genes_df['gene'].map(str.strip)
	
	with open(OUTPUT_FILE, 'w') as out: # without self-hits
		out.write(genes_df[~genes_df['gene'].isin(hits_df['query'])].to_csv(sep='\t', header=False, index=False))
		out.flush()
