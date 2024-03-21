__version__ = "1.0"
__author__ = "Korina Šimičević"

# changes made by MDL, 2023/2024  
import os, sys, getopt
import pandas as pd

def parseArgv(argv):
	print("Program arguments:")

	options = "hi:o:" # options that require arguments should be followed by :
	long_options = ["help", "input_file", "output_file"]
	try:
		opts, args = getopt.getopt(argv, options, long_options)
	except getopt.GetoptError:
		print("Error. Usage: get_self_hits.py -i <input_file> -o <output_file>")
		sys.exit(2)
	
	INPUT_FILE = ""
	OUTPUT_FILE = ""	

	for opt, arg in opts:
		print(opt + "\t" + arg)
		if opt in ("-h", "--help"):
			print("Usage: get_self_hits.py -i <input_file> -o <output_file>")
			sys.exit()
		elif opt in ("-i", "--input_file"):
			INPUT_FILE = arg
		elif opt in ("-o", "--output_file"):
			OUTPUT_FILE = arg
	# end for
	return INPUT_FILE, OUTPUT_FILE
# end parseArgv

# main
if __name__ == "__main__":
	INPUT_FILE, OUTPUT_FILE = parseArgv(sys.argv[1:])
	
	#tmp_path = OUT_PATH + "/hits_" + str(taxId) + ".txt"
	tmp_path = INPUT_FILE
	#print(tmp_path)
	#hits_path = OUT_PATH + "/s_hits_" + str(taxId) + ".txt"
	hits_path = OUTPUT_FILE
	
	df = pd.read_csv(tmp_path, delimiter='\t', names=['query', 'target', 'e-value'])
	df = df[df['query'] == df['target']]
	
	with open(hits_path, 'a') as hitf:
		hitf.write(df.to_csv(sep='\t', header=False, index=False))
		hitf.flush()
