__version__ = "1.0"
__author__ = "Mirjana Domazet-Loso"

#
# clean.py - clean proteins sequences
#	--> in case of multiple splicing variants, keep the protein with the longest sequence, and other variants (proteins with the same name) delete
#	--> delete proteins shorter than MIN_LEN and longer than MAX_LEN a.a.
#	--> delete proteins containing "*" in the middle of the sequence
#	--> Later: (after the mmseqs search has been done) delete proteins from the input list no_self_hits.txt (i.e. protein seq-s with no self-hits obtained using mmseqs search)

#	--> Run: python3 clean.py -t 9606 -i PATH_TO_INPUT_FOLDER/Homo_sapiens.fa -r PATH_TO_FOLDER_REPORT/report.txt -o PATH_TO_FOLDER_EDITED_FASTA
# 
import os, sys, getopt

MIN_LEN = 10
MAX_LEN = 10000

def parseArgv(argv):
	print("Program arguments:")

	options = "hi:r:t:o:" # options that require arguments should be followed by :
	long_options = ["help", "input_file_name", "report_file_name", "input_file_taxId", "output_folder"]
	try:
		opts, args = getopt.getopt(argv, options, long_options)
	except getopt.GetoptError:
		print("Error. Usage: clean.py -i <input_file_name> [-r <report_file_name>] -t <input_file_taxId> -o <output_folder>")
		sys.exit(2)
	
	INPUT_FILE_NAME = ""
	REPORT_FILE_NAME = ""
	INPUT_FILE_TAXID = ""
	OUTPUT_FOLDER = ""	

	for opt, arg in opts:
		print(opt + "\t" + arg)
		if opt in ("-h", "--help"):
			print("Usage: clean.py -i <input_file_name> -t <input_file_taxId> [-r <report_file_name>] -o <output_folder>")
			sys.exit()
		elif opt in ("-i", "--input_file_name"):
			INPUT_FILE_NAME = arg
		elif opt in ("-r", "--report_file_name"):
			REPORT_FILE_NAME = arg
		elif opt in ("-t", "--input_file_taxId"):
			INPUT_FILE_TAXID = arg
		elif opt in ("-o", "--output_folder"):
			OUTPUT_FOLDER = arg
	# end for
	return INPUT_FILE_NAME, REPORT_FILE_NAME, INPUT_FILE_TAXID, OUTPUT_FOLDER
# end parseArgv

# auxilary function - print geneId-s from the list to a report file
def printList(rFile, headerLine, listGeneId):

	rFile.write(headerLine + "\n")
	text = ""
	for i in range(len(listGeneId)):
		text += listGeneId[i] + "\n"
	# end for
	if text != "":
		rFile.write(text)
	rFile.write("\n")

# end printList

# auxiliary function - resolve the previous line
def resolvePrevGene(seqLen, listShort, listLong, listStop, geneId, geneHeader, dictSeq, dictMul, listGeneId, seq):
	# delete proteins shorter than MIN_LEN and longer than MAX_LEN a.a.
	# or containing "*" in the middle of the sequence
	if seqLen <= MIN_LEN:
		listShort.append(geneHeader + "\t" + geneId)				
	elif seqLen >= MAX_LEN:
		listLong.append(geneHeader + "\t" + geneId)
	elif "*" in seq[:(len(seq)- 2)]: # stop codon found before the last char, -1 to skip "\n", since each line ends with "\n"
		listStop.append(geneHeader + "\t" + geneId)

	#if not (seqLen <= MIN_LEN or seqLen >= MAX_LEN or "*" in seq[:(len(seq)- 1)]):
	else:
		# in case of multiple splicing variants, keep the protein with the longest sequence, and other variants (proteins with the same name) delete
		#if geneId in dictSeq and len(dictSeq[geneId][1]) < len(seq):
		addGeneId = False
		if geneId in dictSeq:
			if geneId in dictMul: # count the number of splicing variants
				dictMul[geneId] += 1 # found 3rd or later splicing variant
			else:
				dictMul[geneId] = 2 # found 2nd splicing variant
						
			if len(dictSeq[geneId][1]) < len(seq): # found the longest splicing variant so far --> delete the previous one
				listGeneId.remove(geneId)
				addGeneId = True
		
		if geneId not in dictSeq or addGeneId:
			listGeneId.append(geneId) # add the first or the longest variant according to its position in the original file
			dictSeq[geneId] = (geneHeader, seq)			
# end resolvePrevGene

# read original fasta file and construct list of geneId-s and a dict with its sequences and gene headers
def readInputFile(inFile, reportFileName, listGeneId, dictSeq):
	
	firstLine = True
	listShort = []
	listLong = []
	dictMul = {}
	listStop = []
	seq = ""
	seqLen = 0
	geneId = ""
	geneHeader = ""

	# e.g. extract ENSG00000283951.1 from a header: ">ENSP00000491436.1 pep chromosome:GRCh38:CHR_HSCHR19KIR_7191059-1_CTG3_1:54850452:54867240:1 gene:ENSG00000283951.1 transcript:ENST00000640782.1 gene_biotype:protein_coding ..."
	for line in inFile: # e.g. 526226	2
		if ">" in line: # gene name
			# resolve the previous protein sequence (if exists)
			resolvePrevGene(seqLen, listShort, listLong, listStop, geneId, geneHeader, dictSeq, dictMul, listGeneId, seq)
				
			# a new protein seq.
			startPos = line.find("gene:") + 5
			endPos = line.find("transcript:") - 1
			if startPos >= 0 and endPos >= startPos:
				geneId = line[startPos:endPos] # get all characters from start to end - 1
				# Note: there can be >= 1 genes with the same geneId, but different gene headers (ENSP... number in the above example)
			else:
				# in case geneId is empty string, then set geneID to the data in the 1st column; 
				# this can happen, if the file's source is NCBI, since they don't have the above naming convention
				geneId = line.strip().split()[0][1:] # skip ">"
			
			geneHeader = line.strip().split()[0][1:] # sometimes the same as geneId
			seq = ""
			seqLen = 0
			firstLine = False
		
		else: # sequence line(s)
			seq += line.strip() + "\n" # proteinFullId = sprintf(">pid|%019d|tx|%d|\t%s", proteinId, taxId, substr($1, 2));
			seqLen += len(line.strip())
	# end for

	# resolve the last protein sequence
	resolvePrevGene(seqLen, listShort, listLong, listStop, geneId, geneHeader, dictSeq, dictMul, listGeneId, seq)

	# store information on removed protein sequences in the report file (if it exists)
	if reportFileName != "":
		rFile = open(reportFileName, "w")
		printList(rFile, "Too short protein sequences (<= 10 a.a.):", listShort)
		printList(rFile, "Too long protein sequences (>= 10000 a.a.):", listLong)
		printList(rFile, "Protein sequences with the stop codon in the middle of the seq.:", listStop)

		# print geneId-s with multiple splicing variants
		rFile.write("Multiple splicing variants:\n")
		text = ""
		for geneId in dictMul:
			text += geneId + "\t" + str(dictMul[geneId]) + "\n"
		if text != "": 
			rFile.write(text)
		rFile.close()
# end readInputFile


# main
if __name__ == "__main__":
	INPUT_FILE_NAME, REPORT_FILE_NAME, INPUT_FILE_TAXID, OUTPUT_FOLDER = parseArgv(sys.argv[1:])

	# read fasta file
	dictSeq = {}
	listGeneId = []
	inFile = open(INPUT_FILE_NAME, "r")
	reportFileName = REPORT_FILE_NAME
	readInputFile(inFile, REPORT_FILE_NAME, listGeneId, dictSeq)
	inFile.close()
	
	# print new geneIds, old geneIds and seq-s (in case of multiple splicing variants, keep the protein with the longest sequence)
	taxId = INPUT_FILE_TAXID
	outFile = open(OUTPUT_FOLDER + "/" + taxId + ".faa", "w")
	
	text = ""
	for i in range(len(listGeneId)):
		# sprintf(">pid|%019d|tx|%d|", proteinId, taxId);
		geneId = listGeneId[i]
		geneId_i = int(taxId + str(i + 1)); # starting with 1, ...
		newGeneId = ">pid|{:019d}|tx|{:d}|".format(geneId_i, int(taxId)) # sign-aware zero-padding
		(geneHeader, seq) = dictSeq[geneId]
		text += newGeneId + "\t" + geneHeader + "\n" + seq
		
		if i % 500 == 0: # print in batches
			outFile.write(text)
			text = ""
	# end for
	outFile.write(text) # print last batch
	outFile.close()



