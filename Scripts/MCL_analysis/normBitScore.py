__version__ = "1.0"
__author__ = "Mirjana Domazet-Loso"
# 2023

#
# normBitScore.py - normalize bit scores for each species pair (query, subject); 
#                   similar to Alg. on p. 10 of the paper "OrthoFinder: solving fundamental biases in whole genome comparisons dramatically 
# 															improves orthogroup inference accuracy" (Emms and Kelly, 2015);
#                   output: abc input to the mcl program

import os, sys, getopt
import math
import numpy as np

BINSIZE = 1000

def parseArgv(argv):
	print("Program arguments:")

	options = "hi:o:" # options that require arguments should be followed by :
	long_options = ["help", "input_file_name", "output_file_name"]
	try:
		opts, args = getopt.getopt(argv, options, long_options)
	except getopt.GetoptError:
		print("Error. Usage: normBitScore.py -i <input_file_name> -o <output_file_name>")
		sys.exit(2)
	
	INPUT_FILE_NAME = ""
	OUTPUT_FILE_NAME = ""	

	for opt, arg in opts:
		print(opt + "\t" + arg)
		if opt in ("-h", "--help"):
			print("Usage: normBitScore.py -i <input_file_name> -o <output_file_name>")
			sys.exit()
		elif opt in ("-i", "--input_file_name"):
			INPUT_FILE_NAME = arg
		elif opt in ("-o", "--output_file_name"):
			OUTPUT_FILE_NAME = arg
	# end for
	return INPUT_FILE_NAME, OUTPUT_FILE_NAME
# end parseArgv

# compute and print normalised bit-score for each hit (query_pid, subject_pid)
def printNormBitScore(dictQS, outFile, q, s):
	#	count rows for q, h
	#	sort by L_q_h
	#	split rows in bins of size = 10^3 (or 10^4??) according to L_q_h
	#	for each bin
	#		select the top 5% S'(i.e. B_q_h) 
	#		compute a, b using lin. regression (least squares) in log10(B_q_h) = a x log10(L_q_h) + b
	#		normalise each B_q_h --> B_q_h' = B_q_h / (10^b x L_q_h^a)
	lenDict = len(dictQS)
	listAll = sorted(dictQS.items(), key=lambda item: item[1][0], reverse=False) # sort by L_q_h asc
	numBins = math.ceil(lenDict / BINSIZE)
	#listAll = list(dictSort.items()) # (pidQ, pidS), (pidQLen * pidSLen, S', S_norm' = 0)
	
	lastBinMerged = False
	for i in range(numBins):
		if (i + 1) * BINSIZE <= lenDict: # except the last bin			
			# if the last bin is smaller than 0.5 * BINSIZE, then merge it with the preceding one
			if i == numBins - 2 and lenDict - (numBins - 1) * BINSIZE < 0.5 * BINSIZE:
				listBin = listAll[i * BINSIZE :]
				lastBinMerged = True
			else:
				listBin = listAll[i * BINSIZE : (i + 1) * BINSIZE]
				
		elif numBins == 1 or not lastBinMerged: # the last bin (which was not merged with the previous one)
			listBin = listAll[i * BINSIZE : ]
		else: 
			continue
		numTop = int(0.05 * len(listBin))
		
		# correction for small bins
		if numTop < 0.01 * BINSIZE: 
			numTop = len(listBin)

		listBin.sort(key=lambda x: x[1][1], reverse = True) # sort by S' desc
		listTop = listBin[:numTop] # select the top 5% S'(i.e. B_q_h) 
		
		# compute a, b using lin. regression (least squares) in log10(B_q_h) = a x log10(L_q_h) + b
		listX = []
		listY = []
		for j in range(len(listTop)):
			listX.append(listTop[j][1][0]) # L_q_h
			listY.append(listTop[j][1][1]) # B_q_h
		x = np.log10(np.array(listX))
		y = np.log10(np.array(listY))
		A = np.vstack([x, np.ones(len(x))]).T # transposed matrix --> 1st col: x values, 2nd col: 1s
		a, b = np.linalg.lstsq(A, y, rcond=None)[0]
		
		#print("a: " + str(a) + " b: " + str(b) + " q: " + q + " s: " + s + "\n")

		# update listBin: normalise each B_q_h --> B_q_h' = B_q_h / (10^b x L_q_h^a)
		try:
			for j in range(len(listBin)):
				#listBin[j][1][2] = float(listBin[j][1][1]) / ( pow(10, b) * pow(listBin[j][1][0], a) )
				listBin[j][1][2] = float(listBin[j][1][1]) / pow(10, b) / pow(listBin[j][1][0], a)
		except Exception as e:
			print("a: " + str(a) + " b: " + str(b) + " q: " + q + " s: " + s + "\n")
			print("The error is: ", e)
			sys.exit(2)

		# again sort by L_q_h asc
		listBin.sort(key=lambda x: x[1][0], reverse = False) # sort by L_q_h asc

		listAll[i * BINSIZE : i * BINSIZE + len(listBin)] = listBin
		listBin.clear()
	# end for
	
	# print the output --> (queryId, subjId, S'_norm)
	text = ""
	for i in range(len(listAll)):
		text += listAll[i][0][0] + "\t" + listAll[i][0][1] + "\t" + str(listAll[i][1][2]) + "\n"
		if i % 1000 == 0: # print in batches
			outFile.write(text)
			outFile.flush()
			text = ""
	# end for
	outFile.write(text) # print the last batch

# end printNormBitScore

if __name__ == "__main__":
	# ---------- gene length correction / the adapted OrthoFinder's alg. (Emms and Kelly, 2015 )-----------
	# for each query_taxID # q
	#	for each subj_taxID # h
	# 		read all rows for q, h (since the data are sorted by query_taxID, subj_taxID)
	#		while reading compute for each row: L_q_h (L_q x L_h)
	#		count rows for q, h
	#		sort by L_q_h
	#		split rows in bins of size = 10^3 (or 10^4??) according to L_q_h
	#		for each bin
	#			select the top 5% S'(i.e. B_q_h) 
	#			compute a, b using lin. regression (least squares) in log10(B_q_h) = a x log10(L_q_h) + b
	#			normalise each B_q_h --> B_q_h' = B_q_h / (10^b x L_q_h^a)

	INPUT_FILE_NAME, OUTPUT_FILE_NAME = parseArgv(sys.argv[1:])
	
	# in input file: rows are sorted by query_taxID (numerically), subject_taxID (n), query_pid (n), subject_pid (n)
	inFile = open(INPUT_FILE_NAME, "r")
	outFile = open(OUTPUT_FILE_NAME, "w")
	prevQ = ""
	prevS = ""
	dictQS = {}
	cntSelfHits = 0
	try: 
		for line in inFile: 
			# e.g. pid|0000000000000139259|tx|1392|	pid|0000000000398314972|tx|3983|	1.175E-38	149	180	180	187
			# query,target,evalue,bits,alnlen,qlen,tlen
			arr = line.strip().split("|")
			q = arr[3] # e.g. 1392
			s = arr[7]

			if q != prevQ or s != prevS: # new (query, subj) pair
				# print the results for the previous query
				if prevQ != "": # not the first line
					printNormBitScore(dictQS, outFile, q, s)
					dictQS.clear()
			pidQ = "pid|" + arr[1] + "|tx|" + q
			pidS = "pid|" + arr[5] + "|tx|" + s
			
			# skip self-hits
			if pidQ != pidS:
				pidQLen = arr[8].strip().split("\t")[3]
				pidSLen = arr[8].strip().split("\t")[4]
				dictQS[(pidQ, pidS)] = [float(pidQLen) * float(pidSLen), float(arr[8].strip().split("\t")[1]), 0.0] # dict of lists; better than dict of tuples
			else:
				cntSelfHits += 1
			
			prevQ = q
			prevS = s
		# end for		
		
		# print the last (query, subj) pair
		printNormBitScore(dictQS, outFile, q, s)
		
		inFile.close()
		outFile.close()
		print("Self-hits: " + str(cntSelfHits))
	except Exception as e:
		print("q: " + q + " s: " + s + "\n")
		print("The error is: ", e)
		sys.exit(2)
# end __main__














