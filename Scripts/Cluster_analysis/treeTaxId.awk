# Get the phylogeny tree for a tax id (print all parents up to the root)
# Author: Mirjana Domazet-Lošo, 2017-2022
#
# Usage: awk -v t=535026 -f treeTaxId.awk nodes.txt > Parents/535026parents.txt
#
BEGIN {  
	if (ARGC != 2){
		print "usage: awk -v t=<taxId> -f treeTaxId.awk <nodes-file> > <outputFile>";
		exit;
	}
	################
	# EDIT LIST - for each node in nodes.txt, add a row in the form: list[nodeID] = parent_nodeID, e.g.
	# list[9606] = 131567
	################
}
END {
	fileNames = "names.txt";
	while (1) {
		if (t in list) { # child: t
			# find childID parentID parent name
			
			cmd = "awk '{ if($1 == \"" t "\") { print $1 \"\t\" $2; exit; } }' " fileNames; 
			#print cmd;
			system(cmd);				
			#print t "\t" list[t]; #  print child parent
			t = list[t]; # new child		
		}
		if (t == 131567) { # stop when the root (Cellular_organisms) in encountered
			print "131567\tCellular_organisms"; #  print root
			break;
		}
	}
}
	
