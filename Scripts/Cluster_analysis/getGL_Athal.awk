# Count gene gain/lost for a focal species (given taxID, e.g. 559292 S. cer.)
# Author: MDL, 21 Feb, 2017; MDL modified Feb 15, 2019
# Rules:
	# if (foundT) { # cluster contains genes of the focal species
		# (1a) if cluster contains only members (genes) of the focal species, then I'll suppose that they originated at its minPS (== maxPS), and are not lost (clusterRepresentative is from the focal species)
		# (2) if cluster also contains other species' genes, then the gene family was gained in minPS and still exists (not lost) (clusterRepresentative doesn't have to be from the focal species)
	
	# (1b) if cluster contains only members (genes) of one (non-focal) species, 
	# then I'll suppose that they originated at the species PS (or its lineage) --> CAN BE IGNORED FROM THE POINT OF FOCAL SPECIES
	
	# (3) multiple species, and belonging to different PS; from t's point of view: gene family gained at minPS, gene family lost at (maxPS+1)

	# (4) multiple species, and belonging to the same PS --> CAN BE IGNORED FROM THE POINT OF FOCAL SPECIES
		# (4a) if LCA of these species equals the cluster's minPS (== maxPS), than from t's point of view: gene family gained at minPS, gene family lost at (maxPS+1)
		# (4b) if LCA of these species does not equal the cluster's minPS (== maxPS), than UNRESOLVED (i.e. for example their LCA with the focal species is "Cell. organisms", but their LCA is "Viridiplantae")
	
	# Only cases (foundT) (i.e. 1a and 2) and (3) are USED FOR COMPUTING GENE FAMILIES GAINED/LOST FROM THE FOCAL SPECIES PERSPECTIVE

# Example: time awk -v t=3702 -v P=24 -f getGL_Athal.awk PATH_TO_ANALYSIS_FOLDER/results_0_8/dbAllPlus_clu_0_8_Athal_out_EXT.txt > PATH_TO_ANALYSIS_FOLDER/results_0_8/res_dbAllPlus_clu_0_8_Athal_out_EXT.txt
#

# fill a 
function get_a(a ) {
	# all A. thal.'s parents up to the root were precomputed using treeTaxId.awk
	a[1 ] = 3702	# Arabidopsis_thaliana
	a[2 ] = 3701	# Arabidopsis
	a[3 ] = 980083	# Camelineae
	a[4 ] = 3700	# Brassicaceae
	a[5 ] = 3699	# Brassicales
	a[6 ] = 91836	# Malvids
	a[7 ] = 7000006	# Eurosids
	a[8 ] = 7000005	# Superrosids
	a[9 ] = 7000003	# Core_eudicots
	a[10] = 71240	# Eudicots
	a[11] = 3398	# Magnoliophyta
	a[12] = 58024	# Spermatophyta
	a[13] = 7000001	# Embriophyta
	a[14] = 35493	# Streptophyta
	a[15] = 33090	# Viridiplantae
	a[16] = 3000009	# Archaeplastida
	a[17] = 3000007	# Diaphoretickes
	a[18] = 3000005	# Excavata/Diaphoretickes
	a[19] = 2759	# Eukaryota
	a[20] = 3000004	# Asgard_archaea/Eukaryota
	a[21] = 3000003	# TACK/Asgard_archaea/Eukaryota
	a[22] = 3000002	# Euryarchaeota/TACK/Asgard_archaea/Eukaryota
	a[23] = 3000001	# DPANN/Euryarchaeota/TACK/Asgard_archaea/Eukaryota
	a[24] = 131567	# Cellular_organisms
}

# compute taxID from an id of the form "pgi|0000000003120175935|ti|312017|pi|0|"
function getTaxId(id) {
	#Example: id == "pgi|0000000003120175935|ti|312017|pi|0|" ---> taxid=312017
	pos = index(id, "ti|");
	s = substr(id, pos + 3); # skip "ti|"
	# find next "|", which ends taxid
	pos2 = index(s, "|");
	taxid = substr(s, 1, pos2 - 1);
    return taxid;
} # end function

# compute PS of the previous taxid and update geneListG/geneListL for the focal species
function getClusterGainLost(foundT, maxPS, minPS, arr_gained, arr_lost, geneListG, geneListL, clusterRepresentative, cntCluMembers, taxid, n, PS, arr_unresolved, cnt_arr_unresolved, cnt) {
	
	++cntGetClusterGainLost; # DEBUG purpose - how many times the function was called
	if (foundT) { # cluster contains genes of the focal species
		# (1a) if cluster contains only members (genes) of the focal species, then I'll suppose that they originated at its minPS (== maxPS), and are not lost (clusterRepresentative is from the focal species)
		# (2) if cluster also contains other species' genes, then the gene family was gained in minPS and still exists (not lost) (clusterRepresentative doesn't have to be from the focal species)
		
		++ arr_gained[minPS];
		# GENE FAMILY represented by clusterRepresentative
		geneListG[clusterRepresentative] = minPS;
		geneListL[clusterRepresentative] = -1; # flag - meaning not lost
		
		#DEBUG - it works now, but only one example - try with more than one example!!!
		print "found B. sub. (cases 1a or 2)\t" clusterRepresentative " gained = " minPS
		++ cntFoundT; # DEBUG purpose - how many times foundT
	} # end if
	
	else { # if the cluster doesn't contain the focal species
		# (1b) if cluster contains only members (genes) of one (non-focal) species, 
		# then I'll suppose that they originated at the species PS (or its lineage), and CAN BE IGNORED FROM THE POINT OF FOCAL SPECIES!!!
		if (cntCluMembers[taxid] == cnt) {
			# nothing done with this gene family
			#DEBUG
			print "non-focal species (case 1b)\t" clusterRepresentative "cntMembers = " cnt
			++ cntCase1b; # DEBUG purpose
		}
		# (3) multiple species, and belonging to different PS; from t's point of view: gene family gained at minPS, gene family lost at (maxPS+1)
		else if (minPS != maxPS) { 
			++ arr_gained[minPS];
			++ arr_lost[maxPS + 1];
			geneListG[clusterRepresentative] = minPS; # Is procedure as in (4) needed here? No, since that would mean the multiple origin of the gene family
			geneListL[clusterRepresentative] = maxPS + 1;
			
			++ cntCase3; # DEBUG purpose
			print "Case 3\t" clusterRepresentative " gained = " minPS " lost = "  maxPS + 1;
		}
		# (4) multiple species, and belonging to the same PS; 
		# (4a) if LCA of these species equals the cluster's minPS (== maxPS), than from t's point of view: gene family gained at minPS, gene family lost at (maxPS+1)
		# (4b) if LCA of these species does not equal the cluster's minPS (== maxPS), than UNRESOLVED
		else {
			# find LCA of the species in the cluster !! (it could be more than one species)
			# find all taxa in this cluster
			cntTaxaCluster = 0;
			for (i = 1; i <= n; i++) { 
				if (cntCluMembers[arr_taxID[i]]) {
					taxaInCluster[cntTaxaCluster ++ ] = arr_taxID[i];
				}
			} # end for

			# find lca of the taxa in the cluster
			foundExpectedLCA = 0;
			for (i = 0; i < cntTaxaCluster && !foundExpectedLCA; i++) {
				for (j = i + 1; j < cntTaxaCluster && !foundExpectedLCA; j++) { 
					lowest_LCA = LCA[taxaInCluster[i], taxaInCluster[j]];
					for (k = 0; k < PS; k++) {
						if (a[k] == lowest_LCA && PS - k == minPS) { # Can it be that both conditions are not satisfied?					
							foundExpectedLCA = 1;
							break; 
						}
					}
				}
			} # end for

			# (4a) if LCA of these species equals the cluster's minPS (== maxPS), than from t's point of view: gene family gained at minPS, gene family lost at (maxPS+1)
			if (foundExpectedLCA) {
				#++ arr_gained[minPS];
				#++ arr_lost[maxPS + 1];
				#geneListG[clusterRepresentative] = minPS;
				#geneListL[clusterRepresentative] = maxPS + 1;
				++ cntCase4a; # DEBUG purpose
				print "Case 4a\t" clusterRepresentative " gained = " minPS " lost = "  maxPS + 1; 
			}
			# (4b) if LCA of these species does not equal the cluster's minPS (== maxPS), than UNRESOLVED
			else {
				arr_unresolved[cnt_arr_unresolved] = clusterRepresentative;
				++ cnt_arr_unresolved;		
				++ cntCase4b; # DEBUG purpose		
				print "Case 4b\t" clusterRepresentative " unresolved"; 				
			}
		} # end else
	}
} # end function

BEGIN {  
	if (ARGC != 2){
		print "usage: awk -v t=<taxId> -v P=<taxPS> -f getGL.awk <targetDB_cluster_file>";
		exit;
	}
	
	PS = P; # for S.cer.
	get_a(a ); # compute array a containing the focal species parental taxa stored in the corresponding ps, i.e. a[i] is the taxonID of the taxon at the phylostratum = PS-i (a[0] = 7227, PS of 7227 is 48)
	# DEBUG - OK
	# for (i = 0; i < PS; i++) print i "\t" a[i];
		
	MIN_PS = 1000; # minimal PS
	MAX_PS = 0; # maximal PS
	minPS = MIN_PS;
	maxPS = MAX_PS;			
	
	foundT = 0; # flag: cluster contains focal species
	cnt_arr_unresolved = 0;
	cntCluster = 0;
	
	# flags used to read the LCA matrix at the beginning of the file
	firstRow = 1;
	secondRow = 0;
	thirdRow = 0;
	readLCA = 0;	
}
{
	# ----------------------- READING LCA MATRIX -----------------------
	# the first row is n (number of taxa), followed by two rows containing n taxa's id-s and names
	# then follows the matrix with n x n values LCA between each pair of taxa (the two leftmost columns in the matrix are again taxa's id-s and names)
	if (firstRow) {
		n = $1;
		firstRow = 0;
		secondRow = 1;
	}
	else if (secondRow) {	# taxa's id-s
		for (i = 1; i <= NF; i++) arr_taxID[i] = $i;
		numTaxId = NF;
		secondRow = 0;
		thirdRow = 1;
	}
	else if (thirdRow) {	# taxa's names
		for (i = 1; i <= NF; i++) arr_taxNames[i] = $i;
		thirdRow = 0;
		readLCA = 1;
	}
	else if (readLCA) { # read LCA matrix (matrix data start from the third column)
	    #13037 (Danaus_ple)         -1    1708736    1708736    1708696 itd.
		# LCA[ taxID1, taxID2]
		for (i = 1; i <= n; i++) {			
			LCA[$1, arr_taxID[i]] = $(i + 2);
		}
		if (NR == n + 3) {
			readLCA = 0; # the last row has been read
		}
	}
	
	# ----------------------- PROCESSING CLUSTERS -----------------------
	# Ignoring situations which cannot be reconstructed:
	# (i) that a gene was lost and gained somewhere in the past, but there is no trace of it (---> possible correction model should be added in the future)
	# (ii) Wagner's parsimony - that a gene can be gained in multiple branches (--> we use Dollo's parsimony, which does not include this model)

	# Possible cases:
	# (1) If the focal species (e.g. 7227) is included as a cluster member, 
	# 		then the corresponding gene was gained at the phylostratum with the smallest value in the cluster, and it was never lost in the focal species part of the tree	
	# (2) If the focal species is not included in a cluster
	# 		(2a) if the LCA(taxa whose ps is equal to maxPS) == maxPS, then the gene family can be considered to be lost for the focal species in ps = maxPS + 1??
	# 		(2a) if the LCA(taxa whose ps is equal to minPS) == minPS, then the gene family can be considered to be gained for the focal species in ps = minPS
	# (3) Unresolved cases
	# 		(3a) if the LCA(taxa whose ps is equal to minPS) != minPS, we cannot tell whether the gene family originated in minPS or in some other PS, which is LCA(taxa whose ps is equal to minPS),
	# 		but is not in the focal species' lineage ????????? UNRESOLVED CASE --> count these cases!!
	# 		(3b) if the LCA(taxa whose ps is equal to maxPS) != maxPS, we cannot tell whether the gene family was lost in maxPS or in some other PS, which is LCA(taxa whose ps is equal to maxPS),
	# 		but is not in the focal species' lineage ????????? UNRESOLVED CASE --> count these cases!!
	# (4) Ignored case:
	# 		Only one gene or multiple genes from the same taxon in a gene family --> I'll suppose that such a gene (or possibly gene duplicates) is (are) unique for the taxon (or the lineage)
	# (5) All other cases: ?
	# Possible gene family statuses: G (gained), L (lost), PG (possibly gained - unresolved case 3a), PL (possibly lost - unresolved case 3b), E (still existant in the focal speices), I (ignored), ? (all other cases)
		
	# Example of a cluster (i.e. gene family)
	#(1) Cluster representative: pgi|0000000000072226106|ti|7222|pi|0|	Number of members in the cluster: 3
	#	Cluster members: 
	#	(1) pgi|0000000000072226106|ti|7222|pi|0|
	#	Drosophila_grimshawi	42
	#	(2) pgi|0000000000693194578|ti|69319|pi|0|
	#	Microplitis_demolitor	26
	#	(3) pgi|0000000000034740847|ti|34740|pi|0|
	#	Heliconius_melpomene	28

	if ($0 == "") { # skip empty rows
	}
	
	# --------- NEW CLUSTER --------- 
	#(2) Cluster representative: pgi|0000000005412619063|ti|54126|pi|0|	Number of members in the cluster: 1
	#Cluster members: 
	#(1) pgi|0000000005412619063|ti|54126|pi|0|
	#Pristionchus_pacificus	15

	else if ($0 ~ /Cluster representative/) { # new cluster
		# compute PS of the previous taxid and update geneListG/geneListL for the focal species (except for the first cluster, when minPS == MIN_PS)
		if (minPS != MIN_PS) {
			-- cnt;
			getClusterGainLost(foundT, maxPS, minPS, arr_gained, arr_lost, geneListG, geneListL, clusterRepresentative, cntCluMembers, taxid, n, PS, arr_unresolved, cnt_arr_unresolved, cnt);
			}
		else { # do nothing - this is for the first cluster representative, i.e. nothing to be computed, since it has no predecessor
		}			

		# new cluster representative
		numMembers = $NF;
		clusterRepresentative = $4; # e.g. pgi|0000000002251641234|ti|225164|pi|0|
		# DEBUG - ok
		print "\n(" ++cntCluster ") clusterRepresentative = " clusterRepresentative "\tnumMembers = " numMembers;
		
	} # end cluster representative
	
	# --------- CLUSTER MEMBERS --------- 
	else { # cluster members
		if ($0 ~ /Cluster members/) { # skip row
			# prepare variables for this new cluster
			foundT = 0;
			cnt = 1; # count number of members in the cluster (>= 1)
			minPS = MIN_PS;
			maxPS = MAX_PS;
			secondRowCM = 0;
			# set count-s for all taxa to 0
			for (i = 1; i <= numTaxId; i++) cntCluMembers[arr_taxID[i]] = 0; 				
			# set find for the first member
			find = "(1)" # used for recognizing taxid row, e.g. (1) pgi|0000000005412619063|ti|54126|pi|0|	
		} 
		else {
			# index(in, find)			
			if (index($0, find) > 0) { # taxid row, e.g. (1) pgi|0000000005412619063|ti|54126|pi|0|				
				taxid = getTaxId($2)
				++ cntCluMembers[taxid]; # number of taxon's genes in the cluster
				if (taxid == t) { # determine whether the cluster contains focal taxon (e.g. 7227)
					foundT = 1;
				}
				secondRowCM = 1;
			}
			else if (secondRowCM) { # get PS from this row; e.g. Heliconius_melpomene	28
				secondRowCM = 0;
				ps = $2; # this is computed ok
				if (minPS > ps) { # new min PS
					minPS = ps;
					minTaxid = taxid;	
				}
				if (maxPS < ps) { # new max PS
					maxPS = ps;
					maxTaxid = taxid;	
				}
				# DEBUG
				print find "\t" taxid "\t" ps			
				# prepare for the next member
				++ cnt; 
				find = "(" cnt ")" # used for recognizing taxid row, e.g. (1) pgi|0000000005412619063|ti|54126|pi|0|
			}		
		} # end else		
	} # end else		
	
}
END {
	# DEBUG - test n and lca read correctly (168 - 192) - OK
	#print "n = " n;	
	#for (i = 1; i <= n; i++) {
	#	for (j = 1; j <= n; j++) 
	#		printf ("%11d", LCA[arr_taxID[i], arr_taxID[j]]);
	#		printf("\n");
	#}
	
	# compute PS of the last taxid and update geneListG/geneListL for the focal species
	-- cnt;
	getClusterGainLost(foundT, maxPS, minPS, arr_gained, arr_lost, geneListG, geneListL, clusterRepresentative, cntCluMembers, taxid, n, PS, arr_unresolved, cnt_arr_unresolved, cnt);

	print "\n\n************ Focal Species Gene Families Gained/Lost per Phylostratum **************"
	print "PS\tGained\tLost\tTotal"; # involves cases foundT (1a, 2), case 3, and case 4a
	# unresolved genes not included
	for (i = 1; i <= PS; i++) {
		print i "\t" arr_gained[i] "\t" arr_lost[i];
	}
	
	print "\n\n************ List of Genes (Phylostratum Gained / Phylostratum Lost) **************"	
	print "Gene Family\t\t\t\t\t\t\t\tPS_Gained\t\tPS_Lost"
	for (representative in geneListG)  {
		# cluster representative --> gene family (geneList is geneFamilyList)
		print representative "\t" geneListG[representative] "\t\t\t\t" geneListL[representative];
	}
	
	print "\n\n************ List of Unresolved Genes **************"
	for (i = 0; i < cnt_arr_unresolved; i++) print arr_unresolved[i];
	
	# DEBUG - works fine
	#for (i = 1; i <= numTaxId; i++) print arr_taxID[i] "\t" arr_taxNames[i];

	print "cntGetClusterGainLost = " cntGetClusterGainLost;
	print "cntFoundT = " cntFoundT; # case1a and case2
	print "cntCase1b = " cntCase1b;
	print "cntCase3 = " cntCase3;
	print "cntCase4a = " cntCase4a;
	print "cntCase4b = " cntCase4b;
}





