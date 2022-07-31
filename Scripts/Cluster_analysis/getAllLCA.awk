# Compute all LCA for all pairs of taxa
# Author: MDL, 2017-2018
# Usage: awk -f getAllLCA.awk allParents.txt > allLCA.txt
#
BEGIN {	
	newTaxon = 1;
}
{
	# read all taxa
	# 	foreach new taxon read all its parents up to the root
	if (newTaxon) {
		taxon = $1; # taxId
		names[taxon] = $2;
		cnt[taxon] = 0;
		newTaxon = 0;
		i = 0;
	}
	else { # taxon's parents
		parent = $1; #taxId
		a[taxon, ++i] = parent;
		++ cnt[taxon];
		if (parent == 131567) { # Cellular_organisms
			newTaxon = 1; # get ready for the next taxon
		}
	}
}
END {
	for (t in names) { # for each taxon find its LCA with all other taxa t2 != t
		for (t2 in names) {
			if (t == t2) {
				LCA[t, t2] = -1;
			}
			else {
				found = 0; # flag for found LCA
				for (i = 1; i <= cnt[t]; i++) { # for all its parents				
					for (j = 1; j <= cnt[t2]; j++) {
						if (a[t, i] == a[t2, j]) {
							found = 1;
							break;
						}				
					}
					if (found) {
						LCA[t, t2] = a[t, i]; # or: a[t2, j]
						break;
					}
				} # end for	3
			}
		} # end for-2		
	} # end for-1
	
	# print all LCA
	printf("%23s", ""); # print header row and column names
	for (t in names) { # taxID
		printf("%11.10s", t);
	}
	printf("\n");
	printf("%23s", ""); # print header row and column names
	for (t in names) { # taxname
		printf("%11.10s", names[t]);
	}
	printf("\n");

	# print rows with row names
	for (t2 in names) {
		printf("%10.10s (%10.10s)", t2, names[t2]);
		for (t in names) {
			printf("%11.10s", LCA[t2, t]);
		}		
		printf("\n");
	}	
	printf("\n");
}

