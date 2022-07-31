"""
The cluster_function_enrichment_analysis_v4.py module performs enrichment analysis of functions within a given cluster by
performing a hypergeometric test.The standard (quant, sample, hit and total) values used in a hypergeometric test for
function enrichment in a cluster were:
•	quant		- number of annotations containing the tested function in cluster
•	sample		- number of annotations in cluster
•	hit		    - total number of annotations containing the tested function
•	total		- total number of annotations

The paths and options can be set on the start of the script by changing variables, see comments for instructions.

    --> Run (GO): time python3 cluster_function_enrichment_analysis_v4.py -s H -G -c 0.8 -S > clu_funct_analysis_HS_2022_c_0_8_go.out (running time ~ 15 minutes)
    --> Run (COG): time python3 cluster_function_enrichment_analysis_v4.py -s H -C -c 0.8 > clu_funct_analysis_HS_2022_c_0_8_cog.out (running time ~ 1 minute)

"""
__version__ = "1.0"
__author__ = "Tin Siroki"
__author__ = "Mirjana Domazet-Loso" # newer versions
# MDL changes made 2021-2022

import getopt # parse arguments
import sys # parse arguments

from os import walk
from scipy.stats import hypergeom
import statsmodels.stats.multitest as multi
import math


# path to folder containing files with functional annotations of genes (the files are generated with eggnog-mapper tool), and other auxiliary paths are hard-coded in __main__
#CATEGORIES_IN_PATH/INPUT_GO_PATH, CLUSTER_IN_PATH, HYPER_OUT_PATH, ANNOTATION_OUT
# File containing GO-basic functions
GO_OBO_PATH = "PATH_TO_THE_ANALYSIS_FOLDER/Data/go-basic.obo"

# Other options are set as the programs's arguments 

INPUT_SPECIES = {"H" : "Hsap_2022"}

# MDL, added March 15, 2021
def parseArgv(argv):
    print("Program arguments:")

    options = "hs:GCc:S" # options that require arguments should be followed by :
    long_options = ["help", "species", "GO_function", "COG_function", "c_value", "SLIM"]
    try:
        opts, args = getopt.getopt(argv, options, long_options)
    except getopt.GetoptError:
        print("Error. Usage: cluster_function_enrichment_analysis_mdl_v4.py -s <species> -G[-C] -c <c_value> [-S]")
        sys.exit(2)
    
    SPECIES = "" 
    GO = True
    C_VALUE = ""
    SLIM = False
    
    for opt, arg in opts:
        print(opt + "\t" + arg)
        if opt in ("-h", "--help"):
            print("Usage: cluster_function_enrichment_analysis_mdl_v4.py -s <species> -G[|-C] -c <c_value> [-S]")
            sys.exit()
        elif opt in ("-s", "--species"):
            if (arg in "ADHS"):		
                SPECIES = arg
            else:		
                print("Available species: H(omo sapiens)")
                sys.exit()
        elif opt in ("-G", "--GO_function"):
            GO = True
        elif opt in ("-C", "--COG_function"):
            GO = False
        elif opt in ("-c", "--c_value"): # 0.0, 0.4, 0.8
            C_VALUE = arg
        elif opt in ("-S", "--SLIM"):
            SLIM = True
    # end for
	
    # if COG, SLIM cannot be set
    if not GO: 
        SLIM = False
	
    return SPECIES, GO, C_VALUE, SLIM
# end parseArgv

if __name__ == "__main__":

    # MDL added main arguments (March 15, 2021)
    SPECIES, GO, C_VALUE, SLIM = parseArgv(sys.argv[1:])
    STR_PARAMETER_C = "c_0_" + str(C_VALUE)[2].rstrip()

    CATEGORIES_IN_PATH = "PATH_TO_THE_ANALYSIS_FOLDER/Data/eggnogResults/"
    
    CLUSTER_IN_PATH = "PATH_TO_THE_ANALYSIS_FOLDER/Clusters_2022/" + STR_PARAMETER_C + "/" + INPUT_SPECIES[SPECIES] + "_clusters.txt" # "Hsap_2022_clusters.txt"

    if GO:
        if SLIM: # go-basic
            dir_slim = "go-basic/"
            INPUT_GO_PATH = "PATH_TO_THE_ANALYSIS_FOLDER/Data/annotations_slim_gen_go-basic.gaf"
        else: #
            dir_slim = ""
            INPUT_GO_PATH = "PATH_TO_THE_ANALYSIS_FOLDER/Data/all_go_terms.gaf"
            
        HYPER_OUT_PATH = "PATH_TO_THE_ANALYSIS_FOLDER/results_2022/cluster_hyper_go/" \
                            + dir_slim + STR_PARAMETER_C + "/" + INPUT_SPECIES[SPECIES] + "_go.txt" # "/HomoSapiens_2020_go.txt"
        ANNOTATION_OUT = "PATH_TO_THE_ANALYSIS_FOLDER/results_2022/cluster_hyper_go/" \
                            + dir_slim + STR_PARAMETER_C + "/" + INPUT_SPECIES[SPECIES] + "_go_annotations.txt" # "/HomoSapiens_2020_go_annotations.txt"
        print(INPUT_GO_PATH + "\n" + HYPER_OUT_PATH)
					
    else: # COG
        SLIM = False
        HYPER_OUT_PATH = "PATH_TO_THE_ANALYSIS_FOLDER/results_2022/cluster_hyper_cog/" \
                            + STR_PARAMETER_C + "/" + INPUT_SPECIES[SPECIES] + "_cog.txt" # "/HomoSapiens_2020_cog.txt"
        ANNOTATION_OUT = "PATH_TO_THE_ANALYSIS_FOLDER/results_2022/cluster_hyper_cog/" \
                        + STR_PARAMETER_C + "/" + INPUT_SPECIES[SPECIES] + "_cog_annotations.txt" # "/HomoSapiens_2020_cog_annotations.txt"
	############
	
    gene_2_cat = {}
    file_paths = []
    go_2_slim = {}
    alt_2_orig = {}

    if GO:
        # inputs GO mappings to slim set
        with open(INPUT_GO_PATH) as f:
            for line in f:
                if line.startswith("!"):
                    continue

                line_split = line.split()
                go_id = line_split[0].strip("dummy_")
                go_id = "GO:" + go_id

                if go_id not in go_2_slim:
                    go_2_slim[go_id] = []

                go_2_slim[go_id].append(line_split[3])

        # inputs alternative function ids
        with open(GO_OBO_PATH, "r") as f:
            orig_id = ""
            for line in f:
                if line.startswith("id:"):
                    orig_id = line.split("id: ")[1].strip()
                if line.startswith("alt_id"):
                    alt_id = line.split("alt_id: ")[1].strip()
                    alt_2_orig[alt_id] = orig_id
    print("Reading files from folder: " + CATEGORIES_IN_PATH) # MDL test March, 2021
    for (dirpath, dirnames, filenames) in walk(CATEGORIES_IN_PATH):
         for f in filenames:
            if f.endswith(".annotations"):
                #file_paths.append(dirpath + "\\" + f)
                file_paths.append(dirpath + "/" + f) # MDL changes made March, 2021
                print("Including file:" + dirpath + "/" + f) # MDL test March, 2021

    file_counter = 0
    not_mapped_terms = set()
    alter_count = 0

    for file_path in file_paths:

        with open(file_path, "r") as f:
            print("inputing categories from file " + str(file_counter))
            file_counter += 1
            for line in f:
                if line.startswith("pgi"):
                    line_split = line.split("\t")

                    if GO:
                        cat_str = line_split[5].strip()
                    else:
                        cat_str = line_split[11].strip()

                    if cat_str == "" or (not GO and (cat_str == "R" or cat_str == "S")):
                        continue

                    cats = cat_str.split(",")
                    id_elements = line_split[0].split("|")
                    gene_id = id_elements[1] + "#" + id_elements[3]

                    if gene_id not in gene_2_cat:
                        gene_2_cat[gene_id] = set()

                    for cat in cats:
                        cat = cat.strip()
                        cat = cat.strip("\t")
                        if len(cat) == 0:
                            continue

                        if not GO and (cat == "R" or cat == "S"):
                            continue

                        if SLIM:
                            if cat in alt_2_orig:
                                cat = alt_2_orig[cat]
                                alter_count += 1

                            if cat not in go_2_slim:
                                not_mapped_terms.add(cat)
                                continue

                            slim_cats = go_2_slim[cat]

                            for sl_cat in slim_cats:
                                gene_2_cat[gene_id].add(sl_cat)

                        else:
                            if cat in alt_2_orig:
                                cat = alt_2_orig[cat]
                                alter_count += 1
                            if GO:
                                cat = cat.strip("GO:")
                            gene_2_cat[gene_id].add(cat)
            #end for
            print("\tFinished inputing categories from file " + str(file_counter - 1), flush = True)

    print("Number of annotations not mapped to slim GO: " + str(len(not_mapped_terms)), flush = True)
    print("Number of alter: " + str(alter_count), flush = True)
    print(not_mapped_terms, flush = True)

    for gene_id in gene_2_cat:
        gene_2_cat[gene_id] = list(gene_2_cat[gene_id])

    print("Inputed categories for " + str(len(gene_2_cat)) + " genes.", flush = True)

    clus_2_count_sample = {}
    cluster_name = ""
    cluster_sample_count = []
    sample_count = 0
    cluster_counter = 0
    cat_2_tot_count = {}
    total_count = 0

    first = True
    clus_2_res = {}

    all_cats = set()
    all_genes = set()

    no_annot_gene_count = 0
	# MDL here starts where changes need to be made!!!!!!!!!!!!!!!!! counting genes/annotations per cluster instead of counting annotations per genes!!
    cat_2_clus_count = {} # for each category count the number of clusters in which it has been assigned to at least one gene in the cluster
    cat_2_clus = {} # for each category of the current cluster set "1"
    clus_2_num_functions = {} # the number of different functions assigned to a cluster
    total_clu_count = 0
    ####
    
    with open(ANNOTATION_OUT, "w") as ann_out:
        with open(CLUSTER_IN_PATH) as f:
            for line in f:
                if line.strip() == "":
                    continue

                # new cluster(sample)
                if line.startswith("("):
                    if cluster_counter % 1000 == 0:
                        print("Analyzing cluster: " + str(cluster_counter), flush = True)
                    cluster_counter += 1
                    if SLIM or not GO:
                        ann_out.write("\n")
                        ann_out.write(line.split()[1] + "\n")
                    # end if SLIM or not GO
                    
                    if not first:
                        clus_2_res[cluster_name] = clus_2_count_sample.copy()
                        cluster_sample_count.append(sample_count)

                        # MDL changes made March 2021
                        clus_2_num_functions[cluster_name] = len(cat_2_clus)
                        total_clu_count += clus_2_num_functions[cluster_name]
                        for c in cat_2_clus:
                            if c not in cat_2_clus_count:
                                cat_2_clus_count[c] = 0
                            cat_2_clus_count[c] += 1
                        ###
                        
                    cluster_name = line.split()[1]
                    clus_2_count_sample = {}
                    first = False
                    sample_count = 0
                    
                    # MDL changes made March 2021
                    cat_2_clus.clear()
                    ###
                else:
                    g_id_split = line.split("|")
                    g_id = g_id_split[1] + "#" + g_id_split[3]
                    if g_id in gene_2_cat:
                        line_cats = gene_2_cat[g_id]
                        if g_id in all_genes:
                            print("Duplicate gene in input cluster file!!!")
                            print(g_id)
                        all_genes.add(g_id)

                        if len(line_cats) == 0:
                            print("No category for gene:")
                            print(g_id)

                        for c in line_cats:
                            if c not in clus_2_count_sample:
                                clus_2_count_sample[c] = 0
                            if c not in cat_2_tot_count:
                                cat_2_tot_count[c] = 0
                            
                            # MDL added March 23, 2021
                            if c not in cat_2_clus: cat_2_clus[c] = 1
                            ###
                            
                            clus_2_count_sample[c] += 1
                            sample_count += 1
                            if SLIM or not GO:
                                ann_out.write(g_id + "###" + c + "\n")
                            # MDL, added March 25, 2021    
                            #else: # if go-basic: for each gene, print list of ordinal numbers corresponding to go-basic functions
                                                       
                            all_cats.add(c)

                            cat_2_tot_count[c] += 1
                            total_count += 1
            '''
            storing last cluster
            '''
            clus_2_res[cluster_name] = clus_2_count_sample.copy()
            cluster_sample_count.append(sample_count)
            
            # MDL changes made March 2021
            clus_2_num_functions[cluster_name] = len(cat_2_clus)
            total_clu_count += clus_2_num_functions[cluster_name]
            for c in cat_2_clus:
                if c not in cat_2_clus_count:
                    cat_2_clus_count[c] = 0
                cat_2_clus_count[c] += 1
            ###

    print("Total number of genes with one or more categories: " + str(len(all_genes)), flush = True)
    print("Total number of different categories: " + str(len(all_cats)), flush = True)

    p_values = []
    total_sample_count = sum(cluster_sample_count)

    if total_sample_count != total_count:
        raise Exception("Total count and total sample count must be equal!!!")
    else:
        print("OK")

    # MDL delete lists/dicts which are not needed any more in order to release some memory (march 17, 2021)
    clus_2_count_sample.clear()
    gene_2_cat.clear()
    alt_2_orig.clear()
    go_2_slim.clear()

    cluster_counter = 0
    print(cluster_sample_count[:10], flush = True)
    if not p_values:
        p_adjusted = "NA"
    else:    
        p_adjusted = multi.multipletests(p_values, method="fdr_bh")[1]
    p_counter = 0

    print(cluster_sample_count[:10], flush = True)
    cluster_counter = 0
    with open(HYPER_OUT_PATH, "w") as out:
        # start mdl
        for c in clus_2_res:
            if cluster_counter % 1000 == 0:
            #if cluster_counter % 100 == 0:
                print("Storing cluster " + str(cluster_counter), flush = True)
            out.write(c + "\n")

            sample_count = cluster_sample_count[cluster_counter]
            temp_out = []

            for cat in clus_2_res[c]:
                #quant = clus_2_res[c][cat]
                #sample = sample_count
                #hit = cat_2_tot_count[cat]
                #total = total_count
                # MDL, changes made March 23, 2021
                quant = 1 # a cluster can be assigned a function at most once
                sample = clus_2_num_functions[c] # the number of functions assigned to the cluster
                hit = cat_2_clus_count[cat] # the number of clusters with the category cat
                total = total_clu_count
                
                ####
                
                # MDL, changes made March 22, 2021
                if GO and "GO" not in cat:
                    cat_print = "GO:" + str(cat)
                else:
                    cat_print = cat
                    
                out_line = cat_print + "\t" + str(quant) + "\t" + str(sample) + "\t" \
                           + str(hit) + "\t" + str(total) #+ "\t" + \
                           #str(p_values[p_counter]) + "\t" + str(p_adjusted[p_counter])

                if quant != 0 and (sample - quant) != 0 and (total - hit - sample + quant) != 0 and (hit - quant) != 0:
                    odds_sample = quant / (sample - quant)
                    odds_rest = (hit - quant) / (total - hit - sample + quant)
                    real_log_odds = math.log2(odds_sample / odds_rest)
                    out_line += "\t" + "%.2f" % real_log_odds
                    temp_out.append((real_log_odds, out_line))
                else:
                    out_line += "\tnan"
                    temp_out.append(("nan", out_line))
                #p_counter += 1

            #temp_out.sort()
            cluster_counter += 1

            for t in temp_out:
                out.write(t[1] + "\n")

            out.write("\n")

    print("\nFinished.", flush = True)
            