"""
The cluster_function_enrichment_summary_go.py	 module performs processing of function enrichment per cluster for GO terms.
The module gives a function (one_fun) or functions (multi_fun) to a cluster.
A function is assigned to a cluster if at least one gene in the cluster is annotated with the corresponding function.
Alternatively, if topPercent of genes are assigned the function.
The paths and options can be set on the start of the script by changing variables, see comments for instructions.
"""

"""
MDL changes made March, 2021
    --> Run (GO): time python3 cluster_function_enrichment_summary_go_v4.py -s H -c 0.8 > clu_funct_summary_HS_2022_c_0_8_go.out 
    (running time ~ 4 minutes)
"""

__version__ = "1.0"
__author__ = "Tin Siroki"
__author__ = "Mirjana Domazet-Loso"
# changes made by MDL 2021-2022

import getopt # parse arguments
import sys # parse arguments

import operator
import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt
import numpy as np

# cut-off p_value used ass significant
P_ADJ_CUTOFF = 0.05

CATEGORY_INFO_PATH = "PATH_TO_THE_ANALYSIS_FOLDER/Data/go.obo"
PHYL_NAMES_FOLDER_PATH = "PATH_TO_THE_ANALYSIS_FOLDER/PS_names_short/"

# Include no annotation for cluster without significant function
INCLUDE_NO_ANNOTATION = True

INPUT_SPECIES = {"H" : "Hsap_2022"}

# MDL, added March 26, 2021
def parseArgv(argv):
    print("Program arguments:")

    options = "hs:c:Sp" # options that require arguments should be followed by :
    long_options = ["help", "species", "c_value", "SLIM", "plot"]
    try:
        opts, args = getopt.getopt(argv, options, long_options)
    except getopt.GetoptError:
        print("Error. Usage: cluster_function_enrichment_analysis_mdl_v2.py -s <species> -c <c_value> [-S] [-p]")
        sys.exit(2)
    
    SPECIES = "" 
    C_VALUE = ""
    ALL = True
    PLOT = False
    
    for opt, arg in opts:
        print(opt + "\t" + arg)
        if opt in ("-h", "--help"):
            print("Usage: cluster_function_enrichment_analysis_mdl_v2.py -s <species> -c <c_value> [-N] [-p]")
            sys.exit()
        elif opt in ("-s", "--species"):
            if (arg in "ABDHS"):		
                SPECIES = arg
            else:		
                print("Available species: H(omo sapiens)")
                sys.exit()
        elif opt in ("-c", "--c_value"): # e.g. 0.0, 0.4, 0.8
            C_VALUE = arg
        elif opt in ("-S", "--SLIM"):
            ALL = False
        elif opt in ("-p", "--plot"):
            PLOT = True # plot pdf
    # end for
	
    return SPECIES, C_VALUE, ALL, PLOT
# end parseArgv

# MDL added February, 2021
def get_cluster_top_function_v2(functions, gene_count, topPercent=0.3):
    """
    Returns the best significant function that is present in more than topPercent of the genes in the cluster
    function per cluster, empty list is returned if cluster does not contain significant functions.
    :param functions:
    :param gene_count:
    :param topPercent: return those functions which are assigned to at least topPercent genes (default=0.3) 
    : 					i.e. return those functions annotated to at least 30% of genes in a cluster
    :return:
    """
    #min_p_adj = float('+inf')
    top_function = ""
    for fun in functions:
        #if fun["p_adj"] < min_p_adj and fun["hit"] >= (gene_count/2):
        if fun["hit"] >= topPercent * gene_count:
            top_function = fun["info_all"]
            #min_p_adj = fun["p_adj"]

    res = []
    if top_function != "":
        res.append(top_function)
    return res

# MDL added February, 2021
def get_all_clusters_function_v2(functions, gene_count):
    """
    Returns all significant functions that are present in more than half of the genes in the cluster
    function per cluster, empty list is returned if cluster does not contain significant functions.
    :param functions:
    :param gene_count:
    :return:
    """
    fun_filtered = []
    for fun in functions:
        #if fun["p_adj"] < P_ADJ_CUTOFF and fun["hit"] >= (gene_count/2):
        fun_filtered.append(fun)

    #fun_filtered.sort(key=operator.itemgetter("p_adj"))
    res = []
    for el in fun_filtered:
        res.append(el["info_all"])
    return res



def make_gain_loss_count_summary(out_path, out_counts_phyl, out_10_perc, out_all_tab, out_all_annotations,
                               out_total_category_counts, pdf_path, pdf_path_one, path_phyl_2_name,
                               in_path, cluster_in_path, gain_function, plot, only_top=True):
    cluster_2_info = {}
    cluster_2_hyper = {}
    go_id_2_name = {}
    total_count = 0
    cat_2_total_count = {}
    phyl_2_name = {}
    phyl_counter = 1

    with open(path_phyl_2_name, "r") as f:
        first = True
        for line in f:
            if first:
                first = False
                continue
            line_split = line.split("\t")
            phyl_2_name[phyl_counter] = line_split[1].strip()
            phyl_counter += 1

    min_ps = 1
    max_ps = len(phyl_2_name)

    with open(CATEGORY_INFO_PATH) as f:
        id = ""
        for line in f:
            if line.startswith("id: "):
                id = line.strip("id: ").strip()
            if line.startswith("name: "):
                name = line[6:].strip()
                go_id_2_name[id] = name
    go_id_2_name["no_annotation"] = "Cluster without assigned function"

    with open(cluster_in_path, "r") as f:
        for line in f:
            if line == "" or line == "\r" or line == "\n":
                continue

            line = line.strip()
            line_split = line.split()

            if len(line_split) < 2:
                continue

            if line_split[1].startswith("pid"):
                cl_in = {}
                cl_in["id_line"] = "\t".join(line_split[1:]).strip()
                cl_in["ps_start"] = int(line_split[3])
                #cl_in["ps_end"] = int(line_split[4]) + 1
                if line_split[4] == "-": cl_in["ps_end"] = "-"
                else: cl_in["ps_end"] = int(line_split[4])	# MDL changed 28/01/2021: removed (+ 1)	
                cl_in["gene_count"] = int(line_split[5].strip())
                cluster_2_info[line_split[1]] = cl_in.copy()

    print("Number of clusters: " + str(len(cluster_2_info)))

    with open(in_path, "r") as f:

        cluster_id = ""
        for line in f:
            if line == "" or line == "\r" or line == "\n":
                continue
            line_split = line.split("\t")

            if line_split[0].startswith("pid"):
                cluster_id = line_split[0].strip()
                cluster_2_hyper[cluster_id] = []
            else:
                #if ALL and not "GO" in line_split[0]: continue
                fun = {}
                fun["id"] = line_split[0]
                #print(str(line_split[1]))
                fun["hit"] = int(line_split[1])				
				# MDL changes made for v3; March 22, 2021
                #fun["p_value"] = float(line_split[5])
                #fun["p_adj"] = float(line_split[6])
                #fun["log_ratio"] = float(line_split[7].strip())
                fun["log_ratio"] = float(line_split[5].strip())
                fun["info_all"] = line.strip()

                cluster_2_hyper[cluster_id].append(fun.copy())

    print("Number of clusters hyper: " + str(len(cluster_2_hyper)))

    annotated_clusters = 0
    no_annot_cluster_count = 0
    phyl_2_cat = {}
    all_fun = set()

    with open(out_path, "w") as out:
        with open(out_all_annotations, "w") as out_ann:
            first = True

            for cl_id in cluster_2_hyper:
                if only_top:
                    #top_functions = get_cluster_top_function(cluster_2_hyper[cl_id],
                    # MDL changed February, 2021
                    top_functions = get_cluster_top_function_v2(cluster_2_hyper[cl_id],					
                                                             cluster_2_info[cl_id]["gene_count"])
                else:
                    # top_functions = get_all_clusters_function(cluster_2_hyper[cl_id],
                    # MDL changed February, 2021
                    top_functions = get_all_clusters_function_v2(cluster_2_hyper[cl_id],
                                                        cluster_2_info[cl_id]["gene_count"])

                ps = None
                if gain_function == "gain":
                    ps = cluster_2_info[cl_id]["ps_start"]
                else:
                    ps = cluster_2_info[cl_id]["ps_end"]
                    #if ps > max_ps: # mdl changes January, 2021
                    if ps == "-":
                        continue

                if INCLUDE_NO_ANNOTATION:
                    if len(top_functions) == 0:
                        no_annot_cluster_count += 1
                        top_functions.append("no_annotation")

                if len(top_functions) > 0: # MDL changes should be made here!!!!!!!!!!!!
                    if not first:
                        out.write("\n")

                    out.write(cluster_2_info[cl_id]["id_line"] + "\n")
                    annotated_clusters += 1
                    first = False

                    for fun in top_functions:
                        function_id = fun.split("\t")[0]
                        all_fun.add(function_id)

                        if function_id not in cat_2_total_count:
                            cat_2_total_count[function_id] = 0

                        cat_2_total_count[function_id] += 1 # 
                        total_count += 1
                        out_ann.write(cl_id + "###" + function_id + "\t" + str(ps) + "\n")

                        if ps not in phyl_2_cat:
                            phyl_2_cat[ps] = {}

                        if function_id not in phyl_2_cat[ps]:
                            phyl_2_cat[ps][function_id] = 0

                        phyl_2_cat[ps][function_id] += 1
                        out.write(fun + "\n")

    print("Number of cluster without annotations: " + str(no_annot_cluster_count))

    for i in range(1, max_ps + 1):
        if i not in phyl_2_cat:
            phyl_2_cat[i] = {}

    tot_counts = []
    with open(out_total_category_counts, "w") as out:
        out.write("category\ttotal_count\n")
        for cat in cat_2_total_count:
            line = cat + "\t" + str(cat_2_total_count[cat])
            tot_counts.append(cat_2_total_count[cat])
            out.write(line + "\n")

    if sum(tot_counts) != total_count:
        raise Exception("Total counts error!!!")

    #print("Total annotation count: " + str(total_count))
    # MDL changed March, 2021
    print("Total annotation count: " + str(total_count))

    with open(out_counts_phyl, "w") as out:
        with open(out_10_perc, "w") as out_10:
            phyl_keys = list(phyl_2_cat.keys())
            phyl_keys.sort()

            first = True
            for ps in phyl_keys:
                if not first:
                    out.write("\n\n\n")
                    out_10.write("\n\n\n")
                first = False

                counts_sorted = sorted(phyl_2_cat[ps].items(), key=lambda kv: kv[1], reverse=True)
                total = sum(x[1] for x in counts_sorted)

                out.write("Phylostratum " + str(ps) + "\tcat_count=" + str(len(counts_sorted)) + "\ttotal=" + str(
                    total) + ":\n")

                for el in counts_sorted:
                    # MDL changes made March 22, 2021
                    #out.write(el[0] + "\t" + go_id_2_name[el[0]] + "\t" + str(el[1]) + "\n")
                    if el[0] in go_id_2_name:
                    	out.write(el[0] + "\t" + go_id_2_name[el[0]] + "\t" + str(el[1]) + "\n")
                    else:
                    	out.write(el[0] + "\t" + "UNKNOWN_FUNCTION" + "\t" + str(el[1]) + "\n")
						
                top_len = int(round(len(counts_sorted) * 0.1))

                out_10.write("Phylostratum " + str(ps) + "\tcat_count=" + str(top_len) + "\ttotal=" + str(
                    total) + ":\n")

                for el in counts_sorted[:top_len]:
                    #out_10.write(el[0] + "\t" + go_id_2_name[el[0]] + "\t" + str(el[1]) + "\n")
                    # MDL changes made March 22, 2021
                    if el[0] in go_id_2_name:
                    	out_10.write(el[0] + "\t" + go_id_2_name[el[0]] + "\t" + str(el[1]) + "\n")
                    else:
                    	out_10.write(el[0] + "\t" + "UNKNOWN_FUNCTION" + "\t" + str(el[1]) + "\n")

    phyl_keys = list(phyl_2_cat.keys())
    phyl_keys.sort()
    all_fun = list(all_fun)
    cat_ps = [[0 for i in range(len(phyl_keys))] for j in range(len(all_fun))]

    with open(out_all_tab, "w") as out:

        head = "category_id\tname\t"

        for ps in phyl_keys:
            head += str(ps) + "\t"

        head = head.strip() + "\n"
        out.write(head)

        for ps in phyl_keys:
            j = int(ps) - 1

            for cat in phyl_2_cat[ps]:
                i = all_fun.index(cat)
                cat_ps[i][j] = phyl_2_cat[ps][cat]

        for i in range(len(all_fun)):
            # MDL changes made March 22, 2021
            #line = all_fun[i] + "\t" + go_id_2_name[all_fun[i]] + "\t"
            if all_fun[i] in go_id_2_name:
                line = all_fun[i] + "\t" + go_id_2_name[all_fun[i]] + "\t"
            else:
                line = all_fun[i] + "\t" + "UNKNOWN FUNCTION" + "\t"
				
            for j in range(len(phyl_keys)):
                line += str(cat_ps[i][j]) + "\t"
            line = line.strip()

            out.write(line + "\n")

    if not PLOT:
        return

    pdf = matplotlib.backends.backend_pdf.PdfPages(pdf_path_one)
    plt.style.use('ggplot')
    x_axis = range(0, len(phyl_keys))
    plt.figure()

    avg_phyl = np.zeros(len(phyl_keys))

    for i in range(len(all_fun)):
        for j in range(len(phyl_keys)):
            avg_phyl[j] += cat_ps[i][j]

    avg_phyl = avg_phyl.tolist()
    avg_phyl = [x / len(all_fun) for x in avg_phyl]

    output_string = ""
    if gain_function == "gain":
        output_string = "Number of gained gene families per function per phylostratum"
    else:
        output_string = "Number of lost gene families per function per phylostratum"

    plt.figure()
    y = cat_ps[0]
    x_ticks = []

    for j in range(len(y)):
        # label = str(j + 1) + " " + phyl_2_name[j + 1]
        # MDL changed February, 2021
        label = phyl_2_name[j + 1] + " - ps" + str(j + 1) # MDL changed February 8, 2021

        if len(label) > 40:
            label = label[:40] + "\n" + label[40:]
        x_ticks.append(label)

    for i in range(len(all_fun)):
        y = cat_ps[i]
        print("Plotting to pdf class " + str(i))
        plt.plot(x_axis, y, color='black', linewidth=0.2)

    plt.plot(x_axis, avg_phyl, color='orange', linewidth=0.8)

    plt.title(output_string, fontsize=11)
    plt.xlabel("Phylostratum", fontsize=10)
    #plt.ylabel("Gene family(cluster) count", fontsize=10)
    plt.ylabel("Gene family (cluster) count", fontsize=10) # MDL, changed February, 2021
    plt.grid(b=True, axis='both', linestyle='-', linewidth=0.6)
    # plt.xticks(np.arange(len(y)), tuple(x_ticks), fontsize=5, rotation=80)
    plt.xticks(np.arange(len(y)), tuple(x_ticks), fontsize=5, rotation=90) # MDL changed February, 2021
    plt.yticks(fontsize=6)

    plt.tight_layout()
    pdf.savefig()
    plt.close()
    pdf.close()

    pdf = matplotlib.backends.backend_pdf.PdfPages(pdf_path)
    plt.style.use('ggplot')
    x_axis = range(0, len(phyl_keys))

    for i in range(len(all_fun)):
        y = cat_ps[i]
        x_ticks = []

        for j in range(len(y)):
            # label = str(j + 1) + " " + phyl_2_name[j + 1]
            # MDL changed February, 2021
            label = phyl_2_name[j + 1] + " ps" + str(j + 1) # MDL changed February 8, 2021            
            if len(label) > 40:
                label = label[:40] + "\n" + label[40:]
            x_ticks.append(label)

        print("Plotting to pdf class " + str(i))

        output_string = ""
        # MDL changed March 22, 2021
        if all_fun[i] in go_id_2_name:
            function_name = go_id_2_name[all_fun[i]]
        else:
            function_name = "UNKNOWN_FUNCTION"			
        if gain_function == "gain":
            #output_string = "gain: " + all_fun[i] + " " + go_id_2_name[all_fun[i]]
            output_string = "GF Gain: " + all_fun[i] + " - " + function_name + " (GO)"
        else:
            #output_string = "loss: " + all_fun[i] + " " + go_id_2_name[all_fun[i]]
            output_string = "GF Loss: " + all_fun[i] + " - " + function_name + " (GO)"

        if len(output_string) > 60:
            space_index = 0
            for i in range(60, 0, -1):
                if output_string[i] == " ":
                    space_index = i
                    break
            output_string = output_string[:space_index] + "\n" + output_string[space_index + 1:]

        plt.figure()
        plt.title(output_string, fontsize=11)
        plt.plot(x_axis, y, color='black', linewidth=1.2)
        plt.xlabel("Phylostratum", fontsize=10)
        #plt.ylabel("Gene family(cluster) count", fontsize=10)
        plt.ylabel("Gene family (cluster) count", fontsize=10) # MDL changed February, 2021
        plt.grid(b=True, axis='both', linestyle='-', linewidth=0.6)
        #plt.xticks(np.arange(len(y)), tuple(x_ticks), fontsize=5, rotation=80)
        plt.xticks(np.arange(len(y)), tuple(x_ticks), fontsize=5, rotation=90) # MDL changed
        plt.yticks(fontsize=6)
        # MDL added February 18, 2021
        #plt.grid()

        plt.tight_layout()
        pdf.savefig()
        plt.close()

    pdf.close()

    print("Number of annotated clusters: " + str(annotated_clusters))


if __name__ == "__main__":

    # MDL added main arguments (March 26, 2021)
    SPECIES, C_VALUE, ALL, PLOT = parseArgv(sys.argv[1:])
    STR_PARAMETER_C = "c_0_" + str(C_VALUE)[2].rstrip()

    INPUT_CLUSTERS_FOLDER_PATH = "/export/home/mdloso/PROJECTS/ProteinsT/All_DB/Clusters_2022/"  + STR_PARAMETER_C + "/"
                                    # use for input of parsed mmseqs clusters !!! I use my results, since Tin's results had some errors
    
    
    min_max = ["gain", "loss"]
    #top_fun = [True, False]
    top_fun = [False]
    #for focal_species_key in INPUT_SPECIES:
    if True:
        focal_species = INPUT_SPECIES[SPECIES]
        print(focal_species)
        # MDL changes maded 05/02/2021
        #path_phyl_2_name = "/home/tsiroki/cluster_connection_and_funct_analysis_mdl/PS_names_upd/" + focal_species + ".txt"
        path_phyl_2_name = PHYL_NAMES_FOLDER_PATH + focal_species + ".txt"
        cluster_in_path = INPUT_CLUSTERS_FOLDER_PATH + focal_species + "_clusters.txt"
        INPUT_FOLDER_PATH = "PATH_TO_THE_ANALYSIS_FOLDER/results_2022/" + \
                                "cluster_hyper_go/" + STR_PARAMETER_C + "/"
        OUT_FOLDER_PATH = "PATH_TO_THE_ANALYSIS_FOLDER/results_2022/" + \
                                "cluster_hyper_go_summary_v2/" + STR_PARAMETER_C + "/"

        for mm in min_max:

            for top_f in top_fun:
                top_fun_string = ""
                if top_f:
                    #top_fun_string = "one_fun"
                    top_fun_string = "top_30percent_fun"
                else:
                    top_fun_string = "multi_fun"

                if ALL:
                    #in_path = INPUT_FOLDER_PATH + focal_species + "_go_full.txt"
                    in_path = INPUT_FOLDER_PATH + focal_species + "_go.txt" # MDL changed March 22, 2021
                    print("in_path: " + in_path)
                else:
                    #in_path = INPUT_FOLDER_PATH + focal_species + "_go-basic-slim-alter.txt"
                    in_path = INPUT_FOLDER_PATH + focal_species + "_go.txt"

                #out_fold = OUT_FOLDER_PATH + focal_species + "_" + mm + "\\"
                out_fold = OUT_FOLDER_PATH + focal_species + "_" + mm + "/"
                print("out_fold: " + out_fold)                
                
                # MDL changed February, 2021 --> added "_go" to all output files
                out_path = out_fold + "all_" + top_fun_string + "_" + focal_species + "_" + mm + "_go.txt"
                #out_counts_phyl = out_fold + "phyl_" + top_fun_string + "_ " + focal_species + "_" + mm + ".txt"
                out_counts_phyl = out_fold + "phyl_" + top_fun_string + "_" + focal_species + "_" + mm + "_go.txt"
                
                out_10_perc = out_fold + "phyl_" + top_fun_string + "_top_10_" + focal_species + "_" + mm + "_go.txt"
                out_all_tab = out_fold + "phyl_" + top_fun_string + "_table_" + focal_species + "_" + mm + "_go.txt"
                out_all_annotations = out_fold + "phyl_" + top_fun_string + "_annotations_" + focal_species + "_" + mm + "_go.txt"
                out_total_category_counts = out_fold + "phyl_" + top_fun_string + "_total_counts_" + focal_species + "_" + mm + "_go.txt"
                pdf_path = out_fold + "phyl_" + top_fun_string + "_all_plot_" + focal_species + "_" + mm + "_go.pdf"
                pdf_path_one = out_fold + "phyl_" + top_fun_string + "_preview_plot_" + focal_species + "_" + mm + "_go.pdf"

                make_gain_loss_count_summary(out_path, out_counts_phyl, out_10_perc, out_all_tab, out_all_annotations,
                                           out_total_category_counts, pdf_path, pdf_path_one, path_phyl_2_name, in_path,
                                           cluster_in_path, mm, plot=PLOT, only_top=top_f)