"""
The cluster_function_enrichment_summary_cog.py module performs processing of function enrichment per cluster for COG terms.
The module gives a function (one_fun) or functions (multi_fun) to a cluster.
A function is assigned to a cluster if its hypergeometric p-adjusted value is lower than a specified cut-off value
(P_ADJ_CUTOFF) and half or more genes in the cluster are annotated with the corresponding function.
The paths and options can be set on the start of the script by changing variables, see comments for instructions.
"""
# Usage: time python3 cluster_function_enrichment_summary_cog_v2.py > cog_v2.out

__version__ = "1.0"
__author__ = "Tin Siroki"
__author__ = "Mirjana Domazet-Loso" # new versions

import operator
import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt
import numpy as np

P_ADJ_CUTOFF = 0.05
# MDL changes made 28/01/2021
CATEGORY_INFO_PATH = "PATH_TO_THE_ANALYSIS_FOLDER/Data/COGs.txt"

# MDL, parameters for c = 0.0 and c = 0.4
#PARAMETER_C = 0.0
str_parameter_c = "c_0_7"

INPUT_CLUSTERS_FOLDER_PATH = "PATH_TO_THE_ANALYSIS_FOLDER/Clusters_2022/" + str_parameter_c	+ "/"
# use for input of cog analysis results 
INPUT_FOLDER_PATH = "PATH_TO_THE_ANALYSIS_FOLDER/results_2022/cluster_hyper_cog/" + str_parameter_c	+ "/" 							
OUT_FOLDER_PATH = "PATH_TO_THE_ANALYSIS_FOLDER/results_2022/cluster_hyper_cog_summary_v2/" + str_parameter_c + "/"

# MDL, February 2021 -- use short names
PHYL_NAMES_FOLDER_PATH = "PATH_TO_THE_ANALYSIS_FOLDER/Data/PS_names_short/"

INPUT_SPECIES = ["Hsap_2022", "Athal_2022", "Dmel_2022", "Scer_2022"]

INCLUDE_NO_ANNOTATION = True

# if true performs cluster function assignment by voting
# (requires Mirjana's input format, didn't have it for 2020 analysis), used it in the first batch for testing
VOTING = False


#def get_cluster_top_function(functions, gene_count):
#    """
#    Returns the best significant function that is present in more than half of the genes in the cluster (lowest p-adj)
#    function per cluster, empty list is returned if cluster does not contain significant functions.
#    :param functions:
#    :param gene_count:
#    :return:
#    """
#    min_p_adj = float('+inf')
#    top_function = ""
#    for fun in functions:
#        if fun["p_adj"] < P_ADJ_CUTOFF:
#
#            if fun["p_adj"] < min_p_adj and fun["hit"] >= (gene_count/2):
#                top_function = fun["info_all"]
#                min_p_adj = fun["p_adj"]
#
#    res = []
#    if top_function != "":
#        res.append(top_function)
#    return res
#
#
#def get_all_clusters_function(functions, gene_count):
#    """
#    Returns all significant functions that are present in more than half of the genes in the cluster (lowest p-adj)
#    function per cluster, empty list is returned if cluster does not contain significant functions.
#    :param functions:
#    :param gene_count:
#    :return:
#    """
#    fun_filtered = []
#    for fun in functions:
#        if fun["p_adj"] < P_ADJ_CUTOFF and fun["hit"] >= (gene_count/2):
#            fun_filtered.append(fun)
#
#    fun_filtered.sort(key=operator.itemgetter("p_adj"))
#    res = []
#    for el in fun_filtered:
#        res.append(el["info_all"])
#    return res

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


def get_name_string(cat_id, cat_2_name):
    name_string = ""
    cats = cat_id.split(", ")

    for cat in cats:
        name = cat_2_name[cat]
        name_string += name + ", "

    return name_string.strip(", ")


def make_gain_loss_count_summary(out_path, out_counts_phyl, out_10_perc, out_all_tab, out_all_annotations,
                               out_total_category_counts, pdf_path, pdf_path_one, path_phyl_2_name,
                               in_path, cluster_in_path, gain_function, orig_format=True, only_top=True):

    cluster_2_info = {}
    cluster_2_hyper = {}
    cat_id_2_name = {}
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
        cat_id = ""
        print("Reading categories.")
        for line in f:
            if line.startswith("["):
                line_split = line.split("]")
                cat_id = line_split[0].strip("[")
                name = line_split[1].strip()
                cat_id_2_name[cat_id] = name

    with open(cluster_in_path, "r") as f:

        for line in f:
            if line == "" or line == "\r" or line == "\n":
                continue

            line = line.strip()
            line_split = line.split()

            if len(line_split) < 2:
                continue

            if line_split[1].startswith("pgi"):
                cl_in = {}
                cl_in["id_line"] = "\t".join(line_split[1:]).strip()
                cl_in["ps_start"] = int(line_split[3])
                #cl_in["ps_end"] = int(line_split[4]) + 1
                if line_split[4] == "-": cl_in["ps_end"] = "-"
                else: cl_in["ps_end"] = int(line_split[4])	# MDL changed 28/01/2021: removed (+ 1)				
                cl_in["gene_count"] = int(line_split[5].strip())
                cluster_2_info[line_split[1]] = cl_in.copy()

    print("Number of clusters: " + str(len(cluster_2_info)))

    if orig_format:
        with open(in_path, "r") as f:
            cluster_id = ""
            for line in f:
                if line == "" or line == "\r" or line == "\n":
                    continue
                line_split = line.split("\t")

                if line_split[0].startswith("pgi"):
                    cluster_id = line_split[0].strip()
                    cluster_2_hyper[cluster_id] = []
                else:
                    fun = {}
                    fun["id"] = line_split[0]
                    fun["hit"] = int(line_split[1])
                    # MDl changes made March 22, 2021
                    #fun["p_value"] = float(line_split[5])
                    #fun["p_adj"] = float(line_split[6])
                    #fun["log_ratio"] = float(line_split[7].strip())
                    #fun["info_all"] = line.strip()
                    fun["log_ratio"] = float(line_split[5].strip())
                    fun["info_all"] = line.strip()

                    cluster_2_hyper[cluster_id].append(fun.copy())
    else:
        with open(in_path, "r") as f:
            for line in f:
                if line == "" or line == "\r" or line == "\n":
                    continue

                if "pgi" in line:
                    line_split = line.split("\t")
                    cluster_id = line_split[1]
                    category = line_split[6].strip("[")
                    category = category.strip()
                    category = category.strip("*")
                    category = category.strip()
                    category = category.strip("]")
                    cluster_2_hyper[cluster_id] = category


    print("Number of clusters hyper: " + str(len(cluster_2_hyper)))

    annotated_clusters = 0
    no_annot_cluster_count = 0
    phyl_2_cat = {}
    all_fun = set()

    with open(out_path, "w") as out:
        with open(out_all_annotations, "w") as out_ann:
            print("Writing out_path and out_all_annotations" + out_path + " " + out_all_annotations)
            first = True

            for cl_id in cluster_2_hyper:

                if orig_format:
                    if only_top:
                        # top_functions = get_cluster_top_function(cluster_2_hyper[cl_id],
                        # MDL changed February, 2021
                        top_functions = get_cluster_top_function_v2(cluster_2_hyper[cl_id],
                                                                 cluster_2_info[cl_id]["gene_count"])
                    else:
                        # top_functions = get_all_clusters_function(cluster_2_hyper[cl_id],
                        # MDL changed February, 2021
                        top_functions = get_all_clusters_function_v2(cluster_2_hyper[cl_id],
                                                                  cluster_2_info[cl_id]["gene_count"])
                else:
                    top_functions = [cluster_2_hyper[cl_id]]

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
                if len(top_functions) > 0:
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

                        cat_2_total_count[function_id] += 1
                        total_count += 1
                        out_ann.write(cl_id + "###" + function_id + "\t" + str(ps) + "\n")

                        if ps not in phyl_2_cat:
                            phyl_2_cat[ps] = {}

                        if function_id not in phyl_2_cat[ps]:
                            phyl_2_cat[ps][function_id] = 0

                        phyl_2_cat[ps][function_id] += 1
                        out.write(fun + "\n")

    print("Number of cluster without annotations: " + str(no_annot_cluster_count))

    for i in range(1, max_ps+1):
        if i not in phyl_2_cat:
            phyl_2_cat[i] = {}

    tot_counts = []
    with open(out_total_category_counts, "w") as out:
        print("Writing out_total_category_counts")
        out.write("category\ttotal_count\n")
        for cat in cat_2_total_count:
            line = cat + "\t" + str(cat_2_total_count[cat])
            tot_counts.append(cat_2_total_count[cat])
            out.write(line + "\n")

    if sum(tot_counts) != total_count:
        raise Exception("Total counts error!!!")

    print("Total annotation count: " + str(total_count))

    with open(out_counts_phyl, "w") as out:	
        with open(out_10_perc, "w") as out_10:
            print("Writing out_counts_phyl and out_10_perc.")
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
                    if orig_format:
                        out.write(el[0] + "\t" + cat_id_2_name[el[0]] + "\t" + str(el[1]) + "\n")
                    else:
                        out.write(el[0] + "\t" + get_name_string(el[0], cat_id_2_name) + "\t" + str(el[1]) + "\n")

                top_len = int(round(len(counts_sorted) * 0.1))

                out_10.write("Phylostratum " + str(ps) + "\tcat_count=" + str(top_len) + "\ttotal=" + str(
                    total) + ":\n")

                for el in counts_sorted[:top_len]:
                    if orig_format:
                        out_10.write(el[0] + "\t" + cat_id_2_name[el[0]] + "\t" + str(el[1]) + "\n")
                    else:
                        out_10.write(el[0] + "\t" + get_name_string(el[0], cat_id_2_name) + "\t" + str(el[1]) + "\n")

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
            if orig_format:
                line = all_fun[i] + "\t" + cat_id_2_name[all_fun[i]] + "\t"
            else:
                line = all_fun[i] + "\t" + get_name_string(all_fun[i], cat_id_2_name) + "\t"

            for j in range(len(phyl_keys)):
                line += str(cat_ps[i][j]) + "\t"
            line = line.strip()

            out.write(line + "\n")

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
        #label = str(j + 1) + " " + phyl_2_name[j + 1]
        # MDL changed February, 2021
        #label = phyl_2_name[j + 1] + " ps" + str(j + 1) # MDL changed February 8, 2021
        label = phyl_2_name[j + 1] + " - ps" + str(j + 1) # MDL changed February 18, 2021
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
    #plt.xticks(np.arange(len(y)), tuple(x_ticks), fontsize=5, rotation=80)
    plt.xticks(np.arange(len(y)), tuple(x_ticks), fontsize=5, rotation=90) # MDL, changed February, 2021
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
            #label = str(j+1) + " " + phyl_2_name[j+1]
            # MDL changed February, 2021
            #label = phyl_2_name[j + 1] + " ps" + str(j + 1) # MDL changed February 8, 2021
            label = phyl_2_name[j + 1] + " - ps" + str(j + 1) # MDL changed February 8, 2021
            if len(label) > 40:
                label = label[:40] + "\n" + label[40:]
            x_ticks.append(label)

        print("Plotting to pdf class " + str(i))

        output_string = ""

        if gain_function == "gain": # MDL changes made February, 2021
            if orig_format:
                #output_string = "Gain: " + all_fun[i] + " " + cat_id_2_name[all_fun[i]]
                output_string = "GF Gain: " + all_fun[i] + " - " + cat_id_2_name[all_fun[i]] + " (COG)"
            else:
                #output_string = "Gain: " + all_fun[i]# + " " + get_name_string(all_fun[i], cat_id_2_name)
                output_string = "GF Gain: " + all_fun[i] + " (COG)"# + " " + get_name_string(all_fun[i], cat_id_2_name)
        else:
            if orig_format:
                #output_string = "Loss: " + all_fun[i] + " " + cat_id_2_name[all_fun[i]]
                output_string = "GF Loss: " + all_fun[i] + " - " + cat_id_2_name[all_fun[i]] + " (COG)"
            else:
                #output_string = "Loss: " + all_fun[i]# + " " + get_name_string(all_fun[i], cat_id_2_name)
                output_string = "GF Loss: " + all_fun[i] + " (COG)" # + " " + get_name_string(all_fun[i], cat_id_2_name)



        if len(output_string) > 60:
            space_index = 0
            for i in range(60, 0, -1):
                if output_string[i] == " ":
                    space_index = i
                    break
            output_string = output_string[:space_index] + "\n" + output_string[space_index+1:]

        plt.figure()
        plt.title(output_string, fontsize=11)
        plt.plot(x_axis, y, color='black', linewidth=1.2)
        plt.xlabel("Phylostratum", fontsize=10)
        #plt.ylabel("Gene family(cluster) count", fontsize=10)
        plt.ylabel("Gene family (cluster) count", fontsize=10) # MDL changed February, 2021
        plt.grid(b=True, axis='both', linestyle='-', linewidth=0.6)
        #plt.xticks(np.arange(len(y)), tuple(x_ticks), fontsize=5, rotation=80)
        plt.xticks(np.arange(len(y)), tuple(x_ticks), fontsize=5, rotation=90) # MDL changed February, 2021
        plt.yticks(fontsize=6)
        plt.tight_layout()
        pdf.savefig()
        plt.close()

    pdf.close()
    print("Number of annotated clusters: " + str(annotated_clusters))



if __name__ == "__main__":

    min_max = ["gain", "loss"]
    top_fun = [True, False]
    for focal_species in INPUT_SPECIES:
        print(focal_species)
		# MDL changes maded 28/01/2021
        path_phyl_2_name = PHYL_NAMES_FOLDER_PATH + focal_species + ".txt"
        cluster_in_path = INPUT_CLUSTERS_FOLDER_PATH + focal_species + "_clusters.txt"

        for mm in min_max:
            # Mirjana variant
            if VOTING:
                in_path = INPUT_FOLDER_PATH + focal_species + "_annot_dbAllPlus_clu2.txt"
                #out_fold = OUT_FOLDER_PATH + focal_species + "_" + mm + "\\"
                out_fold = OUT_FOLDER_PATH + focal_species + "_" + mm + "/"

                out_path = out_fold + "all_voting_" + focal_species + "_" + mm + "_cog.txt"
                out_counts_phyl = out_fold + "phyl_voting" + focal_species + "_" + mm + "_cog.txt"
                out_10_perc = out_fold + "phyl_voting_top_10_" + focal_species + "_" + mm + "_cog.txt"
                out_all_tab = out_fold + "phyl_voting_all_" + focal_species + "_" + mm + "_cog.txt"
                out_all_annotations = out_fold + "phyl_voting_annotations_" + focal_species + "_" + mm + "_cog.txt"
                out_total_category_counts = out_fold + "phyl_voting_total_counts_" + focal_species + "_" + mm + "_cog.txt"
                pdf_path = out_fold + "phyl_voting_all_plot_" + focal_species + "_" + mm + "_cog.pdf"
                pdf_path_one = out_fold + "phyl_voting_preview_plot_" + focal_species + "_" + mm + "_cog.pdf"

                make_gain_loss_count_summary(out_path, out_counts_phyl, out_10_perc, out_all_tab, out_all_annotations,
                                             out_total_category_counts, pdf_path, pdf_path_one, path_phyl_2_name,
                                             in_path,
                                             cluster_in_path, mm, orig_format=False)

            for top_f in top_fun:
                top_fun_string = ""
                if top_f:
                    #top_fun_string = "one_fun"
                    top_fun_string = "top_30percent_fun" # MDL changed February, 2021
                else:
                    top_fun_string = "multi_fun"

                in_path = INPUT_FOLDER_PATH + focal_species + "_cog.txt"
                out_fold = OUT_FOLDER_PATH + focal_species + "_" + mm + "/"

                # MDL changed February, 2021 --> added "cog" to all output files
                out_path = out_fold + "all_" + top_fun_string + "_" +focal_species + "_" + mm + "_cog.txt"
                out_counts_phyl = out_fold + "phyl_" + top_fun_string + "_" + focal_species + "_" + mm + "_cog.txt"
                out_10_perc = out_fold + "phyl_" + top_fun_string + "_top_10_" + focal_species + "_" + mm + "_cog.txt"
                out_all_tab = out_fold + "phyl_" + top_fun_string + "_table_" + focal_species + "_" + mm + "_cog.txt"
                out_all_annotations = out_fold + "phyl_" + top_fun_string + "_annotations_" + focal_species + "_" + mm + "_cog.txt"
                out_total_category_counts = out_fold + "phyl_" + top_fun_string + "_total_counts_" + focal_species + "_" + mm + "_cog.txt"
                pdf_path = out_fold + "phyl_" + top_fun_string + "_all_plot_" + focal_species + "_" + mm + "_cog.pdf"
                pdf_path_one = out_fold + "phyl_" + top_fun_string + "_preview_plot_" + focal_species + "_" + mm + "_cog.pdf"

                make_gain_loss_count_summary(out_path, out_counts_phyl, out_10_perc, out_all_tab, out_all_annotations,
                                           out_total_category_counts, pdf_path, pdf_path_one, path_phyl_2_name, in_path,
                                           cluster_in_path, mm, only_top=top_f)
