"""
The phylostratum_function_enrichment_analysis.py module performs an enrichment analysis of functions gained and lost per
phylostratum with a hypergeometric test. The standard (quant, sample, hit and total) values used in a
hypergeometric test for function enrichment per phylostratum:
•	quant		- number of annotations containing the tested function in phylostratum
•	sample		- number of annotations in phylostratum
•	hit		- total number of annotations containing the tested function
•	total		-total number of annotations

The paths and options can be set on the start of the script by changing variables, see comments for instructions.
Usage: set PARAMETER_C = 0.8, COG = True
>>> time python3 phylostratum_function_enrichment_analysis_v3.py > phylo_cog_0_8.out 

"""

__version__ = "1.0"
__author__ = "Tin Siroki"
__author__ = "Mirjana Domazet-Loso"
# changes made by MDL, 2020-2022

from os import walk
from scipy.stats import hypergeom
import statsmodels.stats.multitest as multi
import math
import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt
import numpy as np
import copy

GO_OBO_PATH = "PATH_TO_THE_ANALYSIS_FOLDER/Data/go-basic.obo"
CATEGORY_INFO_PATH = "PATH_TO_THE_ANALYSIS_FOLDER/Data/COGs.txt"

# MDL, February 2021 -- use short names
PHYL_NAMES_FOLDER_PATH = "PATH_TO_THE_ANALYSIS_FOLDER/Data/PS_names_short/"

INPUT_SPECIES = ["Hsap_2022"]

# GO
#COG = False
# COG
COG = True

# MDL, parameters for c = 0.0 to 0.8
str_parameter_c = "c_0_7"

# MDL changes made February, 2021
if COG:
	# COG
	INPUT_FOLDER_PATH = "PATH_TO_THE_ANALYSIS_FOLDER/results_2022/cluster_hyper_cog_summary_v2/" + str_parameter_c + "/"
	OUTPUT_FOLDER_PATH = "PATH_TO_THE_ANALYSIS_FOLDER/results_2022/phylostratum_hyper_cog_v2/" + str_parameter_c + "/"
	category_type = " (COG)"
	filename_end = "_cog"
else:
	# GO
	INPUT_FOLDER_PATH = "PATH_TO_THE_ANALYSIS_FOLDER/results_2022/cluster_hyper_go_summary_v2/" + str_parameter_c + "/"
	OUTPUT_FOLDER_PATH = "PATH_TO_THE_ANALYSIS_FOLDER/results_2022/phylostratum_hyper_go_v2/go-basic/" + str_parameter_c + "/"
	category_type = " (GO)"
	filename_end = "_go"

OVER_UNDER_HYPER = True
PLOT = True

# if true performs cluster function assignment by voting
# (requires Mirjana's input format, didn't have it for 2020 analysis), used it in the first batch for testing
VOTING = False


def hypergeometric_over_under(quant, sample, hit, total):
    """
     Performs hypergeometric test for overe and under expression.
     :param quant: number of successes in sample
     :param sample: sample size
     :param hit: number of successes in population
     :param total: population size
     :return: over p-value,
     """
    over_rep = hypergeom.sf(quant - 1, total, hit, sample)
    under_rep = hypergeom.cdf(quant, total, hit, sample)
    return over_rep, under_rep, min(over_rep, under_rep) * 2


def hypergeometric_over(quant, sample, hit, total):
    """
     Performs hypergeometric test for overexpression.
     :param quant: number of successes in sample
     :param sample: sample size
     :param hit: number of successes in population
     :param total: population size
     :return: over_under p-value,
     """
    over_rep = hypergeom.sf(quant - 1, total, hit, sample)
    return over_rep


def hypergeometric_under(quant, sample, hit, total):
    """
     Performs hypergeometric test for underexpression.
     :param quant: number of successes in sample
     :param sample: sample size
     :param hit: number of successes in population
     :param total: population size
     :return: under p-value,
     """
    under_rep = hypergeom.cdf(quant, total, hit, sample)
    return under_rep


# MDL changed February, 2021
# def performe_phylostratum_hyper_test(total_count_path, sample_count_path, out_txt_path,
def perform_phylostratum_hyper_test(total_count_path, sample_count_path, out_txt_path,
                                     out_plot_path, gain_loss=True, go_path=None,
                                     phyl_name_path=None, voting=False, cog=False):
    cat_2_total_count = {}
    total_count = 0
    cat_2_name = {}
    phyl_2_name = {}

    if not cog:
        if go_path:
            with open(go_path) as f:
                id = ""
                for line in f:
                    if line.startswith("id: "):
                        id = line.strip("id: ").strip()
                    if line.startswith("name: "):
                        name = line[6:].strip()
                        cat_2_name[id] = name
            cat_2_name["no_annotation"] = "Cluster without assigned function"
    else:
        with open(go_path) as f:
            cat_id = ""
            for line in f:
                if line.startswith("["):
                    line_split = line.split("]")
                    cat_id = line_split[0].strip("[")
                    name = line_split[1].strip()
                    cat_2_name[cat_id] = name


    phyl_counter = 1
    if phyl_name_path:
        with open(phyl_name_path, "r") as f:
            first = True
            for line in f:
                if first:
                    first = False
                    continue
                line_split = line.split("\t")
                phyl_2_name[phyl_counter] = line_split[1].strip()
                phyl_counter += 1

    with open(total_count_path, "r") as tot_in:
        first = True

        for line in tot_in:
            if first:
                first = False
                continue

            line_split = line.split("\t")
            cat_id = line_split[0]
            cat_tot_count = int(line_split[1])

            if cat_id in cat_2_total_count:
                raise Exception("Duplicate category id in input total counts!!!")

            cat_2_total_count[cat_id] = cat_tot_count
            total_count += cat_tot_count

    categories = []
    phyl_2_cat_counts = {}

    with open(sample_count_path, "r") as f:
        first = True

        for line in f:
            if first:
                first = False
                continue

            line_split = line.split("\t")

            cat_id = line_split[0]
            name = line_split[1]
            categories.append(cat_id)

            for i in range(2, len(line_split)):
                phyl = i - 1
                if phyl not in phyl_2_cat_counts:
                    phyl_2_cat_counts[phyl] = []

                phyl_2_cat_counts[phyl].append(int(line_split[i]))

    ids_remove = []
    for ph in phyl_2_cat_counts:
        if sum(phyl_2_cat_counts[ph]) == 0:
            ids_remove.append(ph)

    for id_r in ids_remove:
        phyl_2_cat_counts.pop(id_r)

    p_values = []

    for phyl in phyl_2_cat_counts:
        cat_counts = phyl_2_cat_counts[phyl]
        print("Calculating hyper test for phylostratum " + str(phyl))

        sample_count = sum(cat_counts)

        for i in range(len(cat_counts)):
            kvant = cat_counts[i]
            cat_id = categories[i]

            if OVER_UNDER_HYPER:
                over, under, p_value = hypergeometric_over_under(kvant, sample_count, cat_2_total_count[cat_id],
                                                             total_count)
            else:
                p_value = hypergeometric_over(kvant, sample_count, cat_2_total_count[cat_id],
                                                             total_count)
            p_values.append(p_value)

    p_adjusted = multi.multipletests(p_values, method="fdr_bh")[1]
    p_counter = 0

    cat_2_plot_data = {}

    with open(out_txt_path, "w") as out:

        for phyl in phyl_2_cat_counts:
            cat_counts = phyl_2_cat_counts[phyl]
            print("Storing phylostratum " + str(phyl))
            out.write("\nPhylostratum " + str(phyl) + "\n")

            sample_count = sum(cat_counts)

            for i in range(len(cat_counts)):
                kvant = cat_counts[i]
                cat_id = categories[i]

                if cat_id not in cat_2_plot_data:
                    cat_2_plot_data[cat_id] = {}
                    cat_2_plot_data[cat_id]["log_ratio"] = []
                    cat_2_plot_data[cat_id]["p_value"] = []
                    cat_2_plot_data[cat_id]["p_adj"] = []
                    if not voting:
                        cat_2_plot_data[cat_id]["name"] = cat_2_name[cat_id]
                    else:
                        cat_2_plot_data[cat_id]["name"] = ""
                        cat_2_name[cat_id] = ""
                    cat_2_plot_data[cat_id]["odds_ratio"] = []
                    cat_2_plot_data[cat_id]["kvant_sum"] = 0
                    cat_2_plot_data[cat_id]["total_count"] = cat_2_total_count[cat_id]
                    cat_2_plot_data[cat_id]["total_sample"] = total_count
                    cat_2_plot_data[cat_id]["sample_sum"] = 0
                    cat_2_plot_data[cat_id]["odds_sample"] = []

                quant = kvant
                sample = sample_count
                hit = cat_2_total_count[cat_id]
                total = total_count

                if quant != 0 and (sample - quant) != 0 and (total - hit - sample + quant) != 0 and (hit - quant) != 0:
                    odds_sample = quant / (sample - quant)
                    odds_rest = (hit - quant) / (total - hit - sample + quant)
                    real_log_odds = math.log2(odds_sample / odds_rest)
                else:
                    real_log_odds = "nan"

                if not OVER_UNDER_HYPER and not real_log_odds == "nan":
                    if real_log_odds < 0:
                        real_log_odds = 0

                p_value = p_values[p_counter]
                p_adj = p_adjusted[p_counter]

                out_line = cat_id + "\t" + cat_2_name[cat_id] + "\t" + str(kvant) + "\t" \
                           + str(sample_count) + "\t" + str(cat_2_total_count[cat_id]) \
                           + "\t" + str(total_count) + "\t" + str(p_value) + "\t" + str(p_adj)\
                           + "\t" + str(real_log_odds)

                out.write(out_line + "\n")

                # adjust real log odds for plotting.
                if real_log_odds == "nan":
                    if quant == 0 or (total - hit - sample + quant) == 0:
                        real_log_odds = -4
                    if sample - quant == 0 or hit - quant == 0:
                        real_log_odds = 4

                cat_2_plot_data[cat_id]["log_ratio"].append(real_log_odds)
                cat_2_plot_data[cat_id]["p_value"].append(p_value)
                cat_2_plot_data[cat_id]["p_adj"].append(p_adj)
                cat_2_plot_data[cat_id]["kvant_sum"] += kvant
                cat_2_plot_data[cat_id]["sample_sum"] += sample_count
                p_counter += 1

    if not PLOT:
        return

    pdf = matplotlib.backends.backend_pdf.PdfPages(out_plot_path)
    plt.style.use('ggplot')
    phyl_keys = phyl_2_cat_counts.keys()
    phyl_keys = list(phyl_keys)
    phyl_keys.sort()

    offset = phyl_keys[0]

    x_axis = range(len(phyl_keys))
    category_counter = 0

    for cat in cat_2_plot_data:

        y_val = cat_2_plot_data[cat]["log_ratio"]
        x_ticks = []

        if cat_2_plot_data[cat]["kvant_sum"] != cat_2_plot_data[cat]["total_count"]:
            print("Duplicate sampling.")
        if cat_2_plot_data[cat]["sample_sum"] != cat_2_plot_data[cat]["total_sample"]:
            print("Duplicate sampling.")

        for j in range(len(y_val)):
            #label = str(j + offset) + " " + phyl_2_name[j + offset]
            label = phyl_2_name[j + offset] + " - ps" + str(j + offset) # MDL changed February, 2021
            if len(label) > 40:
                label = label[:40] + "\n" + label[40:]
            x_ticks.append(label)

        print("Plotting to pdf category " + str(category_counter))

        annotations = []
        counter = 0
        for p_adj in cat_2_plot_data[cat]["p_adj"]:
            annot = ""
            if p_adj < 0.05:
                annot = "*"
            if p_adj < 0.01:
                annot = "**"
            if p_adj < 0.001:
                annot = "***"
            if cat_2_plot_data[cat]["log_ratio"][counter] == -4:
                annot += "-inf"
            if cat_2_plot_data[cat]["log_ratio"][counter] == 4:
                annot += "inf"

            if not OVER_UNDER_HYPER and cat_2_plot_data[cat]["log_ratio"][counter] <= 0:
                annot = ""

            annotations.append(annot)

            counter += 1

        plt.figure()

        if gain_loss:
            #output_string = "gain: " + cat + " " + cat_2_plot_data[cat]["name"]
            output_string = "GF Gain: " + cat + " - " + cat_2_plot_data[cat]["name"] + category_type # MDL changed February, 2021
        else:
            #output_string = "loss: " + cat + " " + cat_2_plot_data[cat]["name"]
            output_string = "GF Loss: " + cat + " - " + cat_2_plot_data[cat]["name"] + category_type # MDL changed February, 2021

        if len(output_string) > 60:
            space_index = 0
            for i in range(60, 0, -1):
                if output_string[i] == " ":
                    space_index = i
                    break
            output_string = output_string[:space_index] + "\n" + output_string[space_index+1:]

        plt.title(output_string, fontsize=11)
        plt.plot(x_axis, y_val, color='black', linewidth=1.2)
        plt.axhline(linewidth=0.5, color='black')

        lab_counter = 0
        for x, y in zip(x_axis, y_val):
            label = annotations[lab_counter]
            plt.annotate(label,  # this is the text
                         (x, y),  # this is the point to label
                         textcoords="offset points",  # how to position the text
                         xytext=(0, 2),  # distance from text to points (x,y)
                         ha='center',
                         fontsize=5)

            lab_counter += 1

        plt.xlabel("Phylostratum", fontsize=10)
        #plt.ylabel("Log ratio", fontsize=10)
        plt.ylabel("Log odds ratio", fontsize=10)
        plt.grid(b=True, axis='both', linestyle='-', linewidth=0.6) # MDL changed February, 2021

        #plt.xticks(np.arange(len(y_val)), tuple(x_ticks), fontsize=5, rotation=80)
        plt.xticks(np.arange(len(y_val)), tuple(x_ticks), fontsize=5, rotation=90) # MDL changed February 8, 2021
        plt.yticks(fontsize=6)

        plt.tight_layout()
        pdf.savefig()
        plt.close()
        category_counter += 1

    pdf.close()
    return copy.deepcopy(cat_2_plot_data)


def plot_gain_loss(gain, loss, out_path, phyl_name_path=None):
    """
    Hypergeometric test for functions in phylostrata (gain and loss).

    :param gain:
    :param loss:
    :param out_path:
    :param phyl_name_path:
    :return:
    """
    phyl_2_name = {}

    pdf = matplotlib.backends.backend_pdf.PdfPages(out_path)
    plt.style.use('ggplot')
    fun_keys = list(gain.keys())
    first_fun = gain[fun_keys[0]]

    phyl_keys = range(1, len(first_fun["log_ratio"]) + 1)
    phyl_keys = list(phyl_keys)

    offset = phyl_keys[0]

    x_axis = range(len(phyl_keys))
    category_counter = 0


    phyl_counter = 1
    if phyl_name_path:
        with open(phyl_name_path, "r") as f:
            first = True
            for line in f:
                if first:
                    first = False
                    continue
                line_split = line.split("\t")
                phyl_2_name[phyl_counter] = line_split[1].strip()
                phyl_counter += 1

    for cat in gain:

        if cat not in loss:
            continue

        gain_length = len(gain[cat]["log_ratio"])
        loss_length = len(loss[cat]["log_ratio"])

        off_phyl = gain_length - loss_length

        loss[cat]["log_ratio"] = ([0] * off_phyl) + loss[cat]["log_ratio"]
        #loss[cat]["log_ratio"] = ([""] * off_phyl) + loss[cat]["log_ratio"] # MDL changed February 9, 2021
        loss[cat]["p_adj"] = ([1] * off_phyl) + loss[cat]["p_adj"]

        y_val_gain = gain[cat]["log_ratio"]
        y_val_loss = loss[cat]["log_ratio"]
        x_ticks = []

        for j in range(len(y_val_gain)):
            #label = str(j + offset) + " " + phyl_2_name[j + offset]
            label = phyl_2_name[j + offset] + " - ps" + str(j + offset) # MDL changed February, 2021
            if len(label) > 40:
                label = label[:40] + "\n" + label[40:]
            x_ticks.append(label)

        print("Plotting to pdf category " + str(category_counter))

        annotations = []
        counter = 0
        for p_adj in gain[cat]["p_adj"]:
            annot = ""
            if p_adj < 0.05:
                annot = "*"
            if p_adj < 0.01:
                annot = "**"
            if p_adj < 0.001:
                annot = "***"
            if gain[cat]["log_ratio"][counter] == -4:
                annot += "-inf"

            if not OVER_UNDER_HYPER and gain[cat]["log_ratio"][counter] <= 0:
                annot = ""
            annotations.append(annot)
            counter += 1

        counter = 0
        for p_adj in loss[cat]["p_adj"]:
            annot = ""
            if p_adj < 0.05:
                annot = "*"
            if p_adj < 0.01:
                annot = "**"
            if p_adj < 0.001:
                annot = "***"
            if loss[cat]["log_ratio"][counter] == -4:
                annot += "-inf"

            if not OVER_UNDER_HYPER and loss[cat]["log_ratio"][counter] <= 0:
                annot = ""
            annotations.append(annot)
            counter += 1

        plt.figure()

        #output_string = "gain loss " + cat + " " + gain[cat]["name"]
        output_string = "GF Gain/Loss: " + cat + " - " + gain[cat]["name"] + category_type # MDL changed February, 2021

        if len(output_string) > 60:
            space_index = 0
            for i in range(60, 0, -1):
                if output_string[i] == " ":
                    space_index = i
                    break
            output_string = output_string[:space_index] + "\n" + output_string[space_index+1:]

        plt.title(output_string, fontsize=11)
        #plt.plot(x_axis, y_val_gain, color='green', linewidth=1.2, label="gain") # MDL changed 8/2/2021
        plt.plot(x_axis, y_val_gain, color='blue', linewidth=1.2, label="gain")
        #plt.plot(x_axis, y_val_loss, color='red', linewidth=1.2, label="loss")
        plt.plot(x_axis[2:], y_val_loss[2:], color='red', linewidth=1.2, label="loss") # MDL changed 9/2/2021
        plt.axhline(linewidth=0.5, color='black')

        lab_counter = 0
        for x, y in zip(x_axis, y_val_gain):
            label = annotations[lab_counter]

            plt.annotate(label,  # this is the text
                         (x, y),  # this is the point to label
                         textcoords="offset points",  # how to position the text
                         # MDL, February 2021; color added
                         color='blue',
                         xytext=(0, 2),  # distance from text to points (x,y)
                         ha='center',
                         fontsize=6)# MDL changed 9/2/2021, fonst size 5->6

            lab_counter += 1
        
        for x, y in zip(x_axis, y_val_loss):
            label = annotations[lab_counter]
            
            plt.annotate(label,  # this is the text
                         (x, y),  # this is the point to label
                         textcoords="offset points",  # how to position the text
                         # MDL, February 2021; color added
                         color='red',
                         xytext=(0, 2),  # distance from text to points (x,y)
                         ha='center',
                         fontsize=6)# MDL changed 9/2/2021, fonst size 5->6
            lab_counter += 1


        plt.xlabel("Phylostratum", fontsize=10)
        #plt.ylabel("Log ratio", fontsize=10)
        plt.ylabel("Log odds ratio", fontsize=10)
        plt.grid(b=True, axis='both', linestyle='-', linewidth=0.6) # MDL changed February, 2021

        #plt.xticks(np.arange(len(y_val_gain)), tuple(x_ticks), fontsize=5, rotation=80)
        plt.xticks(np.arange(len(y_val_gain)), tuple(x_ticks), fontsize=5, rotation=90) # MDL changed February 8, 2021
        plt.yticks(fontsize=6)
        #plt.legend(loc="lower right")       
        plt.legend(loc="lower right", prop={'size': 7}) # MDL added prop/font size February 8, 2021
        
        plt.tight_layout()
        pdf.savefig()
        plt.close()
        category_counter += 1

    pdf.close()


if __name__ == "__main__":
    gain_loss = [True, False]
    if VOTING:
        multi_fun = ["voting", "multi_fun", "one_fun"]
    else:
        #multi_fun = ["multi_fun", "one_fun"]
        #multi_fun = ["multi_fun", "top_30percent_fun"]
        multi_fun = ["multi_fun"]

    if COG:
        cat_path = CATEGORY_INFO_PATH
    else:
        cat_path = GO_OBO_PATH
    gain = None
    loss = None
    
    for focal_species in INPUT_SPECIES:
        path_phyl_2_name = PHYL_NAMES_FOLDER_PATH + focal_species + ".txt"

        for mul in multi_fun:

            for g_l in gain_loss:

                print(focal_species)
                print("App: " + str(g_l))
                print("Multi: " + mul)

                if g_l:
                    total_count_path = INPUT_FOLDER_PATH + focal_species + "_gain/phyl_" + mul + "_total_counts_"\
                                       + focal_species + "_gain" + filename_end + ".txt"
                    sample_count_path = INPUT_FOLDER_PATH + focal_species + "_gain/phyl_" + mul + "_table_" \
                                        + focal_species + "_gain" + filename_end + ".txt"
                    out_path_txt = OUTPUT_FOLDER_PATH + focal_species + "/" + mul + "_gain_" \
                                   + focal_species + filename_end + ".txt"
                    out_path_plot = OUTPUT_FOLDER_PATH + focal_species + "/" + mul + "_gain_plot_" \
                                    + focal_species + filename_end + ".pdf"
                else:
                    total_count_path = INPUT_FOLDER_PATH + focal_species + "_loss/phyl_" + mul + "_total_counts_" \
                                       + focal_species + "_loss" + filename_end + ".txt"
                    sample_count_path = INPUT_FOLDER_PATH + focal_species + "_loss/phyl_" + mul + "_table_" \
                                        + focal_species + "_loss" + filename_end + ".txt"
                    out_path_txt = OUTPUT_FOLDER_PATH + focal_species + "/" + mul + "_loss_" \
                                   + focal_species + filename_end + ".txt"
                    out_path_plot = OUTPUT_FOLDER_PATH + focal_species + "/" + mul + "_loss_plot_" \
                                    + focal_species + filename_end + ".pdf"

                if mul == "voting":
                    if g_l:
                        gain = perform_phylostratum_hyper_test(total_count_path, sample_count_path, out_path_txt,
                                                                out_path_plot, gain_loss=g_l, go_path=cat_path,
                                                                phyl_name_path=path_phyl_2_name, voting=True, cog=COG)
                    else:
                        loss = perform_phylostratum_hyper_test(total_count_path, sample_count_path, out_path_txt,
                                                                out_path_plot, gain_loss=g_l, go_path=cat_path,
                                                                phyl_name_path=path_phyl_2_name, voting=True, cog=COG)

                else:
                    if g_l:
                        gain = perform_phylostratum_hyper_test(total_count_path, sample_count_path, out_path_txt,
                                                                out_path_plot, gain_loss=g_l, go_path=cat_path,
                                                                phyl_name_path=path_phyl_2_name, cog=COG)
                    else:
                        loss = perform_phylostratum_hyper_test(total_count_path, sample_count_path, out_path_txt,
                                                         out_path_plot, gain_loss=g_l, go_path=cat_path,
                                                         phyl_name_path=path_phyl_2_name, cog=COG)

#            plot_gain_loss(gain, loss, OUTPUT_FOLDER_PATH + focal_species + "\\" + mul + "_gain_loss_plot_" \
            plot_gain_loss(gain, loss, OUTPUT_FOLDER_PATH + focal_species + "/" + mul + "_gain_loss_plot_" \
                                    + focal_species + filename_end + ".pdf", phyl_name_path=path_phyl_2_name)
