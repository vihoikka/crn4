import pandas as pd
import os
import argparse
import re
from pprint import pprint
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("-on", "--out_nodes", required=True)  # ok
parser.add_argument("-oe", "--out_edges", required=True)  # ok
parser.add_argument("-om", "--out_merged_all", required=True)  # ok
parser.add_argument("-p", "--phage_data", required=True)  # ok
parser.add_argument("-bp", "--base_path", required=True)  # ok
parser.add_argument("-ah", "--all_hosts", required=True)  # ok
parser.add_argument("-hm", "--host_mastertable", required=True)  # ok
parser.add_argument("-sb", "--spacer_blast_hits", required=True)  # ok
parser.add_argument("-of", "--outfolder", required=True)  # ok
parser.add_argument("-pr", "--phage_rns", required=True)  # ok

# Generate shell command for running the script using the long args (e.g. --out_nodes)
# python3 spacer_hit_analyzer.py --out_nodes nodes.tsv --out_edges edges.tsv --out_merged_all merged_all.tsv --phage_data phage_metadata.tsv --base_path /path/to/ --all_hosts all_hosts.tsv --host_mastertable host_mastertable.tsv --spacer_blast_hits spacer_blast_hits.tsv --outfolder /path/to/


args = parser.parse_args()
out_nodes = args.out_nodes  # ok
out_edges = args.out_edges  # ok
out_merged_all = args.out_merged_all  # ok
phage_metadata = args.phage_data  # ok
base_path = args.base_path  # ok
all_hosts = args.all_hosts  # ok
host_mastertable = args.host_mastertable  # ok
spacer_blast_hits = args.spacer_blast_hits  # ok
outfolder = args.outfolder  # ok
phage_rns = args.phage_rns  # ok

phages = pd.read_csv(
    phage_metadata, sep="\t", header=0
)  # read Milliard phage metadata into panda df, separator tab
blast = pd.read_csv(
    spacer_blast_hits, sep="\t", header=0
)  # read blast hits into panda df, separator tab
bacteria = pd.read_csv(
    all_hosts, sep="\t", header=None
)  # list of all bacteria in dataset regardless of interactions
host_mastertable = pd.read_csv(
    host_mastertable, sep="\t", header=0
)  # read host mastertable into panda df, separator tab
phage_rn = pd.read_csv(
    phage_rns, sep="\t", header=0
)  # this contains phage RN hmm search analysed data
bacteria.columns = ["bacterium"]

### Merging with host data ###

# In the blast hits, each row is a spacer hit. We want to merge this with the host data to get more information on the hosts

# parsing the host column in blast hits. Format example: NC_010803.1_1:21@GCF_000020465.1
blast["host_nt_spacer"] = blast["Query_id"].str.extract(
    r"(.+?(?=@))", expand=False
)  # NC_010803.1_1:21@GCF_000020465.1 -> NC_010803.1_1:21
blast["host_nt"] = blast["host_nt_spacer"].str.extract(
    r"(.*(?=\_))", expand=False
)  # NC_010803.1_1:21 -> NC_010803.1

# make host column, corresponding to our host wildcards
blast["host"] = blast["Query_id"].str.extract(
    r"((?<=\@).*)", expand=False
)  # NC_010803.1_1:21@GCF_000020465.1 -> GCF_000020465.1

# merge host_mastertable and blast hits. Hosts with no spacer hits are also included due to outer join. add suffixes to overlapping columns
merged_1 = pd.merge(
    blast,
    host_mastertable,
    left_on="host",
    right_on="sample",
    how="outer",
    suffixes=("_phage", "_host"),
)

# also merge with phage RN data
# merged_1 = pd.merge(merged_1, phage_rn, left_on = "Subject_accession_ID_version", right_on = "sample", how = "left", suffixes = ("_host", "_phage"))

print(merged_1)

# print the first row of merged_1 with all columns
print(merged_1.iloc[0])

# write merged_1 to file
merged_1.to_csv(out_merged_all, index=False, sep="\t", header=True)


#### cOA analysis ####

list_of_cOA = ["ca3", "ca4", "ca6", "sam-amp"]
coa_results_dict = {}

template_coa_info = {
    "cOA": "-",  # name of cOA
    "total_interactions": 0,
    "total_positives": 0,
    "total_negatives": 0,
    "fraction_positives": 0,
    "fraction_negatives": 0,
}


for signal in list_of_cOA:
    print("Current signal: " + signal)
    # subset merged_1 to the with current signal True
    signal_true = merged_1[merged_1[signal] == True]

    print("Found " + str(len(signal_true)) + " interactions with " + signal + " True")

    # group all interactions from this cOA by phage
    phagecounts = (
        signal_true.groupby(["Subject_id"])["Subject_id"]
        .count()
        .reset_index(name="interactions")
    )  # number of interactions with this cOA

    # change to int
    phagecounts["interactions"] = phagecounts["interactions"].astype(int)

    # order by number of interactions
    phagecounts = phagecounts.sort_values(by="interactions", ascending=False)

    print("phagecounts: " + str(phagecounts))

    # rename Subject_id to phage
    phagecounts = phagecounts.rename(columns={"Subject_id": "phage"}, inplace=False)

    # merge with phage_rn to get more info on phages.
    phagecounts_fullInfo = pd.merge(
        phagecounts, phage_rn, left_on="phage", right_on="sample", how="left"
    )

    # write full info file
    phagecounts_fullInfo.to_csv(
        outfolder + signal + "_phagecounts.tsv", index=False, sep="\t", header=True
    )

    # create another dataframe where interactions are grouped by host
    hostcounts = (
        signal_true.groupby(["host"])["host"].count().reset_index(name="interactions")
    )  # number of interactions with this cOA

    # order by number of interactions
    hostcounts = hostcounts.sort_values(by="interactions", ascending=False)

    print("hostcounts: " + str(hostcounts))

    # write to file
    hostcounts.to_csv(
        outfolder + signal + "_hostcounts.tsv", index=False, sep="\t", header=True
    )

    # store to dict
    coa_results_dict[signal] = {"phagecounts": phagecounts, "hostcounts": hostcounts}

    # from the original merged_1, take all entries that are NOT the current signal
    signal_false = merged_1[merged_1[signal] == False]

    # group all interactions from this cOA by phage
    phagecounts_false = (
        signal_false.groupby(["Subject_id"])["Subject_id"]
        .count()
        .reset_index(name="interactions")
    )  # number of interactions with this cOA

    # change to int
    phagecounts_false["interactions"] = phagecounts_false["interactions"].astype(int)

    # order by number of interactions
    phagecounts_false = phagecounts_false.sort_values(
        by="interactions", ascending=False
    )

    print("phagecounts_false: " + str(phagecounts_false))

    # rename Subject_id to phage
    phagecounts_false = phagecounts_false.rename(
        columns={"Subject_id": "phage"}, inplace=False
    )

    # merge with phage_rn to get more info on phages.
    phagecounts_false_fullInfo = pd.merge(
        phagecounts_false, phage_rn, left_on="phage", right_on="sample", how="left"
    )

    # write full info file
    phagecounts_false_fullInfo.to_csv(
        outfolder + signal + "_phagecounts_false.tsv",
        index=False,
        sep="\t",
        header=True,
    )


# using coa_results_dict, create a mastertable in which each signal molecule is depicted by a column
# each row is a phage, and the columns contain the number of interactions with each signal molecule
phage_signal_spacer_mastertable = pd.DataFrame()
phage_signal_host_mastertable = pd.DataFrame()

counter = 0

for signal, results in coa_results_dict.items():
    if (
        counter == 0
    ):  # if on first iteration we need establish the mastertable from the first signal df
        phagecounts = results["phagecounts"]
        hostcounts = results["hostcounts"]

        phagecounts = phagecounts.rename(
            columns={"Subject_id": "phage", "interactions": signal + "_interactions"},
            inplace=False,
        )

        hostcounts = hostcounts.rename(
            columns={"host": "host", "interactions": signal + "_interactions"},
            inplace=False,
        )

        phage_signal_spacer_mastertable = phagecounts
        phage_signal_host_mastertable = hostcounts
        counter += 1
        # break current iteration of for loop
        continue
    else:
        phagecounts = results["phagecounts"]
        hostcounts = results["hostcounts"]

    # phagecounts is a dataframe with two columns: Subject_id and interactions
    # Subject_id is the phage name
    # interactions is the number of interactions with this cOA
    phagecounts = phagecounts.rename(
        columns={"Subject_id": "phage", "interactions": signal + "_interactions"},
        inplace=False,
    )

    hostcounts = hostcounts.rename(
        columns={"host": "host", "interactions": signal + "_interactions"},
        inplace=False,
    )

    phage_signal_spacer_mastertable = pd.merge(
        phage_signal_spacer_mastertable, phagecounts, on="phage", how="outer"
    )
    phage_signal_host_mastertable = pd.merge(
        phage_signal_host_mastertable, hostcounts, on="host", how="outer"
    )

# create column total interactions
phage_signal_spacer_mastertable["total_interactions"] = (
    phage_signal_spacer_mastertable.iloc[:, 1:].sum(axis=1)
)

# change to int
phage_signal_spacer_mastertable["total_interactions"] = phage_signal_spacer_mastertable[
    "total_interactions"
].astype(int)

# ca3_fraction
phage_signal_spacer_mastertable["ca3_interactions_fraction"] = round(
    phage_signal_spacer_mastertable["ca3_interactions"]
    / phage_signal_spacer_mastertable["total_interactions"],
    2,
)

# ca4_fraction
phage_signal_spacer_mastertable["ca4_interactions_fraction"] = round(
    phage_signal_spacer_mastertable["ca4_interactions"]
    / phage_signal_spacer_mastertable["total_interactions"],
    2,
)

# ca6_fraction
phage_signal_spacer_mastertable["ca6_interactions_fraction"] = round(
    phage_signal_spacer_mastertable["ca6_interactions"]
    / phage_signal_spacer_mastertable["total_interactions"],
    2,
)

# sam-amp_fraction
phage_signal_spacer_mastertable["sam-amp_interactions_fraction"] = round(
    phage_signal_spacer_mastertable["sam-amp_interactions"]
    / phage_signal_spacer_mastertable["total_interactions"],
    2,
)

# merge with phage_rn to get more info on phages.
phage_signal_spacer_mastertable = pd.merge(
    phage_signal_spacer_mastertable,
    phage_rn,
    left_on="phage",
    right_on="sample",
    how="outer",
)

# drop column sample
phage_signal_spacer_mastertable = phage_signal_spacer_mastertable.drop(
    columns=["sample", "Unnamed: 0"]
)

# replace empty cells with 0
phage_signal_spacer_mastertable = phage_signal_spacer_mastertable.fillna(0)

# reorder by total_interactions
phage_signal_spacer_mastertable = phage_signal_spacer_mastertable.sort_values(
    by="total_interactions", ascending=False
)

# save
phage_signal_spacer_mastertable.to_csv(
    outfolder + "phage_signal_spacer_mastertable.tsv",
    index=False,
    sep="\t",
    header=True,
)

# exit()


# For Gephi visualisation, we need the interaction data (edges) and bacterium/phage data (nodes).
# Nodes contain metadata on each node. These can be used in Gephi for filtering or aggregating things

# Gephi nodes. A node is either a bacterium or a phage.
print("Creating data for Gephi")
gephi_nodes_phages_temp = pd.DataFrame()
gephi_nodes_bacteria_temp = pd.DataFrame()

# phage id is the phage name, obtained by copying the Subject_accession_ID_version column from merged_1
gephi_nodes_phages_temp["Id"] = merged_1["Subject_accession_ID_version"]

# bacterium id is the host name, obtained by copying the host column from merged_1
gephi_nodes_bacteria_temp["Id"] = merged_1["host"]

# dropping duplicates from both dataframes
gephi_nodes_phages_temp = gephi_nodes_phages_temp.drop_duplicates()
gephi_nodes_bacteria_temp = gephi_nodes_bacteria_temp.drop_duplicates()

# set type for each node
gephi_nodes_phages_temp["type"] = "phage"
gephi_nodes_bacteria_temp["type"] = "host"

# creating gephi_nodes dataframe by "merging" phage and bacteria nodes. Resulting df contains both concatenated
gephi_nodes = pd.merge(
    gephi_nodes_phages_temp, gephi_nodes_bacteria_temp, on=["type", "Id"], how="outer"
)

# merge with host mastertable to get more info on hosts
gephi_nodes = pd.merge(
    gephi_nodes, host_mastertable, left_on="Id", right_on="sample", how="left"
)  # add info to nodes that are hosts

# merge with phage RN data
gephi_nodes = pd.merge(
    gephi_nodes,
    phage_rn,
    left_on="Id",
    right_on="sample",
    how="left",
    suffixes=("_host", "_phage"),
)

pprint("Gephi nodes:")
pprint(gephi_nodes)

CRISPRCas_types = ["I-", "II-", "III-", "IV-", "V-", "VI-"]
CRISPRCas_types_dict = {
    "I": ["I-A", "I-B", "I-B2", "I-C", "I-D", "I-E", "I-F1", "I-F2", "I-F3", "I-G"],
    "II": ["II-A", "II-B", "II-C"],
    "III": ["III-A", "III-B", "III-C", "III-D", "III-E", "III-F", "III-G", "III-H"],
    "IV": ["IV-A", "IV-B", "IV-C", "IV-D", "IV-E"],
    "V": ["V", "V-K"],
    "VI": ["VI"],
    "other": ["other"],
}

gephi_nodes = pd.merge(
    gephi_nodes, phages, left_on="Id", right_on="Accession", how="left"
)
# gephi_nodes = pd.merge(gephi_nodes, phagecounts, left_on = "Id", right_on="phage", how="left")


# Create new boolean columns for each CRISPR-Cas subtype.
# If any of the subtypes are present in the Subtype column, add True to corresponding type column
gephi_nodes["Subtype"] = gephi_nodes["Subtype"].fillna("")
for type, subtypes in CRISPRCas_types_dict.items():
    gephi_nodes[type] = gephi_nodes.Subtype.apply(set) != (
        gephi_nodes.Subtype.apply(set) - set(subtypes)
    )


gephi_nodes.to_csv(out_nodes, index=False, sep="\t", header=True)

# Gephi edges require Source and Target columns. These represent the arrows in Gephi. Each source and target must
# have equivalent values in the nodes file
gephi_edges = merged_1
gephi_edges = gephi_edges.rename(
    columns={"host": "Source", "Subject_accession_ID_version": "Target"}, inplace=False
)

gephi_edges.to_csv(
    outfolder + "gephi_edges_metadata.tsv", index=False, sep="\t", header=True
)

# remove all other columns except Source and Target
gephi_edges = gephi_edges[["Source", "Target"]]

gephi_edges.to_csv(out_edges, index=False, sep="\t", header=True)
