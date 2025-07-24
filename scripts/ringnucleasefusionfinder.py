import argparse
import os
import pandas as pd
from Bio import SeqIO
import gffutils
import re
import sys

'''
Finds and reports multiple hits in hmm results for ring nuclease domains
'''

parser = argparse.ArgumentParser(description='Preps file for cATyping in the locus-centered Cas10/CorA pipeline')
parser.add_argument('-l', '--locus', help='locus', required=True)
parser.add_argument('-o', '--output_folder', help='output folder', required=True)
parser.add_argument('-r', '--hmm_rows', help='hmm rows', required=False)

args = parser.parse_args()

#assign args to similarly named variables
locus = args.locus
output_folder = args.output_folder
hmm_rows = args.hmm_rows

#open hmmrows in pandas
ring_nuclease_hmm_results = pd.read_csv(hmm_rows, delim_whitespace=True, header = 0, index_col=False)
print(ring_nuclease_hmm_results)

#divide the dataframe into multiple dataframes based on the query_name
ring_nuclease_hmm_results = ring_nuclease_hmm_results.groupby("query_name")
print("Grouped by query_name") 
print(ring_nuclease_hmm_results)

#create new dataframe with columns locus, ringToRing_fusion
ring_nuclease_fusion_df = pd.DataFrame(columns=["locus", "ringToRing_fusion", "protein"])

#iterate through the dataframes
for name, group in ring_nuclease_hmm_results:
    print("Iterating...")
    print("Group: " + str(group))
    print("Name: " + str(name))
    print(group.head())
    #print all values in the column target_name
    print(group["target_name"])
    #remove duplicates rows based on the column target_name
    group = group.drop_duplicates(subset="target_name")
    #if the dataframe contains more than 1 row
    if (len(group) > 1):
        print("Multiple hits found in hmm results for " + str(name))
        print(group)
        print("The targets in group: " + str(group["target_name"]))
        group["target_name"] = group["target_name"].str.split("_").str[0]
        concat = group["target_name"].str.cat(sep="_")
        target_protein = group["query_name"]
        #add the concatenated string to ring_nuclease_fusion_df
        ring_nuclease_fusion_df = ring_nuclease_fusion_df.append({"locus": locus, "ringToRing_fusion": concat, "protein": target_protein}, ignore_index=True)

ring_nuclease_fusion_df.to_csv(os.path.join(output_folder, locus + "_ring_nuclease_fusions.tsv"), index=False, sep = '\t', header = True)