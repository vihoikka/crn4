#    Takes protein_to_effector data from all three search modules (known effectors, new validated effectors and ring nucleases)
#    and by comparing the proteins listed in ring nucleases, creates a table
#    showing to which effectors the ring nucleases are fused to


#create args for         python scripts/ring_fusions.py --known_effectors {input.known_effectors} --validated_new_effectors {input.validated_new_effectors} --ring_nucleases {input.ring_nucleases} --output {output.ring_fusions}

import argparse
import os
import subprocess
import sys
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import ast
import json
import re

# Parse command-line arguments for input, output and msa and faa. Also use_existing_alignment and existing_alignment_path
parser = argparse.ArgumentParser(
    description='Determines fusions between ring nucleases and effectors.')
parser.add_argument('-k', '--known_effectors', type=str, metavar='known_effectors', required=True)
parser.add_argument('-v', '--validated_new_effectors', type=str, metavar='validated_new_effectors', required=True)
parser.add_argument('-r', '--ring_nucleases', type=str, metavar='ring_nucleases', required=True)
parser.add_argument('-o', '--output', type=str, metavar='output', required=True)

args = parser.parse_args()

known_effectors = args.known_effectors
validated_new_effectors = args.validated_new_effectors
ring_nucleases = args.ring_nucleases
output = args.output

def dict_converter(data): # a custom converter for json/dict for pandas
    dict_string = re.sub(r"Seq\((.*?)\)", r"\1", data) #removing the Seq object from the string because it can't be handled by as.literal_eval
    dict = ast.literal_eval(dict_string)
    return dict

# Read in the tables. There are no headers and the second column should be annotated as a dictionary
known_effectors_table = pd.read_csv(known_effectors, sep='\t', header=None, converters={2:dict_converter})
validated_new_effectors_table = pd.read_csv(validated_new_effectors, sep='\t', header=None, converters={2:dict_converter})
ring_nucleases_table = pd.read_csv(ring_nucleases, sep='\t', header=None, converters={2:dict_converter})



#for each, add columns: something, protein_id, effector_dict, locus
columns = ["something", "protein_id", "effector_dict", "locus"]
known_effectors_table.columns = columns
validated_new_effectors_table.columns = columns
ring_nucleases_table.columns = columns


#add column clade to ring_nucleases_table
ring_nucleases_table["rn_clade"] = ring_nucleases_table["effector_dict"].apply(lambda x: x["clade"])


#For each create a new column called "effector" and set it to the value of the dictionary key "effector" from the column "effector_dict"
known_effectors_table["effector"] = known_effectors_table["effector_dict"].apply(lambda x: x["effector"])
ring_nucleases_table["effector"] = ring_nucleases_table["effector_dict"].apply(lambda x: x["effector"])
validated_new_effectors_table["effector"] = validated_new_effectors_table["effector_dict"].apply(lambda x: x["effector"])

#rename the column "effector" in ring_nucleases_table to "ring_nuclease"
ring_nucleases_table = ring_nucleases_table.rename(columns={"effector": "ring_nuclease"})

print("Validated new effectors:")
print(validated_new_effectors_table)
print("Number of validated new effectors: ", len(validated_new_effectors_table))

print("Known effectors:")
print(validated_new_effectors_table)
print("Number of known effectors: ", len(known_effectors_table))

print("Ring nucleases:")
print(ring_nucleases_table)

#merge ring nucleases and known effectors based on protein_id. From known effectors only bring the effector column
merged = pd.merge(ring_nucleases_table, known_effectors_table[["protein_id", "effector"]], on="protein_id", how="left")
#rename effector column known_effector
merged = merged.rename(columns={"effector": "known_effector"})

#merge also with validated new effectors, similar to above
merged = pd.merge(merged, validated_new_effectors_table[["protein_id", "effector"]], on="protein_id", how="left")
#rename effector column validated_new_effector
merged = merged.rename(columns={"effector": "validated_new_effector"})

print("Merged:")
print(merged)
print("Number of merged: ", len(merged))

#create new column effector which takes value from either known_effector or validated_new_effector, depending on which is not null
merged["effector"] = merged["known_effector"].fillna(merged["validated_new_effector"])


#create column "fusion_protein" and set to True for all
merged["fusion_protein"] = True

#create column fusion_components, which is just effector + _ + ring_nuclease
merged["fusion_components"] = merged["effector"] + "_" + merged["ring_nuclease"]

#drop effector_dict
merged = merged.drop(columns=["effector_dict", "something"])

print("Merged:")
print(merged)

#remove rows where fusion_protein is nul or an empty string
merged = merged[merged["fusion_components"].notnull()]
merged = merged[merged["fusion_components"] != ""]

print("Merged after removing empty fusion_components:")
print(merged)

#write to file
merged.to_csv(output, sep='\t', index=False)