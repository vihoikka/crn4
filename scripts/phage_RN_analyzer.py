import argparse
import os
import pandas as pd
from Bio import SeqIO
import gffutils
import re
import sys
import subprocess
from tm_checker import run_TMHMM

"""
Test data:
--cas_operons_file /media/volume/st_andrews/cas10_corA_2/071_cctyper_loci/GCF_028335125.1_2/cas_operons.tsv
--locus GCF_028335125.1_2
--output_folder /media/volume/st_andrews/cas10_corA_2/test_out
--host_genomes_folder /media/volume/st_andrews/cas10_corA_2/06_host_genomes
--mode post_cctyper

Test command:
python3 catyper_prepper.py -i /media/volume/st_andrews/cas10_corA_2/071_cctyper_loci/GCF_028335125.1_2/cas_operons.tsv -l GCF_028335125.1_2 -o /media/volume/st_andrews/cas10_corA_2/test_out -hg /media/volume/st_andrews/cas10_corA_2/06_host_genomes
"""

print("Starting catyper_prepper.py")
# create args inputs for cas_operons_file, locus, output_folder, host_genomes_folder, mode, hmm_rows, hmm_targets, catyper_out
parser = argparse.ArgumentParser(
    description="Preps file for cATyping in the locus-centered Cas10/CorA pipeline"
)

parser.add_argument("-s", "--sample", help="sample", required=True)
parser.add_argument("-o", "--output_folder", help="output folder", required=True)
parser.add_argument("-pf", "--proteins_fasta", help="proteins as faa", required=True)
parser.add_argument("-r", "--hmm_rows", help="hmm rows", required=False)
parser.add_argument(
    "-t", "--hmm_targets", help="hmm targets", required=False
)  # this contains the effectors in the hmm database
parser.add_argument("-c", "--catyper_out", help="catyper output", required=False)
parser.add_argument("-ct", "--catyper_type", help="catyper type", required=False)
parser.add_argument(
    "-ep", "--effector_plot_data", help="effector plot data", required=False
)
parser.add_argument(
    "-rn", "--ring_nuclease", help="ring nuclease True of False", required=False
)
parser.add_argument(
    "-tm", "--tmhmm_model_path", help="tmhmm model path", required=False
)

args = parser.parse_args()

# generate bash script to run this using all arguments
# python

# assign args to similarly named variables
sample = args.sample
output_folder = args.output_folder
proteins_fasta = args.proteins_fasta
hmm_rows = args.hmm_rows
hmm_targets = args.hmm_targets
catyper_out = args.catyper_out
catyper_type = args.catyper_type
effector_plot_data = args.effector_plot_data
ring_nuclease = args.ring_nuclease
tmhmm_model_path = args.tmhmm_model_path

# when looking for effectors in the locus, the effector search range is the number of bases up or downstream of the cctyper defined cas operon boundaries
effector_search_range = 4000

plottable_effector_table = pd.DataFrame(
    columns=[
        "protein_id",
        "start",
        "end",
        "effector",
        "locus",
        "sample",
        "strand",
        "sequence",
    ]
)

# construct a precise effector dictionary using the original msa files that the hmms were derived from
print("Constructing precise effector dictionary")
hmm_target_file = open(hmm_targets, "r").read().splitlines()
effector_precise = {}

# here we generate effector names based on the hmm profile paths
for line in hmm_target_file:
    filename = os.path.basename(
        line
    )  # returns the filename from the path, e.g. ca3_nucc.faa or ca6_csm6-ca6_italicus#csm6-ca6.faa
    effector = re.split("\.", filename)[
        0
    ]  # removes extension. ca3_nucc.fa -> ca3_nucc or ca6_csm6-ca6_italicus#csm6-ca6.faa -> ca6_csm6-ca6_italicus#csm6-ca6
    if "#" not in effector:
        if ring_nuclease == "True":
            effector = re.split("_", effector)[
                0
            ]  # ring nucleae currently lack the cOA type in their name (29.11.2020)
        else:
            effector = re.split("_", effector)[1]  # ca3_nucc -> nucc
    elif "#" in effector:  # ca6_csm6-ca6_italicus#csm6-ca6
        effector = re.split("#", effector)[1]  # csm6-ca6

    effector_precise[effector] = False
    print("Current effector: " + effector)

protein_to_effector = (
    {}
)  # stores proteins as keys and associated effectors list as values
effector_to_protein = (
    {}
)  # stores effectors as keys and associated proteins list as values

for key, value in effector_precise.items():
    effector_to_protein[key] = []

print("Length of effector_precise dictionary: " + str(len(effector_precise)))
print("Effector precise dictionary: " + str(effector_precise))

# construct dictionaries for different types of cOAs
cOA_dict_binary = {
    "ca3": False,
    "ca4": False,
    "ca5": False,
    "ca6": False,
    "sam-amp": False,
    "unk": False,
    "val": False,
    "mem": False,
    "rng": False,
}  # TODO MODIFY WHEN ADDING NEW EFFECTORS
effector_dict = {
    "ca3": 0,
    "ca4": 0,
    "ca5": 0,
    "ca6": 0,
    "sam-amp": 0,
    "unk": 0,
    "val": 0,
    "rng": 0,
    "mem": 0,
}  # TODO MODIFY WHEN ADDING NEW EFFECTORS

# check if hmm result exists
print("Checking if hmm result exists")
if os.stat(hmm_rows).st_size == 0:
    print("No cA info found. This is an interesting locus.")
    results = {}
    results.update(
        cOA_dict_binary
    )  # this contains the cOAs. Everything is false by default
    results.update(effector_precise)  # this contains the effector itself
    results = pd.DataFrame([results])
    results["sample"] = sample
    results["no_RN"] = True
    filename = catyper_out
    print("Writing file " + filename)
    results.to_csv(filename, index=True, sep="\t", header=True)
    sys.exit()  # exit the script

# read in the hmm results if the results file is not empty
print("Reading hmm results")
hmm = pd.read_csv(hmm_rows, delim_whitespace=True, header=0, index_col=False)

protein_sequences = SeqIO.to_dict(SeqIO.parse(proteins_fasta, "fasta"))

dict_of_protein_matches_to_effectors_evalues = {}

min_lengths = (
    {  # minimum lengths for known effectors  #TODO MODIFY WHEN ADDING NEW EFFECTORS
        "crn1": 50,
        "crn2": 50,
        "crn3": 50,
        "crn4a": 50,
        "crn4b": 50,
        "csx16": 50,
        "csx20": 50,
        "csx15": 50,
        "unk01": 50,
        "solosavedRN": 50,
    }
)

max_lengths = {  # maximum lengths for known effectors #TODO MODIFY WHEN ADDING NEW EFFECTORS
    "crn1": 350,  # this is to avoid cross-annotation with csx6, which is usually over 200 AA. Now trying 350, check results manually
    "crn2": 250,
    "crn3": 250,
    "crn4a": 250,
    "crn4b": 250,
    "csx16": 250,
    "csx20": 250,
    "csx15": 250,
    "unk01": 250,
    "solosavedRN": 250,
}

print("...Checking hmm results for cA type...")
results = {}
for index, row in hmm.iterrows():  # check each row in hmm result
    print(row)
    protein_id = row["query_name"]  # returns the protein id
    clade = "noCladeInfo"  # clade is used to differentiate between different clades of the same protein (some HMM profiles consist of multiple clades)
    # check if the hmm result contains the word 'icity'. This refers to CorA and is a remnant from the HMM profile fil

    if (
        "_" in row["target_name"]
    ):  # if the hmm result contains an underscore, it is one of the other effectors
        # if the hmm result contains a hashtag #, then the bit followed by the hashtag is precise_hit
        if "#" in row["target_name"]:
            # name with # are in the format ca6_csm6-ca6_italicus#csm6-ca6_clustered_aligned. Whenever a hastag is present, a clade that differentiates the effector from other clades is present (the last part of the name before the hashtag, in this case italicus)
            # the signal molecule is first fetched splitting by _ and getting first hit
            split_underscore = re.split("_", row["target_name"])
            hit_type = split_underscore[0]
            # the actual effector is fetched splitting by # and getting second hit, and splitting it by _ and getting first hit. For example, the name could be ca6_csm6-ca6_italicus#csm6-ca6
            split_hashtag = re.split(
                "#", row["target_name"]
            )  # split_hashtag becomes ['ca6_csm6-ca6_italicus', 'csm6-ca6_clustered_aligned']
            precise_hit = re.split("_", split_hashtag[1])[
                0
            ]  # precise_hit becomes csm6-ca6
            clade = split_hashtag[0].split("_")[-1]  # clade becomes italicus
            print("Precise hit: " + str(precise_hit) + " Clade: " + str(clade))

            if (
                ring_nuclease == "True"
            ):  # if we are running catyper_prepper10 in run_ring_nuclease mode
                print(
                    "...Ring nuclease detected with hashtag: " + str(row["target_name"])
                )
                print("....The precise hit is " + str(precise_hit))
                hit_type = "rng"
                # precise_hit = re.split("_", row["target_name"])[0]
                # ring nuclease clades are defined by the last _ part of the name prior to hashtag
                try:
                    clade = split_hashtag[0].split("_")[-1]
                except:
                    clade = "noCladeInfo"

        else:  # if no hashtag is present, the precise_hit is simply the effector or ring nuclease after the underscore
            if (
                ring_nuclease == "True"
            ):  # if we are running catyper_prepper10 in run_ring_nuclease mode
                print("...Ring nuclease detected")
                hit_type = "rng"
                precise_hit = re.split("_", row["target_name"])[0]
            else:
                hit = re.split("_", row["target_name"])
                hit_type = hit[0]  # returns ca3, ca4, ca5, ca6 or SAM-AMP
                precise_hit = hit[
                    1
                ]  # returns the hmm target, e.g. nucc, can1, can2, cora...

    print("...Discovered cA type in hmm: " + str(hit_type) + " (" + precise_hit + ")")

    target_score = row["score_fullseq"]  # the hmm score for the target sequence

    # check if the protein is in the cas operon and exceeds length cutoff for given effector. Note that the borders are already expanded prior to this step
    print("Checking length cutoffs and positional filters for " + precise_hit)
    print(
        "Comparing protein length "
        + str(len(protein_sequences[protein_id].seq))
        + " against cutoff min "
        + str(min_lengths[precise_hit])
        + " and max "
        + str(max_lengths[precise_hit])
    )
    protein_length = len(protein_sequences[protein_id].seq)
    if (
        protein_length > min_lengths[precise_hit]
        and protein_length < max_lengths[precise_hit]
    ):
        print(
            "...Protein "
            + protein_id
            + " is long/short enough at "
            + str(protein_length)
            + " AA"
        )
        #
        # if the protein is listed in the custom evalues, check if it has already been detected by another effector in that list
        # check if the protein_id exists as key in the protein_to_effector dictionary
        precise_hit_and_score = {
            "RN": precise_hit,
            "clade": clade,
            "score": target_score,
            "hit_type": hit_type,
            "sample": sample,
            "type": catyper_type,
            "sequence": protein_sequences[protein_id].seq,
        }
        if protein_id in protein_to_effector:
            print(
                "...Protein "
                + protein_id
                + " already exists in protein_to_effector dictionary: "
                + str(protein_to_effector[protein_id])
            )
            print(
                "Checking the score of the existing hit ("
                + str(
                    protein_to_effector[protein_id][0]["RN"]
                    + ") and comparing to new one: "
                    + str(protein_to_effector[protein_id][0]["score"])
                    + " vs new score of "
                    + str(target_score)
                    + " in "
                    + str(precise_hit)
                )
            )
            if (
                target_score > protein_to_effector[protein_id][0]["score"]
            ):  # if the current hit has a higher score than the existing hit, replace the existing hit with the current hit
                print(
                    "...Current hit has a higher score than the existing hit. Replacing existing hit with current hit"
                )
                protein_to_effector[protein_id] = (
                    []
                )  # initiate new key and empty list for current protein
                protein_to_effector[protein_id].append(
                    precise_hit_and_score
                )  # add the precise hit to the list

                effector_to_protein[precise_hit].append(
                    protein_id
                )  # add the protein to the list of proteins associated with the hmm

            else:
                print(
                    "...Current hit has a equal or lower score than the existing hit. Not adding new hit to protein_to_effector dictionary"
                )
        else:
            print(
                "...Protein "
                + protein_id
                + " is not in the protein_to_effector dictionary. Adding it now."
            )
            protein_to_effector[protein_id] = (
                []
            )  # initiate new key and empty list for current protein
            protein_to_effector[protein_id].append(
                precise_hit_and_score
            )  # add the precise hit and its hmm score to the list

        # use the protein id to get the protein sequence from the protein multifasta file
        protein_sequence = protein_sequences[protein_id].seq

    else:
        print("...Protein is not in the cas operon or is too short")

# after going through all the proteins, determine which hit_types are present in the locus
print("...Determining which hit_types (cNTs) etc are present in the locus")
for key, value in protein_to_effector.items():  # iterate through all proteins
    if len(value) > 0:  # if the protein has a hit in the hmm
        protein_id = key
        precise_hit = value[0][
            "RN"
        ]  # get the effector from the protein_to_effector dictionary
        sample = value[0]["sample"]
        sequence = value[0]["sequence"]
        hit_type = value[0]["hit_type"]  # get the hit_type

        cOA_dict_binary[hit_type] = (
            True  # set the hit_type to True in the cOA_dict_binary
        )
        effector_dict[hit_type] += 1  # add 1 to the effector_dict for the hit_type

        effector_to_protein[precise_hit].append(
            protein_id
        )  # add the protein to the list of proteins associated with the hmm

        # update the effector_precise
        effector_precise[precise_hit] = True

        protein_name_locus = sample + "__" + protein_id
        with open(os.path.join(output_folder, precise_hit + ".faa"), "w") as f:
            f.write(">" + protein_name_locus + "\n" + str(sequence) + "\n")

        # add the protein to plottable_effector_table
        plottable_effector_table = plottable_effector_table.append(
            {
                "protein_id": protein_id,
                "effector": precise_hit,
                "sample": sample,
                "type": hit_type,
                "sequence": protein_sequence,
            },
            ignore_index=True,
        )


results.update(cOA_dict_binary)
results.update(effector_precise)

results = pd.DataFrame([results])
results["sample"] = sample
results["no_effectors"] = False

print("Results:")
print(results)

filename = catyper_out
print("Writing file " + filename)
results.to_csv(filename, index=True, sep="\t", header=True)

# For protein_to_effector_df.
# Each protein will have its own row with the associated effector in the 'Effector' column.
protein_to_effector_df = pd.DataFrame(
    [(k, v_i) for k, v in protein_to_effector.items() for v_i in v],
    columns=["Protein", "Effector"],
)
print(protein_to_effector_df)

# add the current locus as a column to the df
protein_to_effector_df["Sample"] = sample

# For effector_to_protein_df.
# Each effector will have its own row for every protein it is associated with in the 'Protein' column
effector_to_protein_df = pd.DataFrame(
    [(k, v_i) for k, v in effector_to_protein.items() for v_i in v],
    columns=["Effector", "Protein"],
)
effector_to_protein_df["Sample"] = sample
print(effector_to_protein_df)

# save the effector_to_protein_df and protein_to_effector_df to the output folder
print(
    "Saving effector_to_protein_df to "
    + str(os.path.join(output_folder, sample + "_effector_to_protein.tsv"))
)
effector_to_protein_df.to_csv(
    os.path.join(output_folder, sample + "_effector_to_protein.tsv"),
    index=True,
    sep="\t",
    header=False,
)
protein_to_effector_df.to_csv(
    os.path.join(output_folder, sample + "_protein_to_effector.tsv"),
    index=True,
    sep="\t",
    header=False,
)

# write plottable_effector_table to file
print("Writing plottable_effector_table to file")
plottable_effector_table.to_csv(effector_plot_data, index=False, sep="\t", header=True)
