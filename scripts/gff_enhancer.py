from BCBio import GFF
from dna_features_viewer import GraphicFeature, GraphicRecord, BiopythonTranslator
import matplotlib.pyplot as plt
import argparse
import pandas as pd
import os
import gffutils

# create args inputs for cas_operons_file, locus, output_folder, host_genomes_folder, mode, hmm_rows, hmm_targets, catyper_out, cctyper_protein_table
parser = argparse.ArgumentParser(description="Visualises expanded CRISPR loci")
parser.add_argument("-i", "--cas_operons_file", help="input file", required=True)
parser.add_argument("-c", "--cctyper_folder", help="cctyper folder", required=True)
parser.add_argument("-l", "--locus", help="locus", required=True)
parser.add_argument("-o", "--output_folder", help="output folder", required=True)
parser.add_argument(
    "-hg", "--host_genomes_folder", help="host genomes folder", required=True
)
parser.add_argument(
    "-v", "--validated_effectors", help="validated effectors", required=True
)
parser.add_argument(
    "-cc", "--cctyper_protein_table", help="cctyper plottable table", required=True
)
parser.add_argument(
    "-k", "--known_effector_table", help="known effectors", required=True
)
parser.add_argument("-rn", "--ring_nucleases", help="ring nucleases", required=True)

args = parser.parse_args()

# generate bash command for the above
# python locus_visualiser.py -i cas_operons.tsv -l locus -o output_folder -hg host_genomes_folder

# assign args to similarly named variables
cas_operons_file = args.cas_operons_file
locus = args.locus
output_folder = args.output_folder
host_genomes_folder = args.host_genomes_folder
cctyper_folder = args.cctyper_folder
validated_effectors = args.validated_effectors
cctyper_protein_table = args.cctyper_protein_table
known_effector_table = args.known_effector_table
ring_nucleases_table = args.ring_nucleases

validated_effectors = pd.read_csv(validated_effectors, sep="\t", header=0)
cctyper_protein_table = pd.read_csv(cctyper_protein_table, sep="\t", header=0)
known_effector_table = pd.read_csv(known_effector_table, sep="\t", header=0)
ring_nucleases_table = pd.read_csv(ring_nucleases_table, sep="\t", header=0)

# from cctyper_protein_table, remove columns start	end	annotation_gff	annotation_cctyper	strand	evalue	sequence
cctyper_protein_table = cctyper_protein_table.drop(
    columns=[
        "start",
        "end",
        "annotation_gff",
        "annotation_cctyper",
        "strand",
        "evalue",
        "sequence",
    ]
)

# merge all three tables by id
merged = validated_effectors.merge(
    cctyper_protein_table,
    on="protein_id",
    how="outer",
    suffixes=("_validated", "_cctyper"),
)
merged = merged.merge(
    known_effector_table, on="protein_id", how="outer", suffixes=("", "_known")
)
merged = merged.merge(
    ring_nucleases_table, on="protein_id", how="outer", suffixes=("", "_ring")
)

# when looking for effectors in the locus, the effector search range is the number of bases up or downstream of the cctyper defined cas operon boundaries
effector_search_range = 20000

# get sample name by splitting the locus name at the underscore and taking the first two parts
sample = locus.split("_")[0] + "_" + locus.split("_")[1]

# # read in the cas operons file
cas_operons_df = pd.read_csv(cas_operons_file, sep="\t", header=0)

# # get contig ID of the CRISPR positive contig
contig = cas_operons_df.iloc[0]["Contig"]  # This is wrong.

# # define start and end coordinates for plotting
cas_operon_start = int(cas_operons_df["Start"][0]) - effector_search_range
cas_operon_end = int(cas_operons_df["End"][0]) + effector_search_range
# features_list = []


def enhance_gff_with_annotations(
    gff_db,
    annotation_table,
    start_col,
    end_col,
    annotation_col,
    new_gff_file,
    cas_operon_start,
    cas_operon_end,
):
    """
    Enhance a GFF file with new annotations from an annotation table.

    Parameters:
    - gff_db: The gffutils database object.
    - annotation_table: The DataFrame containing the annotations.
    - start_col: The column name in the annotation table that corresponds to the start position.
    - end_col: The column name in the annotation table that corresponds to the end position.
    - annotation_col: The column name in the annotation table that contains the annotations.
    - new_gff_file: The path to the new GFF file to write enhanced annotations to.
    - cas_operon_start: The start position of the cas operon.
    - cas_operon_end: The end position of the cas operon.
    """
    if not annotation_table.empty:
        with open(new_gff_file, "a") as new_gff:  # Use "a" to append to the file
            for feature in gff_db.features_of_type("CDS"):
                # Check if the feature is within the cas operon boundaries
                if feature.start >= cas_operon_start and feature.end <= cas_operon_end:
                    # Get the annotation from the annotation table
                    annotation = annotation_table.loc[
                        (annotation_table[start_col] == feature.start)
                        & (annotation_table[end_col] == feature.end),
                        annotation_col,
                    ].values
                    annotation = annotation[0] if len(annotation) > 0 else ""

                    # get the feature's annotation attribute if it exists, otherwise make it Other
                    if "annotation" in feature.attributes:
                        current_gff_annotation = feature.attributes["annotation"]
                    else:
                        feature.attributes["annotation"] = "Other"
                        current_gff_annotation = "Other"

                    if (len(annotation) > 0) & (current_gff_annotation == "Other"):
                        if ("Cas" in annotation) or ("Csm" in annotation):
                            annotation = annotation.split("_")[0]
                            print("Adding annotation to GFF: ", annotation)
                        annotation = annotation
                    else:
                        annotation = "Other"
                        print("Adding Other annotation to GFF: ")

                    # Add the annotation to the feature
                    feature.attributes["annotation"] = annotation
                    new_gff.write(str(feature) + "\n")


# Paths and input files
gff_file = os.path.join(host_genomes_folder, sample, sample + "_features.gff")
new_gff_file = os.path.join(output_folder, locus + "_enhanced.gff")

# Create GFF database in memory
gff_db = gffutils.create_db(
    gff_file,
    dbfn=":memory:",
    force=True,
    keep_order=True,
    merge_strategy="merge",
    sort_attribute_values=True,
)

# Print first rows of data for debugging
# input("Press Enter to continue...")

# Define cas operon boundaries (these should be defined elsewhere in your script)
# cas_operon_start = ...
# cas_operon_end = ...

# Enhance GFF with annotations from various tables
enhance_gff_with_annotations(
    gff_db=gff_db,
    annotation_table=cctyper_protein_table,
    start_col="start_cctyper",
    end_col="end_cctyper",
    annotation_col="annotation_cctyper_cctyper",
    new_gff_file=new_gff_file,
    cas_operon_start=cas_operon_start,
    cas_operon_end=cas_operon_end,
)

enhance_gff_with_annotations(
    gff_db=gff_db,
    annotation_table=known_effector_table,
    start_col="start",
    end_col="end",
    annotation_col="effector",
    new_gff_file=new_gff_file,
    cas_operon_start=cas_operon_start,
    cas_operon_end=cas_operon_end,
)

enhance_gff_with_annotations(
    gff_db=gff_db,
    annotation_table=ring_nucleases_table,
    start_col="start",
    end_col="end",
    annotation_col="effector",
    new_gff_file=new_gff_file,
    cas_operon_start=cas_operon_start,
    cas_operon_end=cas_operon_end,
)

enhance_gff_with_annotations(
    gff_db=gff_db,
    annotation_table=validated_effectors,
    start_col="start",
    end_col="end",
    annotation_col="effector",
    new_gff_file=new_gff_file,
    cas_operon_start=cas_operon_start,
    cas_operon_end=cas_operon_end,
)
