"""
This file is the python script for a snakemake rule. 
The rule brings together all_proteins_against_RN_diamond, phigaro, and contig_to_genome_mapping.
contig_to_genome_mapping is used to map the accession numbers of the diamond blasts to the genomic accessions.
Then, the phigaro results are used to find prophage boundaries and mark whether each hit is within
a prophage region or now. We also use CRISPR-Cas locus coordinates from rule.

Input files are (from snakemake):
    RN_diamond = rules.concatenate_all_proteins_against_RN_diamond.output,
    phigaro = rules.concatenate_phigaro.output,
    contig_to_genome_mapping = rules.concatenate_contig_to_genome_mapping.output,
    cas_operons = rules.concatenate_crispr_locus_coordinates.output

Example run
python scripts/RN_hits_outside_crispr_merger.py --RN_diamond {input.RN_diamond} --phigaro {input.phigaro} --contig_to_genome_mapping {input.contig_to_genome_mapping} --cas_operons {input.cas_operons} --out {output}
"""

import pandas as pd
import argparse


def main():
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--RN_diamond", help="Input RN diamond blast results")
    argparser.add_argument("--phigaro", help="Input phigaro results")
    argparser.add_argument(
        "--contig_to_genome_mapping", help="Input contig to genome mapping file"
    )
    argparser.add_argument("--cas_operons", help="Input CRISPR-Cas locus coordinates")
    argparser.add_argument("--out", help="Output file")

    args = argparser.parse_args()

    RN_diamond = args.RN_diamond
    phigaro = args.phigaro
    contig_to_genome_mapping = args.contig_to_genome_mapping
    cas_operons = args.cas_operons
    out = args.out

    # read in RN diamond blast results
    RN_diamond_df = pd.read_csv(RN_diamond, sep="\t", header=0)

    # the diamond_target columns is a string with RN|<protein_id>. We want to save the RN to a new column
    RN_diamond_df["diamond_target"] = (
        RN_diamond_df["diamond_target"].str.split("|").str[0]
    )

    # read in phigaro results
    phigaro_df = pd.read_csv(phigaro, sep="\t")

    # read in contig to genome mapping
    contig_to_genome_mapping_df = pd.read_csv(contig_to_genome_mapping, sep="\t")

    # read in CRISPR-Cas locus coordinates
    cas_operons_df = pd.read_csv(cas_operons, sep="\t")

    general_evalue = 1e-40

    evalue_rn_dict = {
        "crn1": general_evalue,
        "crn2": general_evalue,
        "crn3": general_evalue,
        "crn4a": general_evalue,
        "crn4b": general_evalue,
        "csx15": general_evalue,
        "csx16": general_evalue,
        "csx20": general_evalue,
        "csx1": general_evalue,
    }

    # Function to get the threshold for a given diamond_target
    def get_threshold(target):
        return evalue_rn_dict.get(target, general_evalue)

    # Apply the threshold filter
    RN_diamond_df = RN_diamond_df[
        RN_diamond_df.apply(
            lambda row: row["evalue"] < get_threshold(row["diamond_target"]), axis=1
        )
    ]

    print(RN_diamond_df)
    print(phigaro_df)
    print(contig_to_genome_mapping_df)

    # add prophage_ prefix to all columns in phigaro_df
    phigaro_df.columns = ["prophage_" + col for col in phigaro_df.columns]

    # merge contig/host maps with prophage results
    prophage_coordinates = pd.merge(
        phigaro_df,
        contig_to_genome_mapping_df,
        left_on="prophage_scaffold",
        right_on="contig",
        how="left",
    )

    prophage_coordinates.to_csv("prophage_coordinates.tsv", sep="\t", index=False)

    # merge with RN diamond blast results
    merged = pd.merge(
        RN_diamond_df, prophage_coordinates, left_on="host", right_on="host", how="left"
    )

    merged.to_csv("merged_rn_prophage.tsv", sep="\t", index=False)

    print(merged)

    # from cas_operons_df, print the row where host GCF_000175295.2
    print(cas_operons_df[cas_operons_df["sample"] == "GCF_000175295.2"])

    # merge with CRISPR-Cas locus coordinates. Add a prefix to the columns in cas_operons_df
    merged = pd.merge(
        merged, cas_operons_df, left_on="host", right_on="sample", how="left"
    )

    merged.to_csv("merged_crispr_coords.tsv", sep="\t", index=False)

    """
    The next part examines each ring nuclease to see if its coordinates are within a CRISPR-Cas locus or within a prophage region.
    Note that a ring nuclease may be present multiple times in the table, if the host contains multiple crispr loci or multiple prophages.
    """

    # create a new column to indicate whether the ring nuclease is within a CRISPR-Cas locus
    merged["within_crispr"] = False
    merged["within_prophage"] = False
    merged["close_crispr"] = False

    crispr_closensess_cutoff = 1000  # the number of base pairs that a ring nuclease can be from a CRISPR-Cas locus and still be considered "close" to it

    for index, row in merged.iterrows():  #
        if row["host"] == "GCF_000175295.2":
            # we want to track this specific sample. Print details about hit:
            print("Hit details:")
            print(row)
            print("CRISPR coordinates: ", row["Start"], row["End"])
        if row["diamond_contig"] == row["Contig"] and (
            row["genomic_start"] >= row["Start"] and row["genomic_end"] <= row["End"]
        ):  # if the ring nuclease coordinates are within the CRISPR-Cas locus coordinates and on the same contig
            if row["host"] == "GCF_000175295.2":
                print("Ring nuclease is within CRISPR-Cas locus")
            merged.at[index, "within_crispr"] = True
        # similar if statement but check if "close" to crispr locus
        elif (
            row["genomic_start"] >= row["Start"] - crispr_closensess_cutoff
            and row["genomic_end"] <= row["End"] + crispr_closensess_cutoff
        ):
            merged.at[index, "close_crispr"] = True
        if (
            row["diamond_contig"] == row["prophage_scaffold"]
            and row["genomic_start"] >= row["prophage_begin"]
            and row["genomic_end"] <= row["prophage_end"]
        ):
            merged.at[index, "within_prophage"] = True

    merged.to_csv(out, sep="\t", index=False)

    # Next, create a summary table that shows the number of ring nucleases within a CRISPR-Cas locus and within a prophage region
    summary = (
        merged.groupby(["host", "within_crispr", "within_prophage"])
        .size()
        .reset_index(name="counts")
    )

    summary.to_csv("rn_summary.tsv", sep="\t", index=False)

    # Next, collapse rows so that each protein exists only once in the table. Currently the reason we have multiple rows for each protein is because of the multiple prophages and CRISPR-Cas loci in the host. These need to be collapsed into a single row.
    # For each protein, we will create a new column that lists the prophages and CRISPR-Cas loci that the protein is within.
    # We will also create a new column that lists the prophages and CRISPR-Cas loci that the protein is not within.
    # Finally, we will create a new column that lists the prophages and CRISPR-Cas loci that the protein is within, but not within a CRISPR-Cas locus.

    # First, create a new column that lists the prophages and CRISPR-Cas loci that the protein is within.
    merged["within_prophage_list"] = ""
    merged["within_crispr_list"] = ""
    merged["close_crispr_list"] = ""
    merged["within_prophage_not_crispr_list"] = ""

    for index, row in merged.iterrows():
        if row["within_prophage"]:
            merged.at[index, "within_prophage_list"] = row["prophage_scaffold"]
        if row["within_crispr"]:
            merged.at[index, "within_crispr_list"] = row["host"]
        if row["within_prophage"] and not row["within_crispr"]:
            merged.at[index, "within_prophage_not_crispr_list"] = row[
                "prophage_scaffold"
            ]
        if row["close_crispr"]:
            merged.at[index, "close_crispr_list"] = row["host"]

    # Next, collapse rows so that each protein exists only once in the table
    collapsed = (
        merged.groupby("diamond_query")
        .agg(
            {
                "diamond_target": "first",
                "pident": "first",
                "diamond_length": "first",
                "mismatch": "first",
                "gapopen": "first",
                "qstart": "first",
                "qend": "first",
                "sstart": "first",
                "send": "first",
                "evalue": "first",
                "diamond_bitscore": "first",
                "genomic_start": "first",
                "genomic_end": "first",
                "host": "first",
                "within_crispr": "first",
                "within_prophage": "first",
                "within_prophage_list": "first",
                "within_crispr_list": "first",
                "within_prophage_not_crispr_list": "first",
                "close_crispr_list": "first",
            }
        )
        .reset_index()
    )

    collapsed.to_csv("collapsed_rn.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
