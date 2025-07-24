import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_folder", help="input folder")
parser.add_argument("-o", "--output_folder", help="output folder")
parser.add_argument("-in", "--interference_cutoff", help="interference cutoff")
args = parser.parse_args()

inputfolder = args.input_folder
outputfolder = args.output_folder
interference_cutoff = args.interference_cutoff

samples = [f for f in os.listdir(inputfolder) if os.path.isdir(os.path.join(inputfolder, f))]

def processInput(sample):
    if os.path.isfile(f"{inputfolder}/{sample}/cas_operons.tab"):
        cas_operons = f"{inputfolder}/{sample}/cas_operons.tab"
    else:
        print(f"No cas_operons file found for {sample}")
        return

    cas_operons_df = pd.read_csv(cas_operons, sep='\t', header = 0)
    cas_operons_df = cas_operons_df[cas_operons_df['Prediction'].str.contains("III-")]
    cas_operons_df = cas_operons_df[cas_operons_df['Genes'].str.contains("cas10", case=False)]
    cas_operons_df = cas_operons_df[~cas_operons_df['Genes'].str.contains("cas10.*cas10", case=False)]
    print(f"Discarded {len(cas_operons_df)} rows where cas10 was not found exactly once")

    print(f"No of complete cas_operons found: {str(len(cas_operons_df))}")

    cas_operons_df["locus_id"] = sample + "_" + cas_operons_df.index.astype(str)
    cas_operons_df['sample'] = sample
    
    for index, row in cas_operons_df.iterrows():
        print(row["locus_id"])
        os.makedirs(f"{outputfolder}/{row['locus_id']}")
        row_transposed = row.to_frame().T
        row_transposed.to_csv(f"{outputfolder}/{row['locus_id']}/cas_operons.tsv", sep='\t', header=True, index=True)

print("Processing input...")
for sample in samples:
    processInput(sample)
