'''
    The snakemake rule that calls this script is shown below:

    Uses the locus list from cOA_RN_explorer and the "unknown" protein from tule unknown_finder.
    It extracts all <150 AA unknown proteins from the loci defined by cOA_RN_explorer
    and writes them to a new .tsv file and .faa file.
    input:
        cA3_loci = rules.cOA_RN_explorer.output.ca3_loci,
        all_unknown_proteins = rules.concatenate_unknowns.output.info
    output:
        cA3_small_unk_proteins_info = base_path + "/92_cA3_small_proteins/cA3_small_proteins.tsv",
        cA3_small_unk_proteins_faa = base_path + "/92_cA3_small_proteins/cA3_small_proteins.faa"
    params:
        outfolder = base_path + "/92_cA3_small_proteins"
    threads: thread_ultrasmall
    shell:
        python scripts/cA3_small_unk_proteins.py --loci "{input.cA3_loci},{input.cA4_loci},{input.cA6_loci},{input.sam_amp_loci}" --unknown_proteins {input.all_unknown_proteins} --output_folder {params.outfolder} --length_cutoff {params.length_cutoff}
'''

import argparse
import pandas as pd
import os
import Bio.SeqIO

argparser = argparse.ArgumentParser()
argparser.add_argument('--loci', help='loci file', required=True)
argparser.add_argument('--unknown_proteins', help='unknown proteins file', required=True)
argparser.add_argument('--output_folder', help='output folder', required=True)
argparser.add_argument('--length_cutoff', help='length cutoff for small proteins', required=True)
argparser.add_argument('--all_loci', help='all loci regardless of signal', required=True)

args = argparser.parse_args()

cOA_locus_list = args.loci.split(",")

results = {}

length_list = [150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 300, 400, 500, 600, 700]

print("Length cutoff is ", str(args.length_cutoff))

# Read in the unknown proteins file
unknown_proteins = pd.read_csv(args.unknown_proteins, sep='\t', header = 0)
#create new column length_AA that contains the length of the protein sequence
unknown_proteins['length_AA'] = unknown_proteins['sequence'].apply(lambda x: len(x))
small_proteins = unknown_proteins[unknown_proteins['length_AA'] <= int(args.length_cutoff)]

length_dist = {}

#for every length in length_list, extract the small proteins with length less than the length into a dictionary
for length in length_list:
    length_dist[str(length)] = unknown_proteins[unknown_proteins['length_AA'] < length]

print("Length distribution of small proteins: ")
for length, proteins in length_dist.items():
    print(length + ": " + str(len(proteins)))

print(small_proteins.head())
print("Number of small proteins in signal loci under preset cutoff: " + str(len(small_proteins)))

for loci in cOA_locus_list:
    #extract coa name from the string "loci". For example the string is base_path + "/91_coa_RN_explorer/ca4_loci.txt", and the coa is then ca4
    signal = str(loci.split("/")[-1].split("_")[0])
    signal = signal.strip()
    # Read in the locus file as list
    with open(loci, 'r') as f:
        locus_list = f.read().splitlines() #e.g. cA3 associated loci
    #remove whitespace from the signal

    print("Finding small unknown proteins in ", signal + " loci") #"loci" is for example the list of cA3 loci

    print("Number of loci: ", len(locus_list))

    print("First three loci: ", locus_list[:3])
    print("Last three loci: ", locus_list[-3:])

    # Filter the small proteins by current cOA associated loci
    small_proteins_cOA = small_proteins[small_proteins['locus_id'].isin(locus_list)]

    results[signal] = small_proteins_cOA

# Write each signal result to a separate tsv file
for signal, result in results.items():
    out_file = os.path.join(args.output_folder, f'{signal}_small_unk_proteins.tsv')
    result.to_csv(out_file, sep='\t', index=False)
    #also create .faa files for the proteins by extracting the sequences from the table and using protein name as header
    out_faa = os.path.join(args.output_folder, f'{signal}_small_unk_proteins.faa')
    print("Output file: ", out_faa)
    with open(out_faa, 'w') as f:
        for index, row in result.iterrows():
            f.write(f'>{row["id"]}\n{row["sequence"]}\n')


#write all unknown proteins under cutoff to a tsv and faa file
out_file = os.path.join(args.output_folder, 'all_small_unk_proteins.tsv')
small_proteins.to_csv(out_file, sep='\t', index=False)
out_faa = os.path.join(args.output_folder, 'all_small_unk_proteins.faa')
print("Output file: ", out_faa)
with open(out_faa, 'w') as f:
    for index, row in small_proteins.iterrows():
        f.write(f'>{row["id"]}\n{row["sequence"]}\n')