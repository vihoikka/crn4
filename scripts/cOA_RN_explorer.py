'''
This script is launched from the ring nuclease snakemake pipeline rule shown below:

rule cOA_RN_explorer:
    Looks at mastertable_2 output to find loci with effectors associated withn given cOA.
    Then examines small proteins in those loci that are potential cA3 ring nucleases.
    For now, we just focus on cA3.
    input:
        mastertable = rules.mastercombiner.output.final_info_table
    output:
        coa_explorer = base_path + "/91_coa_RN_explorer/coa_RN_explorer.tsv"
    threads: thread_small
    shell:
        python scripts/cOA_RN_explorer.py --mastertable {input.mastertable} --output {output.ca3_explorer} --coa ca3

'''

import argparse
import pandas as pd
import os

argparser = argparse.ArgumentParser()
argparser.add_argument('--mastertable', help='mastertable_v2', required=True)
argparser.add_argument('--output', help='output file', required=True)
argparser.add_argument('--coa', help='cOA to search for', required=True)
args = argparser.parse_args()

# Read in the mastertable
mastertable = pd.read_csv(args.mastertable, sep='\t')

#exclude rows where Cas10_length is less than 500
mastertable = mastertable[mastertable['Cas10_length'] > 500]

coa_list = args.coa.split(",")
#in the coa_list, replace _ with -
coa_list = [coa.replace("_", "-") for coa in coa_list]

for coa in coa_list:
    coa_loci = mastertable[mastertable[coa]]
    print(f'Found {len(coa_loci)} loci with effectors associated with {coa}')
    loci = coa_loci['locus'].tolist()
    print(loci)
    #write locus names to newline separated file
    output_file = args.output + "/" + coa + "_loci.txt"
    with open(output_file, 'w') as f:
        for locus in loci:
            f.write(locus + "\n")