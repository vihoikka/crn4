'''
This script takes as input list of phages and outputs their Millard proteomes
Run command:
python scripts/separate_coa_proteomes_phages.py --ca3_phages {input.ca3_phages} --non_ca3_phages {input.non-ca3_phages} --all_phage_proteins {input.all_phage_proteins} --out_ca3_proteins {output.ca3_proteins} --out_non_ca3_proteins {output.non_ca3_proteins}
'''

import argparse
from Bio import SeqIO
from tqdm import tqdm
import pandas as pd

def main():
    print("Starting the cOA proteome separation script")
    parser = argparse.ArgumentParser(description='Separate phage proteomes into CA3 and non-CA3')
    parser.add_argument('--ca3_phages', help='List of CA3 phages')
    parser.add_argument('--non_ca3_phages', help='List of non-CA3 phages')
    parser.add_argument('--all_phage_proteins', help='All phage proteins')
    parser.add_argument('--out_ca3_proteins', help='Output file for CA3 phages')
    parser.add_argument('--out_non_ca3_proteins', help='Output file for non-CA3 phages')
    args = parser.parse_args()

    ca3_phages = pd.read_csv(args.ca3_phages, sep='\t', header=0)
    ca3_phages = set(ca3_phages["phage"])

    non_ca3_phages = pd.read_csv(args.non_ca3_phages, sep='\t', header=0)
    non_ca3_phages = set(non_ca3_phages["phage"])

    #print first 10 phages
    print("First 10 CA3 phages: {}".format(list(ca3_phages)[:10]))
    print("First 10 non-CA3 phages: {}".format(list(non_ca3_phages)[:10]))

    print("Found {} CA3 phages and {} non-CA3 phages".format(len(ca3_phages), len(non_ca3_phages)))
    print("Reading the phage proteins file")

    ca3_proteins = []
    non_ca3_proteins = []

    #read protein fasta file

    input("Press Enter to continue...")


    for record in tqdm(SeqIO.parse(args.all_phage_proteins, 'fasta'), desc='Processing proteins'):
        #print(record.id)
        phage_id = record.id.split(' ')[0].split('_')[0]
        #print("Processing phage: {}".format(phage_id))
        #print("Protein ID: {}".format(record.id))
        if phage_id in ca3_phages:
            ca3_proteins.append(record)
            #print("Added protein from CA3 phage")
        elif phage_id in non_ca3_phages:
            non_ca3_proteins.append(record)
            #print("Added protein from non-CA3 phage")
        #input("Press Enter to continue...")

    print("Finished separating the phage proteomes")

    print("Writing the CA3 phage proteins")
    with open(args.out_ca3_proteins, 'w') as f:
        SeqIO.write(ca3_proteins, f, 'fasta')

    print("Writing the non-CA3 phage proteins")
    with open(args.out_non_ca3_proteins, 'w') as f:
        SeqIO.write(non_ca3_proteins, f, 'fasta')

if __name__ == '__main__':
    main()