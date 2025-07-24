'''
Outputs a table that maps the contigs from this genome to this genome.
Output table has columns host and contig. Host is always the input genome name.
Contig contains each contig ID from the fasta file.

Example run:
python scripts/contig_to_genome_mapping.py --genome {input.genome} --contigs {input.contigs} --out {output.mapping}
'''

import os
from Bio import SeqIO
import argparse
import pandas as pd

def main():
    argparser = argparse.ArgumentParser()
    argparser.add_argument('--genome', help='Input genome name')
    argparser.add_argument('--contigs', help='Input contigs fasta file')
    argparser.add_argument('--out', help='Output mapping file')

    args = argparser.parse_args()

    genome = args.genome
    contigs = args.contigs
    out = args.out

    #create pandas df for output file
    contig_map_df = pd.DataFrame(columns=['host', 'contig'])

    #read in contigs
    contig_ids = []
    for record in SeqIO.parse(contigs, 'fasta'):
        contig_ids.append(record.id)

    #add contig ids to df
    contig_map_df['contig'] = contig_ids

    #add host genome name to df
    contig_map_df['host'] = genome

    #write to output file
    contig_map_df.to_csv(out, sep='\t', index=False)

if __name__ == '__main__':
    main()