import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Rename headers in a fasta file')

parser.add_argument('--input', type=str, help='Input fasta file')
parser.add_argument('--output', type=str, help='Output fasta file')
parser.add_argument('--prefix', type=str, help='Prefix')

args = parser.parse_args()

#using biopython, open fasta and rename all using the prefix and running number

running_number = 0

with open(args.output, 'w') as out:
    for record in SeqIO.parse(args.input, 'fasta'):
        running_number += 1
        record.id = f'{args.prefix}_{running_number}'
        record.description = ''
        SeqIO.write(record, out, 'fasta')
