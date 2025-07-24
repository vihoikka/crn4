'''
Cleans up the output of blastdbcmd by removing fasta headers that contain multiple > characters
python scripts/clean_up_blastdbcmd.py --in {output.all_proteins} --out {output.all_proteins}
'''

import argparse

arg_parser = argparse.ArgumentParser(description='Clean up blastdbcmd output')

arg_parser.add_argument('--in', dest='input', help='Input file')
arg_parser.add_argument('--out', dest='output', help='Output file')


args = arg_parser.parse_args()

with open(args.input, 'r') as in_file:
    with open(args.output, 'w') as out_file:
        sequence = ''
        for line in in_file:
            if line.startswith('>'):
                if sequence:
                    out_file.write(sequence + '\n')
                    sequence = ''
                if line.count('>') == 1:
                    out_file.write(line)
                else:
                    accession = line.split('>')[1].split()[0]
                    out_file.write('>' + accession + '\n')
            else:
                sequence += line.strip()
        if sequence:
            out_file.write(sequence + '\n')