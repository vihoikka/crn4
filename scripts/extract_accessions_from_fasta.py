import argparse
from Bio import SeqIO

def extract_accessions(proteins_file, accessions_file, output_file):
    # Read accessions from the input file
    with open(accessions_file, "r") as f:
        accessions = set(line.strip() for line in f)
    
    #replace whitespaces with underscore
    accessions = {acc.replace(" ", "_") for acc in accessions}

    # Extract sequences with matching accessions from the input proteins file
    with open(proteins_file, "r") as input_f, open(output_file, "w") as output_f:
        for record in SeqIO.parse(input_f, "fasta"):
            if record.id in accessions:
                SeqIO.write(record, output_f, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract sequences with specific accessions from a FASTA file.")
    parser.add_argument("--proteins", required=True, help="Input FASTA file containing protein sequences.")
    parser.add_argument("--accessions", required=True, help="File containing the list of accessions to extract.")
    parser.add_argument("--output", required=True, help="Output FASTA file to store the extracted sequences.")
    args = parser.parse_args()

    extract_accessions(args.proteins, args.accessions, args.output)
