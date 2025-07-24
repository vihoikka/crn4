import argparse
import pandas as pd
import subprocess
import sys
import os
from Bio import SeqIO, SearchIO
from hmmer_to_pandas import parse_hmmer_domtblout

argparser = argparse.ArgumentParser()
#args for input and output
argparser.add_argument('--input', type=str, metavar='input', required=True,
                    help='Multifasta of an effector')
argparser.add_argument('--output_basepath', type=str, required=True, metavar='output_basepath')
argparser.add_argument('--effector', type=str, required=True, metavar='effector name')
argparser.add_argument('--pfam_path', type=str, required=True, metavar='pfam hmm library')

args = argparser.parse_args()

input_effector_fasta = args.input
output_basepath = args.output_basepath
effector = args.effector
pfam_path = args.pfam_path

print("Starting effector hmmer")
print("Path to pfam db is " + str(pfam_path))

E_value_pfam = str(1e-2) #this is for the full sequence, not domain-specific

#read in the fasta file
effector_sequences = SeqIO.parse(open(input_effector_fasta),'fasta')

effector_df = pd.DataFrame(columns=["ID",
                                    "length",
                                    ])

print("Loaded effectors in to df")
print(effector_df)

#each row should be a domain that is found in an effector

#run hmmscan against pfam database
hmmout = output_basepath + "/" + effector + ".hmm"
print("Starting hmmscan...")
hmmscan_commands = ['hmmscan', '--domtblout', hmmout, '--cpu', '5', '-E', E_value_pfam, str(pfam_path), input_effector_fasta]
print(hmmscan_commands)

try:
    result = subprocess.run(hmmscan_commands, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    result.check_returncode()  # if command had non-zero exit status, this will raise CalledProcessError.
except subprocess.CalledProcessError as e:
    print(f"hmmscan command failed with exit status {e.returncode}, stderr output:\n{e.stderr.decode('utf-8')}")
else:
    print(f"hmmscan command output:\n{result.stdout.decode('utf-8')}")
