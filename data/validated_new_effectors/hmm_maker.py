import argparse
import subprocess
import glob
import os

print('This script takes in a folder of fasta files and outputs a folder of hmm files')

print("Reading user inputs...")
#user inputs (using args) for the folder that contains desired fastas, the folder to output the hmm files, the name of the final hmm fil
parser = argparse.ArgumentParser(description='This script takes in a folder of fasta files and outputs a folder of hmm files')
parser.add_argument('-i', '--input', help='Input folder of fasta files', required=True)
parser.add_argument('-n', '--name', help='Name of the final hmm file', required=True)
args = parser.parse_args()

#assign the user inputs to variables
print('Assigning user inputs to variables...')
input_folder = args.input
name = args.name

#create temp folder for concatenated fasta file
if not os.path.exists("intermediate_files"):
    os.makedirs("intermediate_files")

intermediate_path = os.path.join("intermediate_files",name)

if not os.path.exists(intermediate_path):
    os.makedirs(intermediate_path)

concatenated_fasta = os.path.join(intermediate_path,str(name + ".faa"))
concatenated_fasta_clustered = concatenated_fasta.replace(".faa", "_clustered.faa")
concatenated_fasta_clustered_aligned = concatenated_fasta_clustered.replace(".faa", "_aligned.afa")

#concatenate all .faa files in the input folder
print('Concatenating all .faa files in the input folder...')
filenames = glob.glob(os.path.join(input_folder, '*.faa'))
with open(concatenated_fasta, 'w') as outfile:
    print("Found " + str(len(filenames)) + " fasta files in the input folder")
    for fname in filenames:
        print('Writing ' + fname + ' to concatenated fasta file...')
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)

#cluster the concatenated fasta file using cd-hit with 95% identity
print('Clustering concatenated fasta file using cd-hit with 95% identity...')
subprocess.call(['cd-hit', '-i', concatenated_fasta, '-o', concatenated_fasta_clustered, '-c', '0.95', '-n', '5', '-d', '0', '-M', '0', '-T', '40'])

#align the concatenated fasta file using muscle
print('Aligning concatenated fasta file using muscle...')
subprocess.call(['muscle', '-align', concatenated_fasta_clustered, '-output', concatenated_fasta_clustered_aligned, '-threads', '40'])

#build the hmm file using hmmbuild
print('Building hmm file using hmmbuild...')

output_folder_hmm = os.path.join("profiles",name)

#create directory for output hmm files
if not os.path.exists(name):
    os.makedirs(output_folder_hmm)

#build the hmm profile
subprocess.call(['hmmbuild', os.path.join(output_folder_hmm, name + ".hmm"), concatenated_fasta_clustered_aligned])

#hmmpress the hmm file
print('hmmpressing the hmm file...')
subprocess.call(['hmmpress', os.path.join(output_folder_hmm, name + ".hmm")])

print('Done!')