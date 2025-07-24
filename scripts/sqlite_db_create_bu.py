'''
Creates a SQLite database from a CD-HIT cluster file and a FASTA file containing cluster protein sequences.
First, the script creates a table to store the cluster information and a table to store the protein information.
'''

import argparse
import sqlite3
from tqdm import tqdm
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import concurrent.futures

def process_line(line, cursor, cluster_sequences):
        line = line.strip()
        if line.startswith('>Cluster'):
            cluster_id = int(line.split()[1])
            # Insert the cluster ID into the clusters table
            cursor.execute('''
                INSERT OR IGNORE INTO clusters (cluster_id, ca3_phage_count, non_ca3_phage_count, ref_sequence, ref_length, ref_description)
                VALUES (?, ?, ?, ?, ?, ?)
            ''', (cluster_id, 0, 0, "unknown", 0, "unknown"))
        else:
            protein_id = line.split('>')[1].split('...')[0]
            length = int(line.split(',')[0].split('\t')[1].split('aa')[0].strip())
            if "*" not in line:
                identity = line.split('at')[1].strip()
                # Load the description from the fasta file
                description = cluster_sequences.loc[cluster_sequences['protein_id'] == protein_id, 'description']
            else:
                identity = "Reference"
                # Set the cluster_id for this sequence in the cluster_sequences df
                cluster_sequences.loc[cluster_sequences['protein_id'] == protein_id, 'cluster_id'] = cluster_id
                # Get the sequence for this protein
                sequence = cluster_sequences.loc[cluster_sequences['protein_id'] == protein_id, 'sequence']
                # Insert the sequence into the clusters table
                cursor.execute('''
                    UPDATE clusters SET ref_sequence = ? WHERE cluster_id = ?
                ''', (str(sequence), cluster_id))
                # Same for ref_length and ref_description
                ref_length = len(sequence)
                ref_description = cluster_sequences.loc[cluster_sequences['protein_id'] == protein_id, 'description']
                cursor.execute('''
                    UPDATE clusters SET ref_length = ? WHERE cluster_id = ?
                ''', (ref_length, cluster_id))
                cursor.execute('''
                    UPDATE clusters SET ref_description = ? WHERE cluster_id = ?
                ''', (str(ref_description), cluster_id))

            # Insert the protein information into the proteins table
            cursor.execute('''
                INSERT OR IGNORE INTO proteins (protein_id, cluster_id, length, identity)
                VALUES (?, ?, ?, ?)
            ''', (protein_id, cluster_id, length, identity))

def main():
    parser = argparse.ArgumentParser(description='Perform enrichment analysis on CA3 and non-CA3 phages')
    parser.add_argument('--all_proteins_cluster', help='All phage proteins')
    parser.add_argument('--all_proteins_cluster_fasta', help='All phage proteins in fasta format')
    parser.add_argument('--out', help='Output file for sqlite db')
    args = parser.parse_args()

    # Create a connection to the SQLite database
    conn = sqlite3.connect(args.out)
    cursor = conn.cursor()

    # Create a table to store the cluster information
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS clusters (
            cluster_id INTEGER PRIMARY KEY,
            ca3_phage_count INTEGER,
            non_ca3_phage_count INTEGER,
            ref_sequence TEXT,
            ref_length INT,
            ref_description TEXT
        )
    ''')

    # Create a table to store the protein information
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS proteins (
            protein_id TEXT PRIMARY KEY,
            cluster_id INTEGER,
            length INTEGER,
            identity REAL,
            FOREIGN KEY (cluster_id) REFERENCES clusters (cluster_id)
        )
    ''')

    # Open the CD-HIT cluster file
    with open(args.all_proteins_cluster, 'r') as file:
        cluster_id = None
        total_lines = sum(1 for _ in file)

    # List to collect sequence data in dictionary format
    sequence_data = []

    # Populate the sequence list with the dictionaries created from the fasta file
    print("Generating protein sequences dataframe...")
    with open(args.all_proteins_cluster_fasta, 'r') as file:
        total_records = sum(1 for _ in SeqIO.parse(file, 'fasta'))
        file.seek(0)  # Reset file pointer to the beginning
        with tqdm(total=total_records, unit='record') as pbar:
            for record in SeqIO.parse(file, 'fasta'):
                protein_id = record.id
                sequence = str(record.seq)
                description = record.description
                sequence_data.append({'cluster_id': "unknown", 'protein_id': protein_id, 'sequence': sequence, 'description': description, 'length': len(sequence)})
                pbar.update(1)
    
    print("Protein sequences dataframe created.")

    # Create DataFrame from the list
    cluster_sequences = pd.DataFrame(sequence_data, columns=['cluster_id', 'protein_id', 'sequence', 'description', 'length'])

    #print a sample of cluster_sequences
    print(cluster_sequences.head())

    # save dataframe to file
    cluster_sequences.to_csv('cluster_sequences.csv', index=False)

    # next, we open the cd-hit cluster definition file and extract information which protein belongs to which cluster
    print("Going line by line through the cluster file, finding corresponding entry in the cluster_sequences dataframe and inserting into the database...")
    with open(args.all_proteins_cluster, 'r') as file:
        with tqdm(total=total_lines, unit='line') as pbar:
            for line in file:
                line = line.strip()
                if line.startswith('>Cluster'):
                    cluster_id = int(line.split()[1])
                    # Insert the cluster ID into the clusters table
                    cursor.execute('''
                        INSERT OR IGNORE INTO clusters (cluster_id, ca3_phage_count, non_ca3_phage_count, ref_sequence, ref_length, ref_description)
                        VALUES (?, ?, ?, ?, ?, ?)
                    ''', (cluster_id, 0, 0, "unknown", 0, "unknown"))
                else: #all other lines are proteins
                    protein_id = line.split('>')[1].split('...')[0]
                    #print(protein_id)
                    length = int(line.split(',')[0].split('\t')[1].split('aa')[0].strip())
                    if "*" not in line: #this symbol signifies a reference sequence. If lacking, this is a protein within the current cluster
                        identity = line.split('at')[1].strip()
                        #load the description from the fasta file
                        description = cluster_sequences.loc[cluster_sequences['protein_id'] == protein_id, 'description']
                    else: #this is a reference sequence for the current cluster
                        identity = "Reference"
                        
                        #set the cluster_id for this sequence in the cluster_sequences df
                        cluster_sequences.loc[cluster_sequences['protein_id'] == protein_id, 'cluster_id'] = cluster_id
                        
                        #get the sequence for this protein
                        sequence = cluster_sequences.loc[cluster_sequences['protein_id'] == protein_id, 'sequence']

                        #insert the sequence into the clusters table
                        cursor.execute('''
                            UPDATE clusters SET ref_sequence = ? WHERE cluster_id = ?
                        ''', (str(sequence), cluster_id))
                        #same for ref_length and ref_description
                        ref_length = len(sequence)
                        ref_description = cluster_sequences.loc[cluster_sequences['protein_id'] == protein_id, 'description']
                        cursor.execute('''
                            UPDATE clusters SET ref_length = ? WHERE cluster_id = ?
                        ''', (ref_length, cluster_id))
                        cursor.execute('''
                            UPDATE clusters SET ref_description = ? WHERE cluster_id = ?
                        ''', (str(ref_description), cluster_id))
                    
                    #print(length)
                    #print(identity)

                    # Insert the protein information into the proteins table
                    cursor.execute('''
                        INSERT OR IGNORE INTO proteins (protein_id, cluster_id, length, identity)
                        VALUES (?, ?, ?, ?)
                    ''', (protein_id, cluster_id, length, identity))


                pbar.update(1)

    # Commit the changes and close the connection
    conn.commit()
    conn.close()


if __name__ == '__main__':
    main()