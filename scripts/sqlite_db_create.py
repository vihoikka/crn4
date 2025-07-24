'''
Creates a SQLite database from a CD-HIT cluster file and a FASTA file containing cluster protein sequences
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

    # Populate the sequence list with the sequences from the fasta file
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

    # Create DataFrame from the list
    cluster_sequences = pd.DataFrame(sequence_data, columns=['cluster_id', 'protein_id', 'sequence', 'description', 'length'])

    #print a sample of cluster_sequences
    print(cluster_sequences.head())

    #save to tsv
    cluster_sequences.to_csv('cluster_sequences.tsv', sep='\t', index=False)





    input_file = args.all_proteins_cluster

    # Count the total number of lines in the file for tqdm
    with open(input_file, 'r') as file:
        total_lines = sum(1 for line in file)

    # Read the protein cluster info file into a DataFrame with a progress bar
    cluster_data = []
    with open(input_file, 'r') as file:
        cluster_id = None
        with tqdm(total=total_lines, unit='line') as pbar:
            for line in file:
                line = line.strip()
                if line.startswith('>Cluster'):
                    cluster_id = int(line.split()[1])
                    cluster_data.append((cluster_id, None, None, None, None, None))
                else:
                    protein_id = line.split('>')[1].split('...')[0]
                    length = int(line.split(',')[0].split('\t')[1].split('aa')[0].strip())
                    if "*" not in line:
                        identity = line.split('at')[1].strip()
                        cluster_data.append((cluster_id, protein_id, length, identity, None, None))
                    else:
                        identity = "Reference"
                        cluster_data.append((cluster_id, protein_id, length, identity, None, None))
                pbar.update(1)
    # Convert the list to a DataFrame
    cluster_df = pd.DataFrame(cluster_data, columns=['cluster_id', 'protein_id', 'length', 'identity', 'ref_sequence', 'ref_description'])

    #save as tsv
    cluster_df.to_csv('cluster_info.tsv', sep='\t', index=False)

    # Connect to SQLite database
    print("Creating DB")
    conn = sqlite3.connect(args.out)
    cursor = conn.cursor()

    # Insert clusters

    # Insert clusters in chunks
    clusters = cluster_df[['cluster_id']].drop_duplicates().assign(ca3_phage_count=0, non_ca3_phage_count=0, ref_sequence="unknown", ref_length=0, ref_description="unknown")
    #clusters.to_sql('clusters', conn, if_exists='append', index=False, method='multi')
    chunk_size = 1000  # Define a chunk size
    for start in tqdm(range(0, len(clusters), chunk_size), unit='chunk'):
        end = start + chunk_size
        clusters_chunk = clusters.iloc[start:end]
        clusters_chunk.to_sql('clusters', conn, if_exists='append', index=False)


    # Count the total number of lines in the file for tqdm
    with open(input_file, 'r') as file:
        total_lines = sum(1 for line in file)

    # Update reference sequences and descriptions
    ref_sequences = cluster_df[cluster_df['identity'] == 'Reference']
    with tqdm(total=ref_sequences.shape[0], unit='row') as pbar:
        for idx, row in ref_sequences.iterrows():
            protein_id = row['protein_id']
            cluster_id = row['cluster_id']
            try:
                sequence = cluster_sequences.loc[cluster_sequences['protein_id'] == protein_id, 'sequence'].values[0]
            except IndexError:
                print(f"Protein ID {protein_id} not found in cluster sequences")
                print("The entry in full is ", cluster_sequences.loc[cluster_sequences['protein_id'] == protein_id])
                exit()
                
            ref_length = len(sequence)
            ref_description = cluster_sequences.loc[cluster_sequences['protein_id'] == protein_id, 'description'].values[0]

            cursor.execute('''
                UPDATE clusters SET ref_sequence = ?, ref_length = ?, ref_description = ? WHERE cluster_id = ?
            ''', (sequence, ref_length, ref_description, cluster_id))
            pbar.update(1)

    # Insert proteins
    proteins = cluster_df[['protein_id', 'cluster_id', 'length', 'identity']].dropna(subset=['protein_id'])
    proteins.to_sql('proteins', conn, if_exists='append', index=False, method='multi')

    # Commit and close the connection
    conn.commit()
    conn.close()



if __name__ == '__main__':
    main()