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
from concurrent.futures import ThreadPoolExecutor, as_completed
from queue import Queue

# Function to process a single line
def process_line(line, cluster_sequences):
    line = line.strip()
    results = []
    if line.startswith('>Cluster'):
        cluster_id = int(line.split()[1])
        print("Current cluster ID: ", str(cluster_id))
        results.append(('clusters', (cluster_id, 0, 0, "unknown", 0, "unknown")))
    else:  # all other lines are proteins
        protein_id = line.split('>')[1].split('...')[0]
        length = int(line.split(',')[0].split('\t')[1].split('aa')[0].strip())
        if "*" not in line:  # this symbol signifies a reference sequence
            identity = line.split('at')[1].strip()
            #description = cluster_sequences.loc[cluster_sequences['protein_id'] == protein_id, 'description'].values[0]
            results.append(('proteins', (protein_id, cluster_id, length, identity)))
        else:  # this is a reference sequence for the current cluster
            identity = "Reference"
            cluster_sequences.loc[cluster_sequences['protein_id'] == protein_id, 'cluster_id'] = cluster_id
            sequence = cluster_sequences.loc[cluster_sequences['protein_id'] == protein_id, 'sequence'].values[0]
            ref_length = len(sequence)
            ref_description = cluster_sequences.loc[cluster_sequences['protein_id'] == protein_id, 'description'].values[0]
            results.append(('update_clusters', (sequence, cluster_id, ref_length, ref_description)))
            results.append(('proteins', (protein_id, cluster_id, length, identity)))
    return results

# Function to handle database operations
def db_worker(queue, conn):
    cursor = conn.cursor()
    while True:
        item = queue.get()
        if item is None:
            break
        table, data = item
        if table == 'clusters':
            cursor.execute('''
                INSERT OR IGNORE INTO clusters (cluster_id, ca3_phage_count, non_ca3_phage_count, ref_sequence, ref_length, ref_description)
                VALUES (?, ?, ?, ?, ?, ?)
            ''', data)
        elif table == 'update_clusters':
            sequence, cluster_id, ref_length, ref_description = data
            cursor.execute('''
                UPDATE clusters SET ref_sequence = ? WHERE cluster_id = ?
            ''', (str(sequence), cluster_id))
            cursor.execute('''
                UPDATE clusters SET ref_length = ? WHERE cluster_id = ?
            ''', (ref_length, cluster_id))
            cursor.execute('''
                UPDATE clusters SET ref_description = ? WHERE cluster_id = ?
            ''', (str(ref_description), cluster_id))
        elif table == 'proteins':
            cursor.execute('''
                INSERT OR IGNORE INTO proteins (protein_id, cluster_id, length, identity)
                VALUES (?, ?, ?, ?)
            ''', data)
        queue.task_done()
    conn.commit()
    conn.close()

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
    print("Ready to start multithreaded database operations.")

    # Create DataFrame from the list
    cluster_sequences = pd.DataFrame(sequence_data, columns=['cluster_id', 'protein_id', 'sequence', 'description', 'length'])

    #print a sample of cluster_sequences
    print(cluster_sequences.head())

    # save dataframe to file
    cluster_sequences.to_csv('cluster_sequences.csv', index=False)


    # Count the total number of lines in the file for tqdm
    with open(args.all_proteins_cluster, 'r') as file:
        total_lines = sum(1 for line in file)

    # Thread-safe queue for database operations
    db_queue = Queue()
    print("Thread-safe queue created.")

    # Connect to SQLite database
    conn = sqlite3.connect(args.out)
    print("Connected to SQLite database.")

    # Start database worker thread
    db_thread = ThreadPoolExecutor(max_workers=1)
    db_thread.submit(db_worker, db_queue, conn)

    # Process lines with ThreadPoolExecutor
    print("Processing lines with ThreadPoolExecutor...")
    with ThreadPoolExecutor(max_workers=4) as executor:
        futures = []
        with open(args.all_proteins_cluster, 'r') as file:
            with tqdm(total=total_lines, unit='line') as pbar:
                for line in file:
                    futures.append(executor.submit(process_line, line, cluster_sequences))
                    pbar.update(1)

        for future in as_completed(futures):
            results = future.result()
            for result in results:
                db_queue.put(result)

    # Signal the database worker to exit
    db_queue.put(None)
    db_thread.shutdown(wait=True)


if __name__ == '__main__':
    main()