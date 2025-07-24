'''

'''
import argparse
import sqlite3
from Bio import SeqIO

def process_proteins_in_batches(cursor, protein_ids, batch_size, is_ca3):
    count = 0
    cluster_counts = {}
    for i in range(0, len(protein_ids), batch_size):
        batch_protein_ids = protein_ids[i:i+batch_size]
        cursor.execute('''
            SELECT protein_id, cluster_id FROM proteins WHERE protein_id IN ({})
        '''.format(','.join(['?'] * len(batch_protein_ids))), batch_protein_ids)
        cluster_ids = cursor.fetchall()

        for _, cluster_id in cluster_ids:
            cluster_counts[cluster_id] = cluster_counts.get(cluster_id, 0) + 1
            count += 1

    if is_ca3:
        column_name = 'ca3_phage_count'
    else:
        column_name = 'non_ca3_phage_count'

    cursor.executemany('''
        UPDATE clusters SET {} = {} + ? WHERE cluster_id = ?
    '''.format(column_name, column_name), [(count, cluster_id) for cluster_id, count in cluster_counts.items()])

def main():
    parser = argparse.ArgumentParser(description='Perform enrichment analysis on CA3 and non-CA3 phages')
    parser.add_argument('--ca3_proteins')
    parser.add_argument('--non_ca3_proteins')
    parser.add_argument('--out', help='Output file for enrichment results')
    parser.add_argument('--sqlite', help="Path to the SQLite database (.db)")
    args = parser.parse_args()

    # Create a connection to the SQLite database
    conn = sqlite3.connect(args.sqlite)
    cursor = conn.cursor()

    # Fetch cluster IDs for CA3 proteins in batches
    ca3_protein_ids = [record.id for record in SeqIO.parse(args.ca3_proteins, "fasta")]
    process_proteins_in_batches(cursor, ca3_protein_ids, batch_size=500, is_ca3=True)

    # Fetch cluster IDs for non-CA3 proteins in batches
    non_ca3_protein_ids = [record.id for record in SeqIO.parse(args.non_ca3_proteins, "fasta")]
    process_proteins_in_batches(cursor, non_ca3_protein_ids, batch_size=500, is_ca3=False)

    # Commit the changes and close the connection
    conn.commit()
    conn.close()

if __name__ == '__main__':
    main()
