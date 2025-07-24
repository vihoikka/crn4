'''
This script performs enrichment analysis on ca3 vs. non-ca3 phages.

First loads up sqlite database where each protein cluster has cluster_id, ca3_phage_count and non_ca3_phage_count.
These values have already been populated.
The script then creates a contingency table for each cluster and performs Fisher's exact test
to determine if there are any significant differences between the two groups.

'''

import argparse
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
import sqlite3
from PIL import Image
import matplotlib.pyplot as plt
from tqdm.auto import tqdm
import concurrent.futures
from statsmodels.stats.multitest import multipletests

def main():
    parser = argparse.ArgumentParser(description='Perform enrichment analysis on CA3 and non-CA3 phages')
    parser.add_argument('--sqlite', help='CA3 clusters')
    parser.add_argument('--phage_mastertable', help='Phage info table')
    parser.add_argument('--out')
    args = parser.parse_args()

    #load the phage_mastertable as pandas tsv
    phage_mastertable = pd.read_csv(args.phage_mastertable, sep='\t', header = 0)

    #calculate the number of phages with a number greater than 0 in column ca3_interactions
    total_ca3_phages = len(phage_mastertable[phage_mastertable['ca3_interactions'] > 0])
    total_nonca3_phages = len(phage_mastertable[phage_mastertable['ca3_interactions'] == 0])

    # Load the data from the SQLite database
    conn = sqlite3.connect(args.sqlite)
    cursor = conn.cursor()
    cursor.execute('''
        SELECT cluster_id, ca3_phage_count, non_ca3_phage_count, ref_sequence, ref_length, ref_description FROM clusters
    ''')
    data = cursor.fetchall()
    conn.close()

    # Create a pandas DataFrame from the data
    df = pd.DataFrame(data, columns=['cluster_id', 'ca3_phage_count', 'non_ca3_phage_count', 'ref_sequence', 'ref_length', 'ref_description'])
    df.set_index('cluster_id', inplace=True)

    #calculate the total number of phages in each group
    #total_ca3_phages = df['ca3_phage_count'].sum()
    #total_non_ca3_phages = df['non_ca3_phage_count'].sum()

    #add these values to df
    #df['ca3_phage_count_neg'] = total_ca3_phages - df['ca3_phage_count']
    #df['non_ca3_phage_count_neg'] = total_non_ca3_phages - df['non_ca3_phage_count']

    df['ca3_phage_count_neg'] = total_ca3_phages - df['ca3_phage_count']
    df['non_ca3_phage_count_neg'] = total_nonca3_phages - df['non_ca3_phage_count']

    #create a column that shows if the cluster is enriched for CA3 phages
    df['enriched'] = df['ca3_phage_count'] > df['non_ca3_phage_count']

    #create column with the difference
    df['difference'] = df['ca3_phage_count'] - df['non_ca3_phage_count']

    #only take enriched rows
    df = df[df['enriched']]

    contingency_tables = []
    for _, row in df.iterrows():
        table = np.array([
            [row['ca3_phage_count'], row['non_ca3_phage_count']],
            [row['ca3_phage_count_neg'], row['non_ca3_phage_count_neg']]
        ])
        contingency_tables.append(table)

    # Print a few contingency tables for verification
    print("Sample contingency tables:")
    for i in range(min(3, len(contingency_tables))):
        print(contingency_tables[i])

    #save the contingency table to a tsv file
    np.save('contingency_table.npy', contingency_tables)
    #then save one of the tables into a .tsv file
    np.savetxt('contingency_table.tsv', contingency_tables[0], delimiter='\t')

    print("Table complete. Performing Fisher's exact test")

    # Perform Fisher's exact test
    def calculate_p_value(table):
        _, p_value = fisher_exact(table, alternative='greater')
        return p_value

    # Calculate p-values without multithreading
    p_values = []
    for table in tqdm(contingency_tables, desc='Calculating p-values', unit='table'):
        p_value = calculate_p_value(table)
        p_values.append(p_value)

    
    # Add the p-values to the DataFrame
    df['p_value'] = p_values

    corrected_p_values = multipletests(p_values, method='fdr_bh')[1] # correcting using the Benjamini-Hochberg method
    df['p_value_corrected_BH'] = corrected_p_values

    #correct the p-values for multiple testing using Bonferroni correction
    df['p_value_value_corrected_Bonferroni'] = df['p_value'] * len(df) # this is the Bonferroni correction

    #sort the rows by ca3_phage_count from largest to smallest
    df = df.sort_values('difference', ascending=False)

    # Write the results to a CSV file
    df.to_csv(args.out, sep='\t')

    #create a version where non_ca3_phage_count is 0
    df_zero = df[df['non_ca3_phage_count'] == 0]

    # Write the results to a CSV file
    df_zero.to_csv('enriched_clusters_non_ca3_0.tsv', sep='\t')

    #Visualise the results using matplotlib

    # Create a histogram of the p-values
    plt.hist(p_values, bins=20, color='skyblue', edgecolor='black')
    plt.xlabel('P-value')
    plt.ylabel('Frequency')
    plt.title('Histogram of p-values')
    plt.savefig('p_values_histogram.png')

    # Create a scatter plot of the p-values
    plt.scatter(range(1, len(p_values) + 1), p_values, color='skyblue')
    plt.axhline(y=0.05, color='r', linestyle='--')
    plt.xlabel('Cluster ID')
    plt.ylabel('P-value')
    plt.title('Scatter plot of p-values')
    plt.savefig('p_values_scatter.png')

    # Create a volcano plot of the p-values
    plt.scatter(-np.log10(df['p_value']), df['ca3_phage_count'] - df['non_ca3_phage_count'], color='skyblue')
    plt.axhline(y=0, color='r', linestyle='--')
    plt.axvline(x=-np.log10(0.05), color='r', linestyle='--')
    plt.xlabel('-log10(P-value)')
    plt.ylabel('CA3 phage count - non-CA3 phage count')
    plt.title('Volcano plot of p-values')
    plt.savefig('volcano_plot.png')



if __name__ == '__main__':
    main()