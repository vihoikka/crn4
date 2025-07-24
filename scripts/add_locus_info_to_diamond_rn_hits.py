'''
Part of the ring nuclease snakemake pipeline.
Finds genomic coordinates of proteins that are hits in the diamond blast using the gff file.
Uses 3rd party program gffutils.

Example run from the snakemake rule:
python scripts/add_locus_info_to_diamond_rn_hits.py --diamond {output.info} --gff {input.gff}
Diamond is both input and output in this script (we overwrite the input file).
'''

import gffutils
import pandas as pd
import argparse
import os

def main():
    argparser = argparse.ArgumentParser()
    argparser.add_argument('--diamond', help='Input diamond blast results')
    argparser.add_argument('--gff', help='Input gff file')
    argparser.add_argument('--sample', help='Sample wildcard')
    argparser.add_argument('--out', help='Sample wildcard')
    args = argparser.parse_args()

    diamond = args.diamond
    gff = args.gff
    sample = args.sample
    out = args.out

    print("Starting add_locus_info_to_diamond_rn_hits.py")

    #read in diamond blast results if file has contents
    if os.stat(diamond).st_size != 0:
        diamond_df = pd.read_csv(diamond, sep='\t', header = None)
        diamond_df.columns = ['diamond_query', 'diamond_target', 'pident', 'diamond_length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'diamond_bitscore']
    else:
        #create empty dataframe if diamond file is empty
        diamond_df = pd.DataFrame(columns=['diamond_query', 'diamond_target', 'pident', 'diamond_length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'diamond_bitscore'])

    #add columns genomic_start and genomic_end
    diamond_df['genomic_start'] = ''
    diamond_df['genomic_end'] = ''
    diamond_df['diamond_contig'] = ''
    diamond_df['ncbi_annotation'] = ''

    #add host column that has value of sample for any row that exists
    diamond_df['host'] = sample
    
    print("Dataframe constructed")
    print("Saving only best hit for each protein")

    #save only the best hit for each protein
    diamond_df = diamond_df.sort_values(by=['diamond_query', 'diamond_bitscore'], ascending=[True, False])
    diamond_df = diamond_df.drop_duplicates(subset='diamond_query', keep='first')

    print("Best hits saved")

    #for each protein in the diamond blast results, find the genomic coordinates from the gff file.
    #First load the gff file using gffutils

    #create a gff database
    print("Creating gff database")
    db = gffutils.create_db(gff, dbfn=':memory:', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)


    #iterate over the diamond blast results
    print("Iterating over diamond blast results")
    for index, row in diamond_df.iterrows():
        #print("Processing protein", index, "of", len(diamond_df))
        #print("Name: ", row['diamond_query'])
        #get the query protein ID
        protein_id = row['diamond_query']
        #get the coordinates of the protein from the gff file
        for feature in db.features_of_type('CDS', order_by='start'):
            if protein_id in feature.id:
                #print("Found matching protein in gff file: ", protein_id)
                #add the protein coordinates directly to the diamond_df
                diamond_df.at[index, 'genomic_start'] = feature.start
                diamond_df.at[index, 'genomic_end'] = feature.end
                #find the contig from the first column of the gff file
                diamond_df.at[index, 'diamond_contig'] = feature.seqid
                #find the product annotation from the gff file
                diamond_df.at[index, 'ncbi_annotation'] = feature.attributes['product'][0]
                break

    #write the updated diamond_df to the input diamond file
    diamond_df.to_csv(out, sep='\t', index=False)

    print("Finished add_locus_info_to_diamond_rn_hits.py for sample ", sample)

if __name__ == '__main__':
    main()