import os
from collections import defaultdict
import pandas as pd

chr_index = 0
start_index = 1
end_index = 2
band_index = 3


def get_genes_in_band(cytoband_row, gene_lookup_df):
    """Given a row from the cytoband dataframe, figure out which genes are included in the segment"""
    chrom = cytoband_row[chr_index]
    start = int(cytoband_row[start_index])
    end = int(cytoband_row[end_index])
    genes_in_section = gene_lookup_df[(gene_lookup_df['chrom'] == chrom) &
                                      (gene_lookup_df['start'] >= start) &
                                      (gene_lookup_df['end'] <= end) &
                                      (gene_lookup_df['section_genes'].notnull())]['section_genes']

    # Some rows have multiple genes, so this list is still not totally collapsed yet
    nested_gene_list = [glist.split(',') for glist in list(genes_in_section)]

    # We use a list instead of a set here so we can preserve the order of genes as they appear on the segment
    flat_gene_set = []
    for gene_sublist in nested_gene_list:
        for gene in gene_sublist:
            if gene not in flat_gene_set:
                flat_gene_set.append(gene)

    return ','.join(flat_gene_set)


def main():
    package_path = os.path.dirname(os.path.realpath(__file__))

    # read in linear gene lookup table
    gene_lookup_filepath = os.path.join(package_path, 'data/refGene_hg19_linear_lookup.txt')
    gene_lookup_df = pd.read_csv(gene_lookup_filepath, sep='\t', header=0)

    # read in cytoband data
    cytoband_filepath = os.path.join(package_path, 'data/cytoBand_hg19.txt')
    cytoband_df = pd.read_csv(cytoband_filepath, sep='\t', header=None)
    cytoband_df['genes'] = pd.Series([get_genes_in_band(row, gene_lookup_df) for index, row in cytoband_df.iterrows()])

    # Write output
    cytoband_genes_filepath = os.path.join(package_path, 'data/cytoBand_with_genes_hg19.txt')
    cytoband_df.to_csv(cytoband_genes_filepath, "\t", header=False, index=False)


if __name__ == '__main__':
    main()