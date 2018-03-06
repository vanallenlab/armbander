import os
import pandas as pd
from armband import get_genes_in_stretch

chr_index = 0
start_index = 1
end_index = 2
band_index = 3


def main():
    package_path = os.path.dirname(os.path.realpath(__file__))

    # read in linear gene lookup table
    gene_lookup_filepath = os.path.join(package_path, 'data/refGene_hg19_linear_lookup.txt')
    gene_lookup_df = pd.read_csv(gene_lookup_filepath, sep='\t', header=0)

    # read in cytoband data
    cytoband_filepath = os.path.join(package_path, 'data/cytoBand_hg19.txt')
    cytoband_df = pd.read_csv(cytoband_filepath, sep='\t', header=None)
    cytoband_df['genes'] = pd.Series([','.join(get_genes_in_stretch(row[chr_index], row[start_index], row[end_index], gene_lookup_df))
                                      for index, row in cytoband_df.iterrows()])

    # Write output
    cytoband_genes_filepath = os.path.join(package_path, 'data/cytoBand_with_genes_hg19.txt')
    cytoband_df.to_csv(cytoband_genes_filepath, "\t", header=False, index=False)


if __name__ == '__main__':
    main()