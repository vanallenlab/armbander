from collections import defaultdict
import numpy as np
import pandas as pd
from argparse import ArgumentParser
import os


def get_band_info(cytoband_dict, chrom, start):
    """Retrieve the arm and band-level information given a chromosome and start position"""
    bands = cytoband_dict['chr{}'.format(chrom)]
    index = np.searchsorted([b['start'] for b in bands], start)
    band = bands[index - 1]
    return band


def add_band_and_arm_columns(input_df, chr_header, start_header, end_header):
    package_path = os.path.dirname(os.path.realpath(__file__))
    cytoband_filepath = os.path.join(package_path, 'data/cytoBand_hg19.txt')
    cytoband_df = pd.read_csv(cytoband_filepath, sep='\t', header=None)
    chr_index = 0
    start_index = 1
    band_index = 3
    cyto_band_dict = defaultdict(list)

    for index, row in cytoband_df.iterrows():
        chrom = row[chr_index]
        start = row[start_index]
        band = row[band_index]
        cyto_band_dict[chrom].append({'start': start, 'band': band, 'arm': band[0]})

    band_info = [get_band_info(cyto_band_dict, r[chr_header], r[start_header]) for index, r in
                 input_df.iterrows()]
    input_df['band'] = pd.Series([b['band'] for b in band_info], index=input_df.index)
    input_df['arm'] = pd.Series([b['arm'] for b in band_info], index=input_df.index)

    return input_df


def write_to_file(self, output_filename=None, output_dir=None):
    """Write the maf to a file"""
    if not output_dir:
        output_dir = os.path.dirname(self.maf_file_path)
    if not output_filename:
        output_filename = '{}.arm_and_band_info.{}'.format(os.path.basename(self.maf_file_path).split('.')[0],
                                                                        os.path.basename(self.maf_file_path).split('.')[1])
        self.maf.to_csv(os.path.join(output_dir, output_filename), "\t", header=True, index=False)


def main():
    parser = ArgumentParser(description='Add arm and band-level information columns')
    parser.add_argument('input_file', metavar='input_file', type=str)
    parser.add_argument('output_dir', metavar='output_dir', type=str)
    parser.add_argument('--chr_header', metavar='chr_header', type=str, default='Chromosome')
    parser.add_argument('--start_header', metavar='start_header', type=str, default='Start_Position(bp)')
    parser.add_argument('--end_header', metavar='end_header', type=str, default='End_Position(bp)')
    parser.add_argument('--skiprows', metavar='skiprows', type=int, default=0)


    args = parser.parse_args()

    input_file = args.input_file
    output_dir = args.output_dir
    chr_header= args.chr_header
    start_header = args.start_header
    end_header = args.end_header
    skiprows = args.skiprows

    filename = os.path.basename(input_file)

    filename_components = filename.split('.')
    filename_sans_extension = '.'.join(filename_components[0:len(filename_components)-1])
    extension = filename_components[-1]

    input_df = pd.read_csv(input_file, delimiter='\t', skiprows=skiprows, comment='#')
    updated_df = add_band_and_arm_columns(input_df, chr_header, start_header, end_header)

    output_file = os.path.join(output_dir, '{}.arm_band_annotated.{}'.format(filename_sans_extension, extension))
    updated_df.to_csv(output_file,
                      "\t", header=True, index=False)


if __name__ == '__main__':
    main()