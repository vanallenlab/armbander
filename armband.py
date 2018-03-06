from collections import defaultdict
import numpy as np
import pandas as pd
from argparse import ArgumentParser
import os
import sys


def get_band_info(cytoband_dict, chrom, position):
    """Retrieve the arm and band-level information given a chromosome and position"""
    bands = cytoband_dict['chr{}'.format(chrom)]
    index = np.searchsorted([b['start'] for b in bands], position)
    band = bands[index - 1]
    return band


def build_cytoband_dict():
    package_path = os.path.dirname(os.path.realpath(__file__))
    cytoband_filepath = os.path.join(package_path, 'data/cytoBand_hg19.txt')
    cytoband_df = pd.read_csv(cytoband_filepath, sep='\t', header=None)
    chr_index = 0
    start_index = 1
    end_index = 2
    band_index = 3
    cytoband_dict = defaultdict(list)

    for index, row in cytoband_df.iterrows():
        chrom = row[chr_index]
        start = row[start_index]
        end = row[end_index]
        band = row[band_index]
        cytoband_dict[chrom].append({'chrom': chrom, 'start': start, 'end': end, 'band': band, 'arm': band[0]})

    return cytoband_dict

def get_genes_in_band(chrom, band_start, band_end, ref_genes_lookup):


def add_band_and_arm_columns(input_df, chr_header, start_header, end_header, ref_genes_lookup=None):
    cytoband_dict = build_cytoband_dict()
    band_info = [get_band_info(cytoband_dict, r[chr_header], (int(r[end_header]) - int(r[start_header]))/2.0) for index, r in
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


def get_stats(arm_band_df, chr_header, label_column):
    """Generate summary of arm and band level events based on label columns given as parameters"""
    arm_level_summary_dict = defaultdict(int)
    band_level_summary_dict = defaultdict(int)
    label_counts = defaultdict(int)

    for index, row  in arm_band_df.iterrows():
        chrom = row[chr_header]
        label = row[label_column]
        arm_event = "arm_level\t{}{}\t{}".format(chrom, row['arm'], label)
        band_event = "band_level\t{}{}\t{}".format(chrom, row['band'].split('.')[0], label)
        arm_level_summary_dict[arm_event] += 1
        band_level_summary_dict[band_event] += 1
        label_counts["label_level\t\t{}".format(label)] += 1

    return arm_level_summary_dict, band_level_summary_dict, label_counts


def import_ref_genes_linear_lookup():
    package_path = os.path.dirname(os.path.realpath(__file__))
    ref_genes_lookup_filepath = os.path.join(package_path, 'data/refGene_hg19_linear_lookup.txt')
    ref_genes_lookup_df = pd.read_csv(ref_genes_lookup_filepath, sep='\t', header=0)
    return ref_genes_lookup_df


def main():
    parser = ArgumentParser(description='Add arm and band-level information columns')
    parser.add_argument('input_file', metavar='input_file', type=str)
    parser.add_argument('output_dir', metavar='output_dir', type=str)
    parser.add_argument('--chr_header', metavar='chr_header', type=str, default='Chromosome')
    parser.add_argument('--start_header', metavar='start_header', type=str, default='Start_Position(bp)')
    parser.add_argument('--end_header', metavar='end_header', type=str, default='End_Position(bp)')
    parser.add_argument('--skiprows', metavar='skiprows', type=int, default=0)
    parser.add_argument('--summary_label', metavar='summary_label', type=str, default=None)
    parser.add_argument('--include_genes', metavar='include_genes', type=bool, default=False)
    args = parser.parse_args()

    input_file = args.input_file
    output_dir = args.output_dir
    chr_header= args.chr_header
    start_header = args.start_header
    end_header = args.end_header
    skiprows = args.skiprows
    label = args.summary_label
    include_genes = args.include_genes

    filename = os.path.basename(input_file)

    filename_components = filename.split('.')
    filename_sans_extension = '.'.join(filename_components[0:len(filename_components)-1])
    extension = filename_components[-1]

    sys.stdout.write("Reading in input file {}...\n".format(filename))
    input_df = pd.read_csv(input_file, delimiter='\t', skiprows=skiprows, comment='#')

    ref_genes_lookup = None
    if include_genes:
        sys.stdout.write("Reading in ref genes linear lookup table...\n")
        ref_genes_lookup = import_ref_genes_linear_lookup()

    sys.stdout.write("Adding arm-level and band-level information...\n")
    updated_df = add_band_and_arm_columns(input_df, chr_header, start_header, end_header,
                                          ref_genes_lookup=ref_genes_lookup)

    output_file = os.path.join(output_dir, '{}.arm_band_annotated.{}'.format(filename_sans_extension, extension))

    if label:
        arm_level_stats, band_level_stats, label_counts = get_stats(updated_df, chr_header, label)
        stats_file_path = os.path.join(output_dir, '{}.arm_band_stats.{}'.format(filename_sans_extension, extension))
        stats_file = open(stats_file_path, 'w')

        stats_file.write('category\tlocation\tlabel\tcount\n')
        for stats_dict in [arm_level_stats, band_level_stats, label_counts]:
            for k, v in stats_dict.items():
                stats_file.write('{}\t{}\n'.format(k, v))

    updated_df.to_csv(output_file,
                      "\t", header=True, index=False)


if __name__ == '__main__':
    main()