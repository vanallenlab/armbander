import os
from collections import defaultdict
import pandas as pd

START = 'START'
END = 'END'


def build_ref_genes_lookup():
    package_path = os.path.dirname(os.path.realpath(__file__))
    ref_genes_filepath = os.path.join(package_path, 'data/refGene_hg19.txt')
    ref_genes_df = pd.read_csv(ref_genes_filepath, sep='\t', header=0)

    ref_genes_lookup = defaultdict(list)

    i = 0

    segment_by_id = {}
    for index, row in ref_genes_df.iterrows():
        chrom = row['chrom']
        tx_start = int(row['txStart'])
        tx_end = int(row['txEnd'])

        symbol = row['name2']

        ref_genes_lookup[chrom].append({'id': index, 'position': tx_start, 'type': START, 'symbol': symbol})
        ref_genes_lookup[chrom].append({'id': index, 'position': tx_end, 'type': END, 'symbol': symbol})
        segment_by_id[index] = {'start': tx_start, 'end': tx_end, 'symbol': symbol}

        i += 1

    # Make sure all the gene segment starts and ends are sorted in ascending order for each chromosome
    for chromosome in ref_genes_lookup.keys():
        starts_and_ends = ref_genes_lookup[chromosome]
        sorted_starts_and_ends = sorted(starts_and_ends, key=lambda x: x['position'])
        ref_genes_lookup[chromosome] = sorted_starts_and_ends

    entries = []
    for chromosome, segments in ref_genes_lookup.items():
        current_segments = set()
        current_left_bookend = 0

        for segment in segments:
            segment_id = segment.get('id')
            segment_type = segment.get('type')
            position = segment.get('position')
            current_gene_symbols_in_this_section = \
                set([segment_by_id[segment_id].get('symbol') for segment_id in current_segments])
            entry = '{}\t{}\t{}\t{}'.format(chromosome, current_left_bookend, position,
                                            ','.join(list(current_gene_symbols_in_this_section)))

            if segment_type == START:
                if position > current_left_bookend:
                    current_left_bookend = position
                    entries.append(entry)
                else:
                    pass
                current_segments.add(segment_id)

            if segment_type == END:
                # Update our sliding window to be just where we left off
                if current_left_bookend < position:
                    current_left_bookend = position
                    entries.append(entry)

                current_segments.remove(segment_id)

    return entries


def main():
    entries = build_ref_genes_lookup()
    package_path = os.path.dirname(os.path.realpath(__file__))
    lookup_filepath = os.path.join(package_path, 'data/refGene_hg19_linear_lookup.txt')
    with open(lookup_filepath, 'w') as lookup_file:
        lookup_file.write('chrom\tstart\tend\tsection_genes\n')
        for entry in entries:
            lookup_file.write('{}\n'.format(entry))


if __name__ == '__main__':
    main()