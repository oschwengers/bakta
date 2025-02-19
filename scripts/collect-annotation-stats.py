#!/usr/bin/env python3

##################################################################################
#
# This script collects annotation information from multiple JSON files
# and exports key metrics as TSV
#
# Thanks to Ahmed M. A. Elsherbini (https://github.com/AhmedElsherbini), for contributing.
#
##################################################################################


import argparse
import json
import os

from pathlib import Path

import bakta
import bakta.constants as bc


parser = argparse.ArgumentParser(
    prog=f'collect-annotation-stats',
    description='Collect annotation statistics and export as TSV',
    epilog=f'Version: {bakta.__version__}\nDOI: {bc.BAKTA_DOI}\nURL: github.com/oschwengers/bakta\n\nCitation:\n{bc.BAKTA_CITATION}',
    formatter_class=argparse.RawDescriptionHelpFormatter,
    add_help=False
)
parser.add_argument('genomes', metavar='<genomes>', nargs='+', help='Bakta genome annotation files in JSON format')
parser.add_argument('--prefix', '-p', action='store', default='annotation-stats', help='Prefix for output file')
parser.add_argument('--output', '-o', action='store', default=os.getcwd(), help='Output directory (default = current working directory)')
args = parser.parse_args()


prefix = args.prefix
output_path = Path(args.output).resolve().joinpath(f'{prefix}.tsv')
with output_path.open('w') as fh_out:
    fh_out.write(
        '\t'.join(
            [
                'Genome',
                'Taxon',
                'Complete',
                'Translation table',
                '# Sequences',
                'Size',
                'GC',
                'N ratio',
                'N50',
                'Coding ratio',
                'tRNA',
                'tmRNA',
                'rRNA',
                'ncRNA',
                'ncRNA region',
                'CRISPR',
                'CDS',
                'CDS hypothetical',
                'CDS pseudogene',
                'sORF',
                'GAP',
                'oriC',
                'oriV',
                'oriT'
            ]
        )
    )
    fh_out.write('\n')
    for genome_file in args.genomes:
        genome_path = Path(genome_file).resolve()
        try:
            with genome_path.open() as fh_in:
                data = json.load(fh_in)
            stats = [
                genome_path.stem,
                f"{' '.join([t for t in [data['genome'].get('genus', None), data['genome'].get('species', None), data['genome'].get('strain', None)] if t is not None])}",
                'y' if data['genome']['complete'] else 'n',
                f"{data['genome']['translation_table']}",
                f"{len(data['sequences'])}",
                f"{data['stats']['size']}",
                f"{100 * data['stats']['gc']:.1f}",
                f"{100 * data['stats']['n_ratio']:.1f}",
                f"{data['stats']['n50']}",
                f"{100 * data['stats']['coding_ratio']:.1f}",
                f"{len([feat for feat in data['features'] if feat['type'] == bc.FEATURE_T_RNA])}",
                f"{len([feat for feat in data['features'] if feat['type'] == bc.FEATURE_TM_RNA])}",
                f"{len([feat for feat in data['features'] if feat['type'] == bc.FEATURE_R_RNA])}",
                f"{len([feat for feat in data['features'] if feat['type'] == bc.FEATURE_NC_RNA])}",
                f"{len([feat for feat in data['features'] if feat['type'] == bc.FEATURE_NC_RNA_REGION])}",
                f"{len([feat for feat in data['features'] if feat['type'] == bc.FEATURE_CRISPR])}",
                f"{len([feat for feat in data['features'] if feat['type'] == bc.FEATURE_CDS])}",
                f"{len([feat for feat in data['features'] if feat['type'] == bc.FEATURE_CDS and 'hypothetical' in feat])}",
                f"{len([feat for feat in data['features'] if feat['type'] == bc.FEATURE_CDS and 'pseudogene' in feat])}",
                f"{len([feat for feat in data['features'] if feat['type'] == bc.FEATURE_SORF])}",
                f"{len([feat for feat in data['features'] if feat['type'] == bc.FEATURE_GAP])}",
                f"{len([feat for feat in data['features'] if feat['type'] == bc.FEATURE_ORIC])}",
                f"{len([feat for feat in data['features'] if feat['type'] == bc.FEATURE_ORIV])}",
                f"{len([feat for feat in data['features'] if feat['type'] == bc.FEATURE_ORIT])}",
            ]
            output_line = '\t'.join(stats)
            print(output_line)
            fh_out.write(f'{output_line}\n')
        except:
            print(f"Error reading genome {genome_path.stem}")
