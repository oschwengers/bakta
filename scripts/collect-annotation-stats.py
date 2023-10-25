#!/usr/bin/python3

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
                'Translation tale',
                '# Sequences',
                'Size',
                'GC',
                'N ratio',
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
    for genome in args.genomes:
        genome_path = Path(genome).resolve()
        try:
            with genome_path.open() as fh_in:
                genome = json.load(fh_in)
            stats = [
                genome_path.stem,
                f"{' '.join([t for t in [genome['genome'].get('genus', None), genome['genome'].get('species', None), genome['genome'].get('strain', None)] if t is not None])}",
                'y' if genome['genome']['complete'] else 'n',
                f"{genome['genome']['translation_table']}",
                f"{genome['stats']['no_sequences']}",
                f"{genome['stats']['size']}",
                f"{100 * genome['stats']['gc']:.1f}",
                f"{100 * genome['stats']['n_ratio']:.1f}",
                f"{genome['stats']['n50']}",
                f"{100 * genome['stats']['coding_ratio']:.1f}",
                f"{len([f for f in genome['features'] if f['type'] == bc.FEATURE_T_RNA])}",
                f"{len([f for f in genome['features'] if f['type'] == bc.FEATURE_TM_RNA])}",
                f"{len([f for f in genome['features'] if f['type'] == bc.FEATURE_R_RNA])}",
                f"{len([f for f in genome['features'] if f['type'] == bc.FEATURE_NC_RNA])}",
                f"{len([f for f in genome['features'] if f['type'] == bc.FEATURE_NC_RNA_REGION])}",
                f"{len([f for f in genome['features'] if f['type'] == bc.FEATURE_CRISPR])}",
                f"{len([f for f in genome['features'] if f['type'] == bc.FEATURE_CDS])}",
                f"{len([f for f in genome['features'] if f['type'] == bc.FEATURE_CDS and 'hypothetical' in f])}",
                f"{len([f for f in genome['features'] if f['type'] == bc.FEATURE_CDS and 'pseudogene' in f])}",
                f"{len([f for f in genome['features'] if f['type'] == bc.FEATURE_SORF])}",
                f"{len([f for f in genome['features'] if f['type'] == bc.FEATURE_GAP])}",
                f"{len([f for f in genome['features'] if f['type'] == bc.FEATURE_ORIC])}",
                f"{len([f for f in genome['features'] if f['type'] == bc.FEATURE_ORIV])}",
                f"{len([f for f in genome['features'] if f['type'] == bc.FEATURE_ORIT])}",
            ]
            output_line = '\t'.join(stats)
            print(output_line)
            fh_out.write(f'{output_line}\n')
        except:
            print(f"Error reading genome {genome_path.stem}")
