#!/usr/bin/python3

##################################################################################
#
# This script extract requested subregions of a genome
# and exports contained features as regular Bakta outputs.
#
##################################################################################


import argparse
import json
import os

from pathlib import Path

import bakta
import bakta.constants as bc
import bakta.io.gff as gff
import bakta.io.insdc as insdc
import bakta.io.fasta as fasta
import bakta.config as cfg


parser = argparse.ArgumentParser(
    prog=f'extract region',
    description='Extract genomic region with a given range and exports selected features as GFF3, FAA, FFN, EMBL and Genbank',
    epilog=f'Version: {bakta.__version__}\nDOI: {bc.BAKTA_DOI}\nURL: github.com/oschwengers/bakta\n\nCitation:\n{bc.BAKTA_CITATION}',
    formatter_class=argparse.RawDescriptionHelpFormatter,
    add_help=False
)
parser.add_argument('genome', metavar='<genome>', help='Bakta genome annotation in JSON format')
parser.add_argument('--prefix', '-p', action='store', default=None, help='Prefix for output files')
parser.add_argument('--output', '-o', action='store', default=os.getcwd(), help='Output directory (default = current working directory)')
parser.add_argument('--sequence', '-s', action='store', default=None, help='Sequence/Contig (default = first)')
parser.add_argument('--min', '-m', action='store', type=int, default=0, help='Left region border including (default = 0)')
parser.add_argument('--max', '-x', action='store', type=int, default=100_000_000, help='Right region border including (default = 100,000,000)')
args = parser.parse_args()


print('Load annotated genome...')
genome_path = Path(args.genome).resolve()
with genome_path.open() as fh:
    genome = json.load(fh)

sequence_id = args.sequence
if(sequence_id is None):  # take first sequence as default
    sequence_id = genome['sequences'][0]['id']

prefix = args.prefix
if(prefix is None):  # use input file prefix as default
    prefix = genome_path.stem


print('Extract features within selected region...')
features_selected = []
for feat in genome['features']:
    if(feat['sequence'] == sequence_id):
        if(feat['start'] >= args.min  and  feat['stop'] <= args.max):
            features_selected.append(feat)
features_by_sequence = {sequence_id: features_selected}  # needed for GFF3 export
print(f'\t...selected features: {len(features_selected)}')

genome['features'] = features_selected
genome['sequences'] = [sequence for sequence in genome['sequences'] if sequence['id'] == sequence_id]
genome['genus'] = genome['genome']['genus']
genome['species'] = genome['genome']['species']
genome['strain'] = genome['genome']['strain']
genome['taxon'] = f"{genome['genome']['genus']} {genome['genome']['species']} {genome['genome']['strain']}"
cfg.db_info = {
    'major': genome['version']['db']['version'].split('.')[0],
    'minor': genome['version']['db']['version'].split('.')[1],
    'type': genome['version']['db']['type']
}

print('Write selected features...')
output_path = Path(args.output).resolve()
gff3_path = output_path.joinpath(f'{prefix}.gff3')
gff.write_features(genome, features_by_sequence, gff3_path)
print('\t...INSDC GenBank & EMBL')
genbank_path = output_path.joinpath(f'{prefix}.gbff')
embl_path = output_path.joinpath(f'{prefix}.embl')
insdc.write_features(genome, features_selected, genbank_path, embl_path)
print('\t...feature nucleotide sequences')
ffn_path = output_path.joinpath(f'{prefix}.ffn')
fasta.write_ffn(features_selected, ffn_path)
print('\t...translated CDS sequences')
faa_path = output_path.joinpath(f'{prefix}.faa')
fasta.write_faa(features_selected, faa_path)