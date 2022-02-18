import json
import os

from pathlib import Path

from matplotlib.font_manager import json_load

import bakta.utils as bu
import bakta.io.gff as gff
import bakta.io.insdc as insdc
import bakta.io.fasta as fasta
import bakta.config as cfg

parser = bu.init_parser()
parser.add_argument('genome', metavar='<genome>', help='Bakta genome annotation in JSON format')
parser.add_argument('--prefix', '-p', action='store', default=None, help='Prefix for output files')
parser.add_argument('--output', '-o', action='store', default=os.getcwd(), help='Output directory (default = current working directory)')
parser.add_argument('--contig', '-c', action='store', default=None, help='Sequence/Contig (default = first)')
parser.add_argument('--min', '-m', action='store', type=int, default=0, help='Left region border including (default = 0)')
parser.add_argument('--max', '-x', action='store', type=int, default=0, help='Right region border including (default = genome size)')
args = parser.parse_args()


genome_path = Path(args.genome).resolve()
genome = json_load(genome_path)


contig = args.contig
if(contig is None):
    contig = genome['sequences'][0]['id']

prefix = args.prefix
if(prefix is None):
    prefix = genome_path.stem


features_by_contig = {k['id']: [] for k in genome['sequences']}
features_selected = []
for feat in genome['features']:
    contig_features = features_by_contig.get(feat['contig'])
    contig_features.append(feat)
    # if(feat['start'] >= args.min  and  feat['stop'] <= args.max):
    if(feat['contig'] == contig  and  feat['start'] >= args.min  and  feat['stop'] <= args.max):
        # print(f"select feature: {feat}")
        features_selected.append(feat)

genome['features'] = features_selected
genome['contigs'] = genome['sequences']

genome['genus'] = genome['genome']['genus']
genome['species'] = genome['genome']['species']
genome['strain'] = genome['genome']['strain']
genome['taxon'] = f"{genome['genome']['genus']} {genome['genome']['species']} {genome['genome']['strain']}"

cfg.db_info = {
    'major': genome['version']['db'].split('.')[0],
    'minor': genome['version']['db'].split('.')[1]
}

output_path = Path(args.output).resolve()
gff3_path = output_path.joinpath(f'{prefix}.gff3')
gff.write_gff3(genome, features_by_contig, gff3_path)
print('INSDC GenBank & EMBL...')
genbank_path = output_path.joinpath(f'{prefix}.gbff')
embl_path = output_path.joinpath(f'{prefix}.embl')
insdc.write_insdc(genome, features_selected, genbank_path, embl_path)
print('feature nucleotide sequences...')
ffn_path = output_path.joinpath(f'{prefix}.ffn')
fasta.write_ffn(features_selected, ffn_path)
print('translated CDS sequences...')
faa_path = output_path.joinpath(f'{prefix}.faa')
fasta.write_faa(features_selected, faa_path)