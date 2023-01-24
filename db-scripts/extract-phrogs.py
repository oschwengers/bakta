import argparse
import logging

from pathlib import Path

from Bio import SeqIO

parser = argparse.ArgumentParser(
    description='Extract annotated PHROG proteins from PHROG db.'
)
parser.add_argument('--annotation', action='store', help='Path to PHROG annotation file.')
parser.add_argument('--proteins', action='store', help='Path to PHROG fasta file.')
parser.add_argument('--filtered-proteins', action='store', dest='filtered_proteins', help='Path to filtered PHROG fasta file.')
args = parser.parse_args()


annotation_path = Path(args.annotation).resolve()
proteins_path = Path(args.proteins)
filtered_proteins_path = Path(args.filtered_proteins)


logging.basicConfig(
    filename='bakta.db.log',
    filemode='a',
    format='%(name)s - PHROG - %(levelname)s - %(message)s',
    level=logging.DEBUG
)
log = logging.getLogger('PSC')


phrogs = {}
with annotation_path.open() as fh_in:
    for line in fh_in:
        if not line.startswith('phrog'):
            (id, color, product, category) = line.strip().split('\t')
            if product != 'NA':
                phrogs[f'phrog_{id}'] = product
print(f'\tstored PHROGs: {len(phrogs)}')
log.debug('summary: PHROGs=%i', len(phrogs))

with proteins_path.open() as fh_in, filtered_proteins_path.open('w') as fh_out:
    for record in SeqIO.parse(fh_in, 'fasta'):
        id = record.id
        if id in phrogs:
            fh_out.write(f'>{id}\n')
            fh_out.write(str(record.seq).upper())
            fh_out.write('\n')
