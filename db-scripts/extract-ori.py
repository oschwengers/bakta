import argparse
import logging

from pathlib import Path


parser = argparse.ArgumentParser(
    description='Extract oriC sequences from DoriC.'
)
parser.add_argument('--doric', action='store', help='Path to oriC DoriC file.')
parser.add_argument('--fasta', action='store', help='Path to oriC fasta file.')
args = parser.parse_args()


doric_path = Path(args.doric).resolve()
fasta_path = Path(args.fasta)


logging.basicConfig(
    filename='bakta.db.log',
    filemode='a',
    format='%(name)s - EXPERT-SEQ - %(levelname)s - %(message)s',
    level=logging.DEBUG
)
log = logging.getLogger('ORI')


oris = 0
with doric_path.open() as fh_doric, fasta_path.open('w') as fh_fasta:
    for line in fh_doric:
        if(line.startswith('ORI') or line.startswith('pORI')):
            cols = line.strip().split(',')
            ori_id = cols[0]
            sequence = cols[-1].upper()
            fh_fasta.write(f'>{ori_id}\n{sequence}\n')
            oris += 1

print(f'\tstored ORI sequences: {oris}')
log.debug('summary: ORI sequences=%i', oris)
