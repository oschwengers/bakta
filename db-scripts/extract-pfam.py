import argparse
import logging

from pathlib import Path

from xopen import xopen


parser = argparse.ArgumentParser(
    description='Extract Pfam family & domain accessions.'
)
parser.add_argument('--pfam', action='store', help='Path to PfamA dat file.')
parser.add_argument('--family', action='store', help='Path to family accession file.')
parser.add_argument('--non-family', action='store', dest='non_family', help='Path to non-family accession file.')
args = parser.parse_args()


pfam_path = Path(args.pfam).resolve()
family_path = Path(args.family)
non_family_path = Path(args.non_family)


logging.basicConfig(
    filename='bakta.db.log',
    filemode='a',
    format='%(name)s - PFAM - %(levelname)s - %(message)s',
    level=logging.DEBUG
)
log = logging.getLogger('PSC')


non_families = 0
with xopen(str(pfam_path)) as fh_pfam, family_path.open('w') as fh_family, non_family_path.open('w') as fh_non_family:
    entries = fh_pfam.read().split('//')
    for entry_text in entries:
        id = None
        acc = None
        desc = None
        type = None
        entry_text = entry_text.strip()
        if(entry_text == ''):
            continue
        for line in entry_text.splitlines():
            line = line.strip()
            if(line == ''):
                continue
            cols = line.split()
            if(cols[1] == 'ID'):
                id = cols[2]
            elif(cols[1] == 'AC'):
                acc = cols[2]
            elif(cols[1] == 'DE'):
                desc = ' '.join(cols[2:])
            if(cols[1] == 'TP'):
                type = cols[2].lower()
        skip = False
        for blacklist_term in ['viral', 'virus', 'vacuol']:
            if(blacklist_term in desc.lower()):
                skip = True
                break
        if(skip is False):
            if(type == 'family'):
                print(f'{type}\t{acc}\t{id}\t{desc}')
                fh_family.write(f'{acc}\t{id}\t{desc}\n')
            else:
                print(f'{type}\t{acc}\t{desc}')
                fh_non_family.write(f'{acc}\t{id}\t{desc}\n')
                non_families += 1

print(f'\tparsed non-family PFAM entries: {non_families}')
log.debug('summary: # PFAM non-familiy=%i', non_families)
