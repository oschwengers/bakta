import argparse
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

with xopen(str(pfam_path)) as fh_pfam, family_path.open('w') as fh_family, non_family_path.open('w') as fh_non_family:
    entries = fh_pfam.read().split('//')
    for entry_text in entries:
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
            if(cols[1] == 'AC'):
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
                print(f'{type}\t{acc}\t{desc}')
                fh_family.write(f'{acc}\t{desc}\n')
            else:
                print(f'{type}\t{acc}\t{desc}')
                fh_non_family.write(f'{acc}\t{desc}\n')
