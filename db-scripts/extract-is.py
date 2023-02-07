import argparse
import logging

from pathlib import Path

from Bio import SeqIO

parser = argparse.ArgumentParser(
    description='Extract IS transposases from ISFinder.'
)
parser.add_argument('--input', action='store', help='Path to ISFinder file.')
parser.add_argument('--output', action='store', help='Path to transposase file.')
args = parser.parse_args()

input_path = Path(args.input).resolve()
output_path = Path(args.output).resolve()

ises = {}
transposases = 0
with input_path.open() as fh_in:
    for record in SeqIO.parse(fh_in, 'fasta'):
        id = record.id
        seq = str(record.seq).upper()
        descriptions = record.description.split(' ', 1)[1][3:-3].split('~~~')
        if(descriptions[1].lower() == 'transposase'  and  len(seq) > 0):
            (is_name, is_group, is_family, is_orf) = descriptions[0].split('_', 3)
            if(id in ises):
                ise = ises[id]
                orfs = ise['orfs']
                orfs.append((is_orf, seq))
            else:
                ises[id] = {
                    'name': is_name,
                    'group': is_group,
                    'family': is_family,
                    'orfs': [(is_orf, seq)]
                }

with output_path.open('w') as fh_out:
    for id, ise in ises.items():
        orfs = ise['orfs']
        if(len(orfs) == 3 and orfs[0][1][:5] == orfs[2][1][:5]):  # check if ORF 3 starts with ORF 1 -> fusion ORF -> annotate ORF1/2 as transposase orfA/orfB
            fh_out.write(f">{ise['name']}_{ise['group']}_{ise['family']}_ORFA\n{orfs[0][1]}\n")
            fh_out.write(f">{ise['name']}_{ise['group']}_{ise['family']}_ORFB\n{orfs[1][1]}\n")
            fh_out.write(f">{ise['name']}_{ise['group']}_{ise['family']}_ORF\n{orfs[2][1]}\n")
            transposases += 3
        else:
            for idx, transposase in enumerate(orfs):
                (orf, seq) = transposase
                fh_out.write(f">{ise['name']}_{ise['group']}_{ise['family']}_ORF{idx+1}\n{seq}\n")
                transposases += 1

print(f'\tstored transposase sequences: {transposases}')
