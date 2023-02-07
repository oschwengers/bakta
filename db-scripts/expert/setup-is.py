import argparse
import logging
import re

from pathlib import Path

from Bio import SeqIO


IS_RANK = 90
IS_MIN_IDENTITY = 95
IS_MIN_QUERY_COV = 90
IS_MIN_MODEL_COV = 90


parser = argparse.ArgumentParser(
    description='Import IS into aa seq expert system.'
)
parser.add_argument('--expert-sequence', action='store', dest='expert_sequences', required=True, help='Path to Bakta expert sequence file.')
parser.add_argument('--proteins', action='store', dest='proteins', required=True, help='Path to input proteins file.')
args = parser.parse_args()


expert_sequences_path = Path(args.expert_sequences).resolve()
proteins_path = Path(args.proteins).resolve()


logging.basicConfig(
    filename='bakta.db.log',
    filemode='a',
    format='%(name)s - EXPERT-SEQ - %(levelname)s - %(message)s',
    level=logging.DEBUG
)
log = logging.getLogger('IS')


ises = {}
with proteins_path.open() as fh_in:
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
print(f'\tstored IS transposases ids: {len(ises)}')


aa_seqs = 0
with expert_sequences_path.open('w') as fh_out:
    for id, ise in ises.items():
        gene = 'tnp'
        product = f"{ise['family']} family {ise['name']} transposase"
        dbxrefs = [f"IS:{ise['name']}"]
        orfs = ise['orfs']
        if(len(orfs) == 3 and orfs[0][1][:5] == orfs[2][1][:5]):  # check if ORF 3 starts with ORF 1 -> fusion ORF -> annotate ORF1/2 as transposase orfA/orfB
            fh_out.write(f">{ise['name']}-ORFA IS~~~{IS_RANK}~~~{IS_MIN_IDENTITY}~~~{IS_MIN_QUERY_COV}~~~{IS_MIN_MODEL_COV}~~~{gene}~~~{product} ORF A~~~{','.join(dbxrefs)}\n{orfs[0][1]}\n")
            log.info(
                'write seq: id=%s, rank=%i, id=%f, q-cov=%f, s-cov=%f, gene=%s, product=%s, dbxrefs=%s',
                f"{ise['name']}-ORFA", IS_RANK, IS_MIN_IDENTITY, IS_MIN_QUERY_COV, IS_MIN_MODEL_COV, gene, f"{product} ORF A", ','.join(dbxrefs)
            )
            fh_out.write(f">{ise['name']}-ORFB IS~~~{IS_RANK}~~~{IS_MIN_IDENTITY}~~~{IS_MIN_QUERY_COV}~~~{IS_MIN_MODEL_COV}~~~{gene}~~~{product} ORF B~~~{','.join(dbxrefs)}\n{orfs[1][1]}\n")
            log.info(
                'write seq: id=%s, rank=%i, id=%f, q-cov=%f, s-cov=%f, gene=%s, product=%s, dbxrefs=%s',
                f"{ise['name']}-ORFB", IS_RANK, IS_MIN_IDENTITY, IS_MIN_QUERY_COV, IS_MIN_MODEL_COV, gene, f"{product} ORF B", ','.join(dbxrefs)
            )
            fh_out.write(f">{ise['name']}-ORF IS~~~{IS_RANK}~~~{IS_MIN_IDENTITY}~~~{IS_MIN_QUERY_COV}~~~{IS_MIN_MODEL_COV}~~~{gene}~~~{product}~~~{','.join(dbxrefs)}\n{orfs[2][1]}\n")
            log.info(
                'write seq: id=%s, rank=%i, id=%f, q-cov=%f, s-cov=%f, gene=%s, product=%s, dbxrefs=%s',
                f"{ise['name']}-ORF", IS_RANK, IS_MIN_IDENTITY, IS_MIN_QUERY_COV, IS_MIN_MODEL_COV, gene, product, ','.join(dbxrefs)
            )
            aa_seqs += 3
        else:
            for idx, transposase in enumerate(orfs):
                (orf, seq) = transposase
                fh_out.write(f">{ise['name']}-{idx+1} IS~~~{IS_RANK}~~~{IS_MIN_IDENTITY}~~~{IS_MIN_QUERY_COV}~~~{IS_MIN_MODEL_COV}~~~{gene}~~~{product}~~~{','.join(dbxrefs)}\n{seq}\n")
                log.info(
                    'write seq: id=%s, rank=%i, id=%f, q-cov=%f, s-cov=%f, gene=%s, product=%s, dbxrefs=%s',
                    f"{ise['name']}-{idx+1}", IS_RANK, IS_MIN_IDENTITY, IS_MIN_QUERY_COV, IS_MIN_MODEL_COV, gene, product, ','.join(dbxrefs)
                )
                aa_seqs += 1
print(f'\tstored IS sequences: {aa_seqs}')
log.debug('summary: IS sequences=%i', aa_seqs)


print("\nsuccessfully setup BlastRules expert system!")
