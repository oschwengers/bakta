import argparse
import logging
import re
from pathlib import Path

from Bio import SeqIO


FASTA_HEADER_PATTERN = re.compile(r"(VFG\d{6})(?:\(gb\|[\w\.]+?\))? \(([\w\/\.\-\\']{2,20})\) (.+)\[(.+)\(([A-Z]{2,3}\d{3,4})\)\] (?:\[.+\])")
VFDB_RANK = 75
VFDB_MIN_IDENTITY = 90
VFDB_MIN_QUERY_COV = 90
VFDB_MIN_MODEL_COV = 90

parser = argparse.ArgumentParser(
    description='Import VFDB into aa seq expert system.'
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
log = logging.getLogger('VFDB')


aa_seqs = 0
print('import NCBI BlastRule proteins...')
with proteins_path.open() as fh_in, expert_sequences_path.open('a') as fh_out:
    for record in SeqIO.parse(fh_in, 'fasta'):
        m = FASTA_HEADER_PATTERN.match(record.description)
        if(m is None):
            print(f"pattern is none: {record.description}")
            continue
        vfdb_id = m.group(1)
        gene = m.group(2)
        product = m.group(3)
        vfdb_category_id = m.group(5)
        seq = str(record.seq).upper()
        dbxrefs = [
            f"VFDB:{vfdb_id}",
            f"VFDB:{vfdb_category_id}"
        ]
        if(vfdb_id is not None and product is not None):
            fh_out.write(f">{vfdb_id} VFDB~~~{VFDB_RANK}~~~{VFDB_MIN_IDENTITY}~~~{VFDB_MIN_QUERY_COV}~~~{VFDB_MIN_MODEL_COV}~~~{gene}~~~{product}~~~{','.join(dbxrefs)}\n")
            fh_out.write(f"{seq}\n")
            log.info(
                'write seq: id=%s, rank=%i, id=%f, q-cov=%f, s-cov=%f, gene=%s, product=%s, dbxrefs=%s',
                vfdb_id, VFDB_RANK, VFDB_MIN_IDENTITY, VFDB_MIN_QUERY_COV, VFDB_MIN_MODEL_COV, gene, product, ','.join(dbxrefs)
            )
            aa_seqs += 1
print(f'\tstored VFDB sequences: {aa_seqs}')
log.debug('summary: VFDB sequences=%i', aa_seqs)


print("\nsuccessfully setup BlastRules expert system!")