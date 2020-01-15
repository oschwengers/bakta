
import argparse
import hashlib
import pickle
from Bio import SeqIO
from pathlib import Path


parser = argparse.ArgumentParser(
    description='Filter Uniprot\'s UniRef100 XML files to bacterial subsequences and create initial ups db.'
)
parser.add_argument('--taxonomy', action='store', help='Path to NCBI taxonomy node.dmp file.')
parser.add_argument('--db-ups', action='store', help='Path to UPS db file.')
parser.add_argument('--db-psc', action='store', help='Path to PSC db file.')
parser.add_argument('--nrp', action='store', help='Path to NCBI NRP fasta file.')
parser.add_argument('--pcla-proteins', action='store', help='Path to NCBI NRP PCLA proteins file.')
parser.add_argument('--pcla-clusters', action='store', help='Path to NCBI NRP PCLA clusters file.')
args = parser.parse_args()


db_path = Path(args.db).resolve()
with db_path.open('rb') as fh:
    upss = pickle.load(fh)

nrps = {}
ncbi_nrp_path = Path(args.nrp).resolve()
with ncbi_nrp_path.open() as fh:
    for record in SeqIO.parse(fh, 'fasta'):
        seq = str(record.seq)
        # hash = hashlib.blake2b(seq.encode()).hexdigest()
        hash = hashlib.md5(seq.encode()).hexdigest()
        ups = upss.get(hash, None)
        if(ups is not None):
            if(ups['length'] == len(seq)):
                ups['ncbi_nrp'] = record.id
                nrps[record.id] = ups
                print("%s -> %s (%d)" % (hash, record.id, len(seq)))
            else:
                print("wrong AA length! legnth=%d, hash=%s, rnp-id=%s" % (len(seq), hash, record.id))
        else:
            print("UPS not found! hash=%s, nrp-id=%s" % (hash, record.id))

ncbi_pcla_path = Path(args.pcla).resolve()
with ncbi_pcla_path.open() as fh:
    for line in fh:
        if line[0] == '#':  # skip first line
            continue
        pcla_id, nrp_id = line.strip().split('\t')[0:2]

with db_path.open('wb') as fh:
    pickle.dump(upss, fh, pickle.HIGHEST_PROTOCOL)
