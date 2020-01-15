
import argparse
import hashlib
import gzip
import sqlite3
from pathlib import Path
import xml.etree.ElementTree as et


parser = argparse.ArgumentParser(
    description='Filter Uniprot\'s UniRef100 XML files to bacterial subsequences and create initial ups db.'
)
parser.add_argument('--taxonomy', action='store', help='Path to NCBI taxonomy node.dmp file.')
parser.add_argument('--xml', action='store', help='Path to UniRef xml file.')
parser.add_argument('--db', action='store', help='Path to Pickle db file.')
args = parser.parse_args()

taxonomy_path = Path(args.taxonomy).resolve()
xml_path = Path(args.xml).resolve()
db_path = Path(args.db).resolve()


def is_taxon_child(child, LCA, taxonomy):
    parent = taxonomy.get(child, None)
    while(parent is not None and parent != '1'):
        if(parent == LCA):
            return True
        else:
            parent = taxonomy.get(parent, None)
    return False


taxonomy = {}
with taxonomy_path.open() as fh:
    for line in fh:
        cols = line.split('\t|\t', maxsplit=2)
        taxonomy[cols[0]] = cols[1]


upss = {}  # unique protein sequences
ns = 'http://uniprot.org/uniref'
with sqlite3.connect(str(db_path), isolation_level='EXCLUSIVE') as conn:
    conn.execute('PRAGMA page_size = 4096;')
    conn.execute('PRAGMA cache_size = 100000;')
    conn.execute('PRAGMA locking_mode = EXCLUSIVE;')
    conn.execute("PRAGMA mmap_size = %i;" % 20 * 1024 * 1024 * 1024)
    conn.execute('PRAGMA synchronous = OFF;')
    conn.execute('PRAGMA journal_mode = OFF')
    conn.execute('PRAGMA threads = 2;')
    conn.execute('DROP TABLE IF EXISTS ups;')
    conn.execute('''CREATE TABLE ups (
        hash TEXT PRIMARY KEY,
        length INTEGER NOT NULL,
        uniref90_id TEXT,
        uniref100_id TEXT NOT NULL,
        uniparc_id TEXT,
        ncbi_nrp_id TEXT,
        uniprotkb_acc TEXT,
        gene TEXT,
        product TEXT
        ) WITHOUT ROWID;''')
    conn.commit()

    with gzip.open(str(xml_path), mode="rt") as fh:
        ups_entries = []
        i = 0
        for event, elem in et.iterparse(fh):
            # print(elem)
            if(elem.tag == "{%s}entry" % ns):
                # print(elem)
                if('Fragment' in elem.find("./{%s}name" % ns).text):  # skip protein fragments
                    continue
                tax_id = elem.find("./{%s}property[@type='common taxon ID']" % ns).attrib['value']
                # print("tax-id=%s" % tax_id)
                if is_taxon_child(tax_id, '2', taxonomy):
                    # print("tax-id=%s" % tax_id)
                    uniref100_id = elem.attrib['id']
                    seq_representative = elem.find("./{%s}representativeMember/{%s}sequence" % (ns, ns))
                    seq = seq_representative.text + '*'
                    # hash = hashlib.blake2b(seq.encode()).hexdigest()
                    hash = hashlib.md5(seq.encode()).hexdigest()
                    if(hash in upss):
                        raise Exception('duplicated hashes! hash=%s, seq-a-UniRef100_id: %s, seq-b-UniRef100_id: %s' % (hash, upss[hash]), uniref100_id)
                    else:
                        rep_member = elem.find("./{%s}representativeMember/{%s}dbReference" % (ns, ns))
                        ups = (
                            hash,  # aa_hash
                            int(seq_representative.attrib['length']),  # length
                            rep_member.find("./{%s}property[@type='UniRef90 ID']" % ns).attrib['value'][9:],  # uniref90
                            uniref100_id[10:],  # uniref100
                            rep_member.find("./{%s}property[@type='UniParc ID']" % ns).attrib['value'][3:],  # uniparc
                            rep_member.find("./{%s}property[@type='UniProtKB accession']" % ns).attrib['value'],  # uniprot_acc
                        )
                        upss[hash] = uniref100_id
                        ups_entries.append(ups)
                        i += 1

                        if((i % 1000000) == 0):
                            conn.executemany(
                                "INSERT INTO ups (hash, length, uniref90_id, uniref100_id, uniparc_id, uniprotkb_acc) VALUES (?,?,?,?,?,?)",
                                ups_entries
                            )
                            conn.commit()
                            print("\t... %i" % i)
                            ups_entries = []

        conn.executemany(
            "INSERT INTO ups (hash, length, uniref90_id, uniref100_id, uniparc_id, uniprotkb_acc) VALUES (?,?,?,?,?,?)",
            ups_entries
        )
        conn.commit()
        ups_entries = []
