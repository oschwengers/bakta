
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
parser.add_argument('--db', action='store', help='Path to Bakta sqlite3 db file.')
args = parser.parse_args()

taxonomy_path = Path(args.taxonomy).resolve()
xml_path = Path(args.xml).resolve()
db_path = Path(args.db)


def is_taxon_child(child, LCA, taxonomy):
    parent = taxonomy.get(child, None)
    while(parent is not None and parent != '1'):
        if(parent == LCA):
            return True
        else:
            parent = taxonomy.get(parent, None)
    return False


print('parse & store NCBI taxonomy...')
taxonomy = {}
with taxonomy_path.open() as fh:
    for line in fh:
        cols = line.split('\t|\t', maxsplit=2)
        taxonomy[cols[0]] = cols[1]
print("\tstored tax ids: %d" % len(taxonomy))


print('parse & store UPS information...')
upss = {}  # unique protein sequences
ns = 'http://uniprot.org/uniref'
with sqlite3.connect(str(db_path), isolation_level='EXCLUSIVE') as conn:
    conn.execute('PRAGMA page_size = 4096;')
    conn.execute('PRAGMA cache_size = 100000;')
    conn.execute('PRAGMA locking_mode = EXCLUSIVE;')
    conn.execute("PRAGMA mmap_size = %i;" % (20 * 1024 * 1024 * 1024))
    conn.execute('PRAGMA synchronous = OFF;')
    conn.execute('PRAGMA journal_mode = OFF')
    conn.execute('PRAGMA threads = 2;')
    conn.commit()

    with gzip.open(str(xml_path), mode="rt") as fh:
        i = 0
        for event, elem in et.iterparse(fh):
            # print(elem)
            if(elem.tag == "{%s}entry" % ns):
                # print(elem)
                if('Fragment' not in elem.find("./{%s}name" % ns).text):  # skip protein fragments
                    try:
                        common_tax_id = elem.find("./{%s}property[@type='common taxon ID']" % ns).attrib['value']
                    except:
                        common_tax_id = 1
                    rep_member = elem.find("./{%s}representativeMember/{%s}dbReference" % (ns, ns))
                    try:
                        rep_member_tax_id = rep_member.find("./{%s}property[@type='NCBI taxonomy']" % ns).attrib['value']
                    except:
                        rep_member_tax_id = 1
                    
                    if(is_taxon_child(common_tax_id, '2', taxonomy) or is_taxon_child(rep_member_tax_id, '2', taxonomy)):
                        # print("tax-id=%s" % tax_id)
                        uniref100_id = elem.attrib['id']
                        seq_representative = elem.find("./{%s}representativeMember/{%s}sequence" % (ns, ns))
                        seq = seq_representative.text.upper()
                        hash = hashlib.md5(seq.encode()).hexdigest()
                        if(hash in upss):
                            raise Exception('duplicated hashes! hash=%s, seq-a-UniRef100_id: %s, seq-b-UniRef100_id: %s' % (hash, upss[hash]), uniref100_id)
                        else:
                            uniref90_id = rep_member.find("./{%s}property[@type='UniRef90 ID']" % ns)
                            if(uniref90_id is not None):
                                uniref90_id = uniref90_id.attrib['value'][9:]
                            uniparc_id = rep_member.find("./{%s}property[@type='UniParc ID']" % ns)
                            if(uniparc_id is not None):
                                uniparc_id = uniparc_id.attrib['value'][3:]
                            uniprot_acc = rep_member.find("./{%s}property[@type='UniProtKB accession']" % ns)
                            if(uniprot_acc is not None):
                                uniprot_acc = uniprot_acc.attrib['value']
                            ups = (
                                hash,  # aa_hash
                                int(seq_representative.attrib['length']),  # length
                                uniref90_id,
                                uniref100_id[10:],
                                uniparc_id,
                                uniprot_acc
                            )
                            upss[hash] = uniref100_id
                            if(uniref90_id):  # UniRef90 cluster available -> skip gene and product
                                # conn.execute("INSERT INTO ups (hash, length, uniref90_id, uniref100_id, uniparc_id, uniprotkb_acc) VALUES (?,?,?,?,?,?)", ups)
                                print("INSERT INTO ups (hash, length, uniref90_id, uniref100_id, uniparc_id, uniprotkb_acc) VALUES (%s,%d,%s,%s,%s,%s)" % ups)
                            else:
                                try:
                                    prot_name = rep_member.find("./{%s}property[@type='protein name']" % ns).attrib['value']
                                    if prot_name.lower() == 'hypothetical protein':
                                        prot_name = None
                                    elif prot_name.lower() == 'Uncharacterized protein':
                                        prot_name = None
                                except:
                                    prot_name = None
                                ups += (prot_name, )
                                # conn.execute("INSERT INTO ups (hash, length, uniref90_id, uniref100_id, uniparc_id, uniprotkb_acc, product) VALUES (?,?,?,?,?,?,?)", ups)
                                print("INSERT INTO ups (hash, length, uniref90_id, uniref100_id, uniparc_id, uniprotkb_acc, product) VALUES (%s,%d,%s,%s,%s,%s,%s)" % ups)
                            i += 1
                            if((i % 1000000) == 0):
                                # conn.commit()
                                print("\t... %i" % i)
                elem.clear()  # forstall out of memory errors
        # conn.commit()
        print("\tparsed UPS: %d" % i)
print("\nsuccessfully initialized UPS table!")
