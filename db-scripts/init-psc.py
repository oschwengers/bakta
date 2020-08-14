
import argparse
import gzip
import sqlite3
import xml.etree.ElementTree as et
from pathlib import Path


parser = argparse.ArgumentParser(
    description='Filter Uniprot\'s UniRef90 XML files to bacterial subsequences and init pc db.'
)
parser.add_argument('--taxonomy', action='store', help='Path to NCBI taxonomy node.dmp file.')
parser.add_argument('--xml', action='store', help='Path to UniRef xml file.')
parser.add_argument('--db', action='store', help='Path to Bakta sqlite3 db file.')
parser.add_argument('--fasta', action='store', help='Path to PSC fasta file.')
args = parser.parse_args()

taxonomy_path = Path(args.taxonomy).resolve()
xml_path = Path(args.xml).resolve()
db_path = Path(args.db)
fasta_path = Path(args.fasta)


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


print('parse & store PSC information...')
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

    with gzip.open(str(xml_path), mode='rt') as fh_xml, fasta_path.open(mode='wt') as fh_fasta:
        psc_entries = []
        i = 0
        for event, elem in et.iterparse(fh_xml):
            # print(elem)
            if(elem.tag == "{%s}entry" % ns):
                # print(elem)
                if('Fragment' not in elem.find("./{%s}name" % ns).text):  # skip protein fragments
                    tax_property = elem.find("./{%s}property[@type='common taxon ID']" % ns)
                    if(tax_property is not None):
                        tax_id = tax_property.attrib['value']
                        if is_taxon_child(tax_id, '2', taxonomy):
                            # print("tax-id=%s" % tax_id)
                            uniref90_id = elem.attrib['id'][9:]  # remove 'UniRef90_' prefix
                            rep_member = elem.find("./{%s}representativeMember/{%s}dbReference" % (ns, ns))
                            try:
                                prot_name = rep_member.find("./{%s}property[@type='protein name']" % ns).attrib['value']
                            except:
                                prot_name = ''
                            seq = elem.find("./{%s}representativeMember/{%s}sequence" % (ns, ns)).text.upper()
                            # cluster = Cluster(uniref90_id, len(seq), seq, prot_name, uniref100_id, uniparc_id)
                            cluster = (
                                uniref90_id,
                                '',  # gene
                                prot_name,  # product
                                '',  # EC
                                '',  # IS
                                '',  # COG id
                                '',  # COG category
                                ''  # GO ids
                            )
                            fh_fasta.write(">%s\n%s\n" % (uniref90_id, seq))
                            psc_entries.append(cluster)
                            i += 1

                            if((i % 1000000) == 0):
                                conn.executemany(
                                    "INSERT INTO psc (uniref90_id, gene, product, ec_id, is_id, cog_id, cog_category, go_ids) VALUES (?,?,?,?,?,?,?,?)",
                                    psc_entries
                                )
                                conn.commit()
                                print("\t... %i" % i)
                                psc_entries = []
                elem.clear()  # forstall out of memory errors
    conn.executemany(
        "INSERT INTO psc (uniref90_id, gene, product, ec_id, is_id, cog_id, cog_category, go_ids) VALUES (?,?,?,?,?,?,?,?)",
        psc_entries
    )
    conn.commit()
    psc_entries = []
    print("\tparsed PSC: %d" % i)
print("\nsuccessfully initialized PSC table!")
