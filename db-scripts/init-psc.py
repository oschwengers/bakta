
import argparse
import logging
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


logging.basicConfig(
    filename='bakta.db.log',
    filemode='a',
    format='%(name)s - UniRef90 - %(levelname)s - %(message)s',
    level=logging.DEBUG
)
log = logging.getLogger('PSC')


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
        i = 0
        for event, elem in et.iterparse(fh_xml):
            if(elem.tag == "{%s}entry" % ns):
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
                        uniref90_id = elem.attrib['id'][9:]  # remove 'UniRef90_' prefix
                        try:
                            product = rep_member.find("./{%s}property[@type='protein name']" % ns).attrib['value']
                            if product.lower() == 'hypothetical protein':
                                product = None
                            elif product.lower() == 'uncharacterized protein':
                                product = None
                        except:
                            product = None
                        try:
                            uniref50_id = rep_member.find("./{%s}property[@type='UniRef50 ID']" % ns).attrib['value']
                            uniref50_id = uniref50_id[9:]  # remove 'UniRef50_' prefix
                        except:
                            uniref50_id = None
                        seq = elem.find("./{%s}representativeMember/{%s}sequence" % (ns, ns)).text.upper()
                        psc = (
                            uniref90_id,
                            uniref50_id,  # UniRef50
                            product
                        )
                        log.debug(
                            'uniref90-id=%s, common-tax-id=%s, rep-tax-id=%s, uniref50-id=%s, product=%s', 
                            uniref90_id, common_tax_id, rep_member_tax_id, uniref50_id, product
                        )
                        fh_fasta.write(">%s\n%s\n" % (uniref90_id, seq))
                        conn.execute("INSERT INTO psc (uniref90_id, uniref50_id, product) VALUES (?,?,?)", psc)
                        log.info('INSERT INTO psc (uniref90_id, uniref50_id, product) VALUES (%s,%s,%s)', *psc)
                        i += 1
                        if((i % 1000000) == 0):
                            conn.commit()
                            print("\t... %i" % i)
                elem.clear()  # forstall out of memory errors
    conn.commit()
    print("\tparsed PSC: %d" % i)
    log.debug('summary: # PSC=%d', i)
print("\nsuccessfully initialized PSC table!")
