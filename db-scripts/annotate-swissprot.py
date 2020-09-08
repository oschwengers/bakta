import argparse
import logging
import hashlib
import gzip
import sqlite3
from pathlib import Path
import xml.etree.ElementTree as et


parser = argparse.ArgumentParser(
    description='Filter Uniprot\'s UniRef100 XML files to bacterial subsequences and create initial ups db.'
)
parser.add_argument('--taxonomy', action='store', help='Path to NCBI taxonomy node.dmp file.')
parser.add_argument('--xml', action='store', help='Path to SwissProt xml file.')
parser.add_argument('--db', action='store', help='Path to Bakta sqlite3 db file.')
args = parser.parse_args()


taxonomy_path = Path(args.taxonomy).resolve()
xml_path = Path(args.xml).resolve()
db_path = Path(args.db)


logging.basicConfig(
    filename='bakta.db.log',
    filemode='a',
    format='%(name)s - SwissProt - %(levelname)s - %(message)s',
    level=logging.INFO
)
logger_ups = logging.getLogger('UPS')
logger_psc = logging.getLogger('PSC')


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

sp_processed = 0
sp_not_found = 0
sp_lenth_mm = 0

ups_annotated = 0
psc_annotated = 0

print('parse & store SwissProt information...')
ns = 'http://uniprot.org/uniprot'
with sqlite3.connect(str(db_path), isolation_level='EXCLUSIVE') as conn:
    conn.execute('PRAGMA page_size = 4096;')
    conn.execute('PRAGMA cache_size = 100000;')
    conn.execute('PRAGMA locking_mode = EXCLUSIVE;')
    conn.execute("PRAGMA mmap_size = %i;" % (20 * 1024 * 1024 * 1024))
    conn.execute('PRAGMA synchronous = OFF;')
    conn.execute('PRAGMA journal_mode = OFF')
    conn.execute('PRAGMA threads = 2;')
    conn.commit()

    conn.row_factory = sqlite3.Row
    c = conn.cursor()
    with gzip.open(str(xml_path), mode="rt") as fh:
        ups_entries = []
        i = 0
        for event, elem in et.iterparse(fh):
            # print(elem)
            if(elem.tag == "{%s}entry" % ns):
                # print(elem)
                sp_processed += 1
                tax_property = elem.find("./{%s}organism/{%s}dbReference[@type='NCBI Taxonomy']" % (ns, ns))
                if(tax_property is not None):
                    tax_id = tax_property.get('id')
                    if is_taxon_child(tax_id, '2', taxonomy):
                        acc = elem.find("./{%s}accession" % ns).text
                        gene = elem.find("./{%s}gene/{%s}name[@type='primary']" % (ns, ns))
                        gene = '' if gene is None else gene.text
                        product = elem.find("./{%s}protein/{%s}recommendedName/{%s}fullName" % (ns, ns, ns))
                        product = '' if product is None else product.text
                        ec_id = elem.find("./{%s}protein/{%s}recommendedName/{%s}ecNumber" % (ns, ns, ns))
                        ec_id = '' if ec_id is None else ec_id.text
                        seq = elem.find("./{%s}sequence" % ns).text.upper()
                        hash = hashlib.md5(seq.encode()).hexdigest()
                        go_terms = []
                        for go_term in elem.findall("./{%s}dbReference[@type='GO']" % ns):
                            go_terms.append(go_term.get('id'))
                        try:
                            c.execute('SELECT * FROM ups WHERE hash=?', (hash,))
                            r = c.fetchone()
                            if(r is None):
                                print("UPS hash not found! swissprot-id=%s, hash=%s" % (acc, hash))
                                sp_not_found += 1
                                continue
                            elif(r['length'] != len(seq)):
                                print("UPS length mismatch! swissprot-id=%s, hash=%s, ups-length=%d, swissprot-length=%d" % (acc, hash, r['length'], len(seq)))
                                sp_lenth_mm += 1
                                continue
                            else:
                                logger_ups.debug(
                                    'swissprot-id=%s, uniref100-id=%s, hash=%s, uniref90-id=%s, gene=%s, product=%s, ec=%s', 
                                    acc, r['uniref100_id'], hash, r['uniref90_id'], gene, product, ec_id
                                )
                                c.execute('SELECT * FROM psc WHERE uniref90_id=?', (r['uniref90_id'],))
                                r = c.fetchone()
                                if(r is None):  # no PSC found -> annotate UPS
                                    update = False
                                    if(gene):
                                        c.execute('UPDATE ups SET gene=? WHERE hash=?', (gene, hash))
                                        logger_ups.info('UPDATE ups SET gene=%s WHERE hash=%s', gene, hash)
                                        update = True
                                    if(product):
                                        c.execute('UPDATE ups SET product=? WHERE hash=?', (product, hash))
                                        logger_ups.info('UPDATE ups SET product=%s WHERE hash=%s', product, hash)
                                        update = True
                                    if(ec_id):
                                        c.execute('UPDATE ups SET ec_id=? WHERE hash=?', (ec_id, hash))
                                        logger_ups.info('UPDATE ups SET ec_id=%s WHERE hash=%s', ec_id, hash)
                                        update = True
                                    if(update):
                                        ups_annotated += 1
                                else:  # PSC found -> annotate PSC
                                    update = False
                                    if(gene):
                                        c.execute('UPDATE psc SET gene=? WHERE uniref90_id=?', (gene, r['uniref90_id']))
                                        logger_psc.info('UPDATE psc SET gene=%s WHERE uniref90_id=%s', gene, r['uniref90_id'])
                                        update = True
                                    if(product):
                                        c.execute('UPDATE psc SET product=? WHERE uniref90_id=?', (product, r['uniref90_id']))
                                        logger_psc.info('UPDATE psc SET product=%s WHERE uniref90_id=%s', product, r['uniref90_id'])
                                        update = True
                                    if(ec_id):
                                        c.execute('UPDATE psc SET ec_id=? WHERE uniref90_id=?', (ec_id, r['uniref90_id']))
                                        logger_psc.info('UPDATE psc SET ec_id=%s WHERE uniref90_id=%s', ec_id, r['uniref90_id'])
                                        update = True
                                    if(update):
                                        psc_annotated += 1
                            if((sp_processed % 10000) == 0):
                                conn.commit()
                                print("\t... %d" % sp_processed)
                        except sqlite3.Error as e:
                            print("SQL ERROR: could not set annotation for sw-id=%s" % acc)
                            print(e)
                elem.clear()  # forstall out of memory errors

print('\n')
print("SwissProt proteins not found: %d" % sp_not_found)
print("Hash protein length collision: %d" % sp_lenth_mm)

print("UPSs annotated: %d" % ups_annotated)
logger_ups.debug('summary: UPS annotated=%d', ups_annotated)
print("PSCs annotated: %d" % psc_annotated)
logger_psc.debug('summary: PSC annotated=%d', psc_annotated)