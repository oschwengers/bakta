import argparse
import logging
import hashlib
import gzip
import sqlite3
from pathlib import Path
from lxml import etree as et


parser = argparse.ArgumentParser(
    description='Filter Uniprot\'s UniRef100 XML files to bacterial subsequences and create initial ips db.'
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
log_ups = logging.getLogger('UPS')
log_ips = logging.getLogger('IPS')
log_psc = logging.getLogger('PSC')


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
print("\tstored tax ids: %i" % len(taxonomy))

sp_processed = 0
sp_not_found = 0

ips_annotated = 0
psc_annotated = 0

print('parse & store SwissProt information...')
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
    with gzip.open(str(xml_path), mode="rb") as fh:
        ups_entries = []
        i = 0
        for event, elem in et.iterparse(fh, tag='{*}entry'):
            sp_processed += 1
            acc = elem.find('./{*}accession').text
            elem_org = elem.find('./{*}organism')
            tax_property = elem_org.find('./{*}dbReference[@type="NCBI Taxonomy"]')
            tax_id = tax_property.get('id') if tax_property is not None else '1'
            org_name = elem_org.find('./{*}name[@type=scientific]')
            org_name = org_name.get('value').lower() if org_name is not None else ''
            if(is_taxon_child(tax_id, '2', taxonomy) or 'phage' in org_name):
                seq = elem.find('./{*}sequence').text.upper()
                seq_hash = hashlib.md5(seq.encode())
                seq_hash_hexdigest = seq_hash.hexdigest()
                rec_ups = conn.execute('SELECT * FROM ups WHERE hash=?', (seq_hash.digest(),)).fetchone()
                if(rec_ups is None):
                    sp_not_found += 1
                    log_ups.debug('no hash hit: SwissProt-id=%s, hash=%s', acc, seq_hash_hexdigest)
                    continue
                assert rec_ups['length'] == len(seq), "Detected SwissProt / UPS length collision! hash=%s, SwissProt-id=%s, UniParc-id=%s, SwissProt-length=%s, db-length=%s" % (seq_hash_hexdigest, acc, rec_ups['uniparc_id'], len(seq), rec_ups['length'])
                
                if(rec_ups['uniref100_id'] is None):
                    log_ups.debug('no UniRef100-id: SwissProt-id=%s, hash=%s', acc, seq_hash_hexdigest)
                    continue
                uniref100_id = rec_ups['uniref100_id']
                
                rec_ips = conn.execute('SELECT * FROM ips WHERE uniref100_id=?', (uniref100_id,)).fetchone()
                if(rec_ips is None):
                    log_ups.debug('no UniRef100 hit: SwissProt-id=%s, hash=%s', acc, seq_hash_hexdigest)
                    continue
                # parse SwissProt annotations
                gene = elem.find('./{*}gene/{*}name[@type="primary"]')
                if(gene is not None):
                    gene = gene.text
                product = elem.find('./{*}protein/{*}recommendedName/{*}fullName')
                if(product is not None):
                    product = product.text
                ec_ids = []
                for ec_id in elem.findall('./{*}protein/{*}recommendedName/{*}ecNumber'):
                    ec_ids.append(ec_id.text)
                go_terms = []
                for go_term in elem.findall('./{*}dbReference[@type="GO"]'):
                    go_terms.append(go_term.get('id'))[3:]
                uniref90_id = rec_ips['uniref90_id']
                if(uniref90_id is None):  # no PSC found -> annotate IPS
                    update = False
                    if(gene):
                        conn.execute('UPDATE ips SET gene=? WHERE uniref100=?', (gene, uniref100_id))
                        log_ips.info('UPDATE ips SET gene=%s WHERE uniref100_id=%s', gene, uniref100_id)
                        update = True
                    if(product):
                        conn.execute('UPDATE ips SET product=? WHERE uniref100_id=?', (product, uniref100_id))
                        log_ips.info('UPDATE ips SET product=%s WHERE uniref100_id=%s', product, uniref100_id)
                        update = True
                    if(len(ec_ids) > 0):
                        ec_list = ','.join(ec_ids)
                        conn.execute('UPDATE ips SET ec_ids=? WHERE uniref100_id=?', (ec_list, uniref100_id))
                        log_ips.info('UPDATE ips SET ec_ids=%s WHERE uniref100_id=%s', ec_list, uniref100_id)
                        update = True
                    if(len(go_terms) > 0):
                        go_list = ','.join(go_terms)
                        conn.execute('UPDATE ips SET go_ids=? WHERE uniref100_id=?', (go_list, uniref100_id))
                        log_ips.info('UPDATE ips SET go_ids=%s WHERE uniref100_id=%s', go_list, uniref100_id)
                        update = True
                    if(update):
                        ips_annotated += 1
                else:  # PSC found -> annotate PSC
                    update = False
                    if(gene):
                        conn.execute('UPDATE psc SET gene=? WHERE uniref90_id=?', (gene, uniref90_id))
                        log_psc.info('UPDATE psc SET gene=%s WHERE uniref90_id=%s', gene, uniref90_id)
                        update = True
                    if(product):
                        conn.execute('UPDATE psc SET product=? WHERE uniref90_id=?', (product, uniref90_id))
                        log_psc.info('UPDATE psc SET product=%s WHERE uniref90_id=%s', product, uniref90_id)
                        update = True
                    if(ec_id):
                        conn.execute('UPDATE psc SET ec_ids=? WHERE uniref90_id=?', (ec_id, uniref90_id))
                        log_psc.info('UPDATE psc SET ec_ids=%s WHERE uniref90_id=%s', ec_id, uniref90_id)
                        update = True
                    if(len(go_terms) > 0):
                        go_list = ','.join(go_terms)
                        conn.execute('UPDATE psc SET go_ids=? WHERE uniref90_id=?', (go_list, uniref90_id))
                        log_psc.info('UPDATE psc SET go_ids=%s WHERE uniref90_id=%s', go_list, uniref90_id)
                        update = True
                    if(update):
                        psc_annotated += 1
                if((sp_processed % 10000) == 0):
                    conn.commit()
                    print("\t... %i" % sp_processed)
            elem.clear()  # forstall out of memory errors

print('\n')
print("SwissProt proteins not found: %i" % sp_not_found)

print("IPSs annotated: %i" % ips_annotated)
log_ips.debug('summary: IPS annotated=%i', ips_annotated)
print("PSCs annotated: %i" % psc_annotated)
log_psc.debug('summary: PSC annotated=%i', psc_annotated)