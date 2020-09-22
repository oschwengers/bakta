
import argparse
import logging
import hashlib
import gzip
import sqlite3
from pathlib import Path
from lxml import etree as et

from Bio import SeqIO


parser = argparse.ArgumentParser(
    description='Filter Uniprot\'s UniRef100 XML files to bacterial/UPSage sequences and create initial UPS & IPS db.'
)
parser.add_argument('--taxonomy', action='store', help='Path to NCBI taxonomy node.dmp file.')
parser.add_argument('--xml', action='store', help='Path to UniRef xml file.')
parser.add_argument('--uniparc', action='store', help='Path to UniParc fasta file.')
parser.add_argument('--db', action='store', help='Path to Bakta sqlite3 db file.')
args = parser.parse_args()


taxonomy_path = Path(args.taxonomy).resolve()
xml_path = Path(args.xml).resolve()
uniparc_path = Path(args.uniparc).resolve()
db_path = Path(args.db)


logging.basicConfig(
    filename='bakta.db.log',
    filemode='a',
    format='%(name)s - UniRef100 - %(levelname)s - %(message)s',
    level=logging.DEBUG
)
log_ups = logging.getLogger('UPS')
log_ips = logging.getLogger('IPS')


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

uniparc_to_uniref100 = {}

print('parse & store IPS information...')
seq_hashes = set()  # unique protein sequences
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

    i = 0
    i_all = 0
    with gzip.open(str(xml_path), mode='rb') as fh:
        for event, elem in et.iterparse(fh, tag='{*}entry'):
            if('Fragment' not in elem.find('./{*}name').text):  # skip protein fragments
                common_tax_id = elem.find('./{*}property[@type="common taxon ID"]')
                common_tax_id = common_tax_id.get('value') if common_tax_id is not None else 1
                
                rep_member_dbref = elem.find('./{*}representativeMember/{*}dbReference')
                rep_member_organism = rep_member_dbref.find('./{*}property[@type="source organism"]')  # source organism
                rep_member_organism = rep_member_organism.get('value') if rep_member_organism is not None else ''
                
                rep_member_tax_id = rep_member_dbref.find('./{*}property[@type="NCBI taxonomy"]')
                rep_member_tax_id = rep_member_tax_id.get('value') if rep_member_tax_id is not None else 1
                
                # filter for bacterial or phage protein sequences
                if(is_taxon_child(common_tax_id, '2', taxonomy) or is_taxon_child(rep_member_tax_id, '2', taxonomy) or 'phage' in rep_member_organism.lower()):
                    uniref100_id = elem.attrib['id'][10:]
                    seq_representative = elem.find('./{*}representativeMember/{*}sequence')
                    seq = seq_representative.text.upper()
                    length = len(seq)
                    seq_hash = hashlib.md5(seq.encode())
                    seq_hash_hexdigest = seq_hash.hexdigest()
                    assert seq_hash_hexdigest not in seq_hashes, "duplicated hashes! hash=%s, UniRef100-id: %s" % (seq_hash_hexdigest, uniref100_id)
                    seq_hashes.add(seq_hash_hexdigest)
                    uniparc_id = rep_member_dbref.find('./{*}property[@type="UniParc ID"]')
                    if(uniparc_id is not None):
                        uniparc_id = uniparc_id.attrib['value'][3:]
                    conn.execute('INSERT INTO ups (hash, length, uniparc_id, uniref100_id) VALUES (?,?,?,?)', (seq_hash.digest(), length, uniparc_id, uniref100_id))
                    log_ups.info('INSERT INTO ups (hash, length, uniparc_id, uniref100_id) VALUES (%s,%s,%s,%s)', seq_hash_hexdigest, length, uniparc_id, uniref100_id)
                    product = rep_member_dbref.find('./{*}property[@type="protein name"]')
                    if(product is not None):
                        product = product.attrib['value']
                        tmp = product.lower()
                        if(tmp == 'hypothetical protein'):
                            product = None
                        elif(tmp == 'Uncharacterized protein'):
                            product = None
                    
                    uniref90_id = rep_member_dbref.find('./{*}property[@type="UniRef90 ID"]')
                    if(uniref90_id is not None):
                        uniref90_id = uniref90_id.attrib['value'][9:]
                        product = None  # skip IPS product and use the PSC/UniRef90 product annotation
                    
                    ips = (
                        uniref100_id,
                        uniref90_id,
                        product
                    )
                    log_ips.debug('uniref100-id=%s, uniref90-id=%s, product=%s', uniref100_id, uniref90_id, product)
                    conn.execute('INSERT INTO ips (uniref100_id, uniref90_id, product) VALUES (?,?,?)', ips)
                    log_ips.info('INSERT INTO ips (uniref100_id, uniref90_id, product) VALUES (%s,%s,%s)', *ips)
                    i += 1
                    # store UPS members of this IPS and store IDs for later hashing
                    for member_dbref in elem.findall('./{*}member/{*}dbReference'):
                        uniparc_id = member_dbref.find('./{*}property[@type="UniParc ID"]')  # search for UniParc annotation
                        if(uniparc_id is not None):  # use UniParc ID
                            seed_db_type = 'UniParc ID'
                            seed_db_id = uniparc_id.attrib['value']
                        else:  # use DBRef type (either UniParc or UniProtKB)
                            seed_db_type = member_dbref.get('type')
                            seed_db_id = member_dbref.get('id')

                        if(seed_db_type == 'UniParc ID'):
                            if(seed_db_id not in uniparc_to_uniref100):
                                uniparc_to_uniref100[seed_db_id] = uniref100_id
                                log_ups.debug('store member: UniParc-id=%s, UniRef100-id=%s', seed_db_id, uniref100_id)
                        else:
                            print("detected additional seed type! UniRef100-id=%s, seed-type=%s, seed-id=%s" % (uniref100_id, seed_db_type, seed_db_id))
                        member_dbref.clear()  # forstall out of memory errors
                    if((i % 100000) == 0):
                        conn.commit()
            elem.clear()  # forstall out of memory errors
            i_all += 1
            if((i_all % 1000000) == 0):
                print("\t... %i" % i_all)
    conn.commit()
    print("\tstored representative IPS: %i" % i)
    log_ips.debug('representative IPS: # IPS=%i', i)


    print("UniParc (%i)..." % len(uniparc_to_uniref100))
    log_ups.debug('lookup non-representative UniParc member sequences: %s', len(uniparc_to_uniref100))
    i = 0
    i_all = 0
    with gzip.open(str(uniparc_path), mode='rt') as fh_uniparc:
        for record in SeqIO.parse(fh_uniparc, 'fasta'):
            uniparc_id = record.id
            uniref100_id = uniparc_to_uniref100.get(uniparc_id, None)
            if(uniref100_id):
                log_ups.debug('member: UniParc=%s, UniRef100=%s, ', uniparc_id, uniref100_id)
                seq = str(record.seq.upper())
                length = len(seq)
                seq_hash = hashlib.md5(seq.encode())
                seq_hash_hexdigest = seq_hash.hexdigest()
                uniparc_id_short = uniparc_id[3:]
                if(seq_hash_hexdigest in seq_hashes):
                    c = conn.cursor()
                    c.execute('SELECT * FROM ups WHERE hash=?', (seq_hash.digest(),))
                    rec = c.fetchone()
                    c.close()
                    assert rec is not None, "Detected hash duplicate without DB entry! hash=%s, UniParc-id=%s, UniRef100-id: %s" % (seq_hash_hexdigest, uniparc_id, uniref100_id)
                    assert rec['uniref100_id'] == uniref100_id, "Detected hash duplicate with different UniRef100 id! hash=%s, UniParc-id=%s, new UniRef100-id=%s, db UniRef100-id=%s" % (seq_hash_hexdigest, uniparc_id, uniref100_id, rec['uniref100_id'])
                    assert rec['length'] == length, "Detected hash duplicate with different seq length! hash=%s, UniParc-id=%s, UniRef100-id=%s, UniParc-length=%s, db-length=%s" % (seq_hash_hexdigest, uniparc_id, uniref100_id, length, rec['length'])
                    # update existing UPS with UniParc id
                    conn.execute('UPDATE ups SET uniparc_id=? WHERE hash=?', (uniparc_id_short, seq_hash.digest()))
                    log_ups.info('UPDATE ups SET uniparc_id=%s WHERE hash=%s', uniparc_id_short, seq_hash_hexdigest)
                else:
                    conn.execute('INSERT INTO ups (hash, length, uniparc_id, uniref100_id) VALUES (?,?,?,?)', (seq_hash.digest(), length, uniparc_id_short, uniref100_id))
                    log_ups.info('INSERT INTO ups (hash, length, uniparc_id, uniref100_id) VALUES (%s,%s,%s,%s)', seq_hash_hexdigest, length, uniparc_id_short, uniref100_id)
                    seq_hashes.add(seq_hash_hexdigest)
                    uniparc_to_uniref100.pop(uniparc_id)
                    i += 1
                    if((i % 100000) == 0):
                        conn.commit()
            i_all += 1
            if((i_all % 1000000) == 0):
                print("\t... %i" % i_all)
    conn.commit()
    print("\twritten UniParc member sequences: %i" % i)
    log_ups.debug('written UniParc member sequences: %i', i)

print("\nsuccessfully initialized UPS & IPS tables!")
