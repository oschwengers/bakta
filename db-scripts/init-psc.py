
import argparse
import logging
import gzip
import sqlite3
from lxml import etree as et
from pathlib import Path

from Bio import SeqIO


parser = argparse.ArgumentParser(
    description='Filter Uniprot\'s UniRef90 XML files to bacterial subsequences and init pc db.'
)
parser.add_argument('--taxonomy', action='store', help='Path to NCBI taxonomy node.dmp file.')
parser.add_argument('--xml', action='store', help='Path to UniRef90 xml file.')
parser.add_argument('--uniparc', action='store', help='Path to UniParc fasta file.')
parser.add_argument('--db', action='store', help='Path to Bakta sqlite3 db file.')
parser.add_argument('--fasta', action='store', help='Path to PSC fasta file.')
args = parser.parse_args()


taxonomy_path = Path(args.taxonomy).resolve()
xml_path = Path(args.xml).resolve()
uniparc_path = Path(args.uniparc).resolve()
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

uniref90_uniparc_ids = {}

print('parse & store PSC information...')
with sqlite3.connect(str(db_path), isolation_level='EXCLUSIVE') as conn:
    conn.execute('PRAGMA page_size = 4096;')
    conn.execute('PRAGMA cache_size = 100000;')
    conn.execute('PRAGMA locking_mode = EXCLUSIVE;')
    conn.execute("PRAGMA mmap_size = %i;" % (20 * 1024 * 1024 * 1024))
    conn.execute('PRAGMA synchronous = OFF;')
    conn.execute('PRAGMA journal_mode = OFF')
    conn.execute('PRAGMA threads = 2;')
    conn.commit()

    with gzip.open(str(xml_path), mode='rb') as fh_xml, fasta_path.open(mode='wt') as fh_fasta:
        i = 0
        for event, elem in et.iterparse(fh_xml, tag='{*}entry'):
            if('Fragment' not in elem.find('./{*}name').text):  # skip protein fragments
                common_tax_id = elem.find('./{*}property[@type="common taxon ID"]')
                common_tax_id = common_tax_id.get('value') if common_tax_id is not None else 1

                rep_member = elem.find('./{*}representativeMember')
                rep_member_dbref = rep_member.find('./{*}dbReference')

                rep_member_organism = rep_member_dbref.find('./{*}property[@type="source organism"]')  # source organism
                rep_member_organism = rep_member_organism.get('value') if rep_member_organism is not None else ''
                
                rep_member_tax_id = rep_member_dbref.find('./{*}property[@type="NCBI taxonomy"]')
                rep_member_tax_id = rep_member_tax_id.get('value') if rep_member_tax_id is not None else 1
                
                if(is_taxon_child(common_tax_id, '2', taxonomy) or is_taxon_child(rep_member_tax_id, '2', taxonomy) or 'phage' in rep_member_organism.lower()):
                    uniref90_id = elem.attrib['id'][9:]  # remove 'UniRef90_' prefix
                    
                    product = rep_member_dbref.find('./{*}property[@type="protein name"]')
                    if(product is not None):
                        product = product.get('value')
                        if(product.lower() == 'hypothetical protein'):
                            product = None
                        elif(product.lower() == 'uncharacterized protein'):
                            product = None
                    
                    uniref50_id = rep_member_dbref.find('./{*}property[@type="UniRef50 ID"]')
                    if(uniref50_id is not None):
                        uniref50_id = uniref50_id.get('value')[9:]  # remove 'UniRef50_' prefix
                    
                    psc = (
                        uniref90_id,
                        uniref50_id,
                        product
                    )
                    conn.execute('INSERT INTO psc (uniref90_id, uniref50_id, product) VALUES (?,?,?)', psc)
                    log.info('INSERT INTO psc (uniref90_id, uniref50_id, product) VALUES (%s,%s,%s)', *psc)
                    
                    # lookup seed sequence
                    is_seed = rep_member_dbref.find('./{*}property[@type="isSeed"]')
                    if(is_seed is not None):  # representative is seed sequence
                        seq = rep_member.find('./{*}sequence').text.upper()
                        fh_fasta.write(">%s\n%s\n" % (uniref90_id, seq))
                        seed_db_type = rep_member_dbref.get('type')
                        seed_db_id = rep_member_dbref.get('id')
                        log.info('write seed: uniref90-id=%s, type=%s, id=%s, length=%s', uniref90_id, seed_db_type, seed_db_id, len(seq))
                    else:  # search for seed member
                        for member_dbref in elem.findall('./{*}member/{*}dbReference'):
                            if(member_dbref.find('./{*}property[@type="isSeed"]') is not None):
                                uniparc_id = member_dbref.find('./{*}property[@type="UniParc ID"]')  # search for UniParc annotation
                                if(uniparc_id is not None):  # use UniParc ID
                                    seed_db_type = 'UniParc ID'
                                    seed_db_id = uniparc_id.attrib['value']
                                else:  # use DBRef type (either UniParc or UniProtKB)
                                    seed_db_type = member_dbref.get('type')
                                    seed_db_id = member_dbref.get('id')
                                
                                if(seed_db_type == 'UniParc ID'):
                                    uniref90_uniparc_ids[seed_db_id] = uniref90_id
                                else:
                                    print("detected additional seed type! UniRef90-id=%s, seed-type=%s, seed-id=%s" % (uniref90_id, seed_db_type, seed_db_id))
                                break
                            member_dbref.clear()
                    
                    i += 1
                    if((i % 1000000) == 0):
                       conn.commit()
                       print("\t... %i" % i)
                rep_member.clear()
                rep_member_dbref.clear()
            elem.clear()  # forstall out of memory errors
    conn.commit()
print("\tparsed PSC: %d" % i)
log.debug('summary: # PSC=%d', i)


print("UniParc (%i)..." % len(uniref90_uniparc_ids))
log.debug('lookup non-representative UniParc seed sequences: %s', len(uniref90_uniparc_ids))
i = 0
with gzip.open(str(uniparc_path), mode='rt') as fh_uniparc, fasta_path.open(mode='at') as fh_fasta:
    for record in SeqIO.parse(fh_uniparc, 'fasta'):
        uniref90_id = uniref90_uniparc_ids.get(record.id, None)
        if(uniref90_id):
            fh_fasta.write(">%s\n%s\n" % (uniref90_id, str(record.seq).upper()))
            log.info('write seed: uniref90-id=%s, type=UniParc, id=%s, length=%s', uniref90_id, record.id, len(record.seq))
            uniref90_uniparc_ids.pop(record.id)
            i += 1
print("\twritten UniParc seed sequences: %i" % i)
log.debug('written UniParc seed sequences: %i', i)


print("\nsuccessfully initialized PSC table!")
