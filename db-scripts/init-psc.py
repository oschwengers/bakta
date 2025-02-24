import argparse
import logging
import sqlite3

from pathlib import Path


from alive_progress import alive_bar
from Bio import SeqIO
from xml.etree import ElementTree as et
from xopen import xopen


parser = argparse.ArgumentParser(
    description='Filter Uniprot\'s UniRef90 XML files to bacterial subsequences and init pc db.'
)
parser.add_argument('--taxonomy', action='store', help='Path to NCBI taxonomy node.dmp file.')
parser.add_argument('--uniref90', action='store', help='Path to UniRef90 xml file.')
parser.add_argument('--uniparc', action='store', help='Path to UniParc fasta file.')
parser.add_argument('--db', action='store', help='Path to Bakta sqlite3 db file.')
parser.add_argument('--psc', action='store', help='Path to PSC fasta file.')
parser.add_argument('--psc_sorf', action='store', help='Path to sORF PSC fasta file.')
args = parser.parse_args()

MAX_SORF_LENGTH = 30

DISCARDED_PRODUCTS = [
    'hypothetical protein',
    'hypothetical conserved protein',
    'uncharacterized protein',
    'hypothetical membrane protein',
    'hypothetical cytosolic Protein'
]

taxonomy_path = Path(args.taxonomy).resolve()
uniref90_path = Path(args.uniref90).resolve()
uniparc_path = Path(args.uniparc).resolve()
db_path = Path(args.db)
psc_path = Path(args.psc)
psc_sorf_path = Path(args.psc_sorf)


logging.basicConfig(
    filename='bakta.db.log',
    filemode='a',
    format='%(name)s - UniRef90 - %(levelname)s - %(message)s',
    level=logging.DEBUG
)
log_psc = logging.getLogger('PSC')
log_sorf = logging.getLogger('PSC_SORF')


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
print(f'\tstored tax ids: {len(taxonomy)}')

uniref90_uniparc_ids = {}


psc_total = 0
psc_seqs = 0
psc_sorf_seqs = 0
print('parse & store PSC information...')
with sqlite3.connect(str(db_path), isolation_level='EXCLUSIVE') as conn, xopen(str(uniref90_path), mode='rb') as fh_xml, psc_path.open(mode='wt') as fh_fasta_psc, psc_sorf_path.open(mode='wt') as fh_fasta_psc_sorf, alive_bar(enrich_print=False) as bar:
    conn.execute('PRAGMA page_size = 4096;')
    conn.execute('PRAGMA cache_size = 100000;')
    conn.execute('PRAGMA locking_mode = EXCLUSIVE;')
    conn.execute(f'PRAGMA mmap_size = {20 * 1024 * 1024 * 1024};')
    conn.execute('PRAGMA synchronous = OFF;')
    conn.execute('PRAGMA journal_mode = OFF')
    conn.execute('PRAGMA threads = 2;')
    conn.commit()

    i = 0
    for event, elem in et.iterparse(fh_xml, events=('end',)):
        if elem.tag == '{http://uniprot.org/uniref}entry':
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
                    if(product.lower() in DISCARDED_PRODUCTS):
                        product = None

                uniref50_id = rep_member_dbref.find('./{*}property[@type="UniRef50 ID"]')
                if(uniref50_id is not None):
                    uniref50_id = uniref50_id.get('value')[9:]  # remove 'UniRef50_' prefix

                # lookup seed sequence
                is_seed = rep_member_dbref.find('./{*}property[@type="isSeed"]')
                if(is_seed is not None):  # representative is seed sequence
                    if(product is None or '(fragment)' not in product.lower()):  # skip protein fragments
                        conn.execute('INSERT INTO psc (uniref90_id, uniref50_id, product) VALUES (?,?,?)', (uniref90_id, uniref50_id, product))
                        log_psc.info('INSERT INTO psc (uniref90_id, uniref50_id, product) VALUES (%s,%s,%s)', uniref90_id, uniref50_id, product)
                        psc_total += 1
                        seq = rep_member.find('./{*}sequence').text.upper()
                        seed_db_type = rep_member_dbref.get('type')
                        seed_db_id = rep_member_dbref.get('id')
                        if(len(seq) >= MAX_SORF_LENGTH):  # normael CDS
                            fh_fasta_psc.write(f'>{uniref90_id}\n{seq}\n')
                            log_psc.info('write seed: uniref90-id=%s, type=%s, id=%s, length=%s', uniref90_id, seed_db_type, seed_db_id, len(seq))
                            psc_seqs += 1
                        else:  # short ORF
                            fh_fasta_psc_sorf.write(f'>{uniref90_id}\n{seq}\n')
                            log_sorf.info('write seed: uniref90-id=%s, type=%s, id=%s, length=%s', uniref90_id, seed_db_type, seed_db_id, len(seq))
                            psc_sorf_seqs += 1
                else:  # search for seed member
                    for member_dbref in elem.findall('./{*}member/{*}dbReference'):
                        if(member_dbref.find('./{*}property[@type="isSeed"]') is not None):
                            uniparc_id = member_dbref.find('./{*}property[@type="UniParc ID"]')  # search for UniParc annotation
                            seed_product = member_dbref.find('./{*}property[@type="protein name"]')
                            seed_product = seed_product.get('value') if seed_product is not None else None
                            if(seed_product is not None and '(fragment)' not in seed_product.lower()):  # seed sequence is not fragmented as well
                                if(product is not None):
                                    product = product.replace(' (Fragment)', '')  # discard fragmented tag from product for unfragmented seeds
                                    if(product.lower() in DISCARDED_PRODUCTS):
                                        product = None
                                conn.execute('INSERT INTO psc (uniref90_id, uniref50_id, product) VALUES (?,?,?)', (uniref90_id, uniref50_id, product))
                                log_psc.info('INSERT INTO psc (uniref90_id, uniref50_id, product) VALUES (%s,%s,%s)', uniref90_id, uniref50_id, product)
                                psc_total += 1
                                if(uniparc_id is not None):  # use UniParc ID
                                    seed_db_type = 'UniParc ID'
                                    seed_db_id = uniparc_id.attrib['value']
                                else:  # use DBRef type (either UniParc or UniProtKB)
                                    seed_db_type = member_dbref.get('type')
                                    seed_db_id = member_dbref.get('id')
                                if(seed_db_type == 'UniParc ID'):
                                    uniref90_uniparc_ids[seed_db_id] = uniref90_id
                                else:
                                    print(f'detected additional seed type! UniRef90-id={uniref90_id}, seed-type={seed_db_type}, seed-id={seed_db_id}')
                                break
                        member_dbref.clear()
            i += 1
            if((i % 1_000_000) == 0):
                conn.commit()
            elem.clear()  # forstall out of memory errors
            bar()
    conn.commit()
print(f'parsed PSC: {psc_total}')
log_psc.debug('summary: # PSC=%i', psc_total)
print('\n')


print(f'UniParc ({len(uniref90_uniparc_ids)})...')
log_psc.debug('lookup non-representative UniParc seed sequences: %s', len(uniref90_uniparc_ids))
with xopen(str(uniparc_path), mode='rt') as fh_uniparc, psc_path.open(mode='at') as fh_fasta_psc, psc_sorf_path.open(mode='at') as fh_fasta_sorf, alive_bar(enrich_print=False) as bar:
    for record in SeqIO.parse(fh_uniparc, 'fasta'):
        uniref90_id = uniref90_uniparc_ids.get(record.id, None)
        if(uniref90_id):
            seq = str(record.seq).upper()
            if(len(seq) >= MAX_SORF_LENGTH):  # normael CDS
                fh_fasta_psc.write(f'>{uniref90_id}\n{seq}\n')
                log_psc.info('write seed: uniref90-id=%s, type=UniParc, id=%s, length=%s', uniref90_id, record.id, len(record.seq))
                psc_seqs += 1
            else:  # short ORF
                fh_fasta_sorf.write(f'>{uniref90_id}\n{seq}\n')
                log_sorf.info('write seed: uniref90-id=%s, type=UniParc, id=%s, length=%s', uniref90_id, record.id, len(record.seq))
                psc_sorf_seqs += 1
            uniref90_uniparc_ids.pop(record.id)
        bar()
uniparc_total_seqs = psc_seqs + psc_sorf_seqs
print(f'written UniParc seed sequences: {uniparc_total_seqs}')
log_psc.debug('written UniParc seed sequences: %i', uniparc_total_seqs)

print(f'PSC normal seqs: {psc_seqs}')
log_psc.debug('summary: # PSC normal=%i', psc_seqs)
print(f'PSC sORF seqs: {psc_sorf_seqs}')
log_psc.debug('summary: # PSC sORFs=%i', psc_sorf_seqs)

print("\nsuccessfully initialized PSC table!")
