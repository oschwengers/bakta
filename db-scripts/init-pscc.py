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
parser.add_argument('--uniref50', action='store', help='Path to UniRef50 xml file.')
parser.add_argument('--uniparc', action='store', help='Path to UniParc fasta file.')
parser.add_argument('--db', action='store', help='Path to Bakta sqlite3 db file.')
parser.add_argument('--pscc', action='store', help='Path to PSCC fasta file.')
parser.add_argument('--pscc_sorf', action='store', help='Path to sORF PSCC fasta file.')
args = parser.parse_args()

PSCC_MIN_MEMBER_COUNT = 10
MAX_SORF_LENGTH = 30

DISCARDED_PRODUCTS = [
    'hypothetical protein',
    'hypothetical conserved protein',
    'uncharacterized protein',
    'hypothetical membrane protein',
    'hypothetical cytosolic Protein'
]

taxonomy_path = Path(args.taxonomy).resolve()
uniref50_path = Path(args.uniref50).resolve()
uniparc_path = Path(args.uniparc).resolve()
db_path = Path(args.db)
pscc_path = Path(args.pscc)
pscc_sorf_path = Path(args.pscc_sorf)


logging.basicConfig(
    filename='bakta.db.log',
    filemode='a',
    format='%(name)s - UniRef50 - %(levelname)s - %(message)s',
    level=logging.DEBUG
)
log_pscc = logging.getLogger('PSCC')
log_sorf = logging.getLogger('PSCC_SORF')


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


pscc_seqs = 0
pscc_sorf_seqs = 0
pscc_total = 0
uniref50_uniparc_ids = {}
print('parse & store PSCC information...')
with sqlite3.connect(str(db_path), isolation_level='EXCLUSIVE') as conn, xopen(str(uniref50_path), mode='rb') as fh_xml, pscc_path.open(mode='wt') as fh_fasta_pscc, pscc_sorf_path.open(mode='wt') as fh_fasta_pscc_sorf, alive_bar() as bar:
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
            try:
                member_count = elem.find('./{*}property[@type="member count"]')
                member_count = int(member_count.get('value')) if member_count is not None else 1
                common_tax_id = elem.find('./{*}property[@type="common taxon ID"]')
                common_tax_id = common_tax_id.get('value') if common_tax_id is not None else 1
                rep_member = elem.find('./{*}representativeMember')
                rep_member_dbref = rep_member.find('./{*}dbReference')
                rep_member_organism = rep_member_dbref.find('./{*}property[@type="source organism"]')  # source organism
                rep_member_organism = rep_member_organism.get('value') if rep_member_organism is not None else ''
                rep_member_tax_id = rep_member_dbref.find('./{*}property[@type="NCBI taxonomy"]')
                rep_member_tax_id = rep_member_tax_id.get('value') if rep_member_tax_id is not None else 1
                if(is_taxon_child(common_tax_id, '2', taxonomy) or is_taxon_child(rep_member_tax_id, '2', taxonomy) or 'phage' in rep_member_organism.lower()):
                    uniref50_id = elem.attrib['id'][9:]  # remove 'UniRef50_' prefix
                    product = rep_member_dbref.find('./{*}property[@type="protein name"]')
                    if(product is not None):
                        product = product.get('value')
                        if(product.lower() in DISCARDED_PRODUCTS):
                            product = None
                    # print(f"\n{i}")
                    # print(f"id={uniref50_id}")
                    # print(f"member_count={member_count}")
                    # print(f"common_tax_id={common_tax_id}")
                    # print(f"rep_member={rep_member}")
                    # print(f"rep_member_dbref={rep_member_dbref}")
                    # print(f"rep_member_organism={rep_member_organism}")
                    # print(f"rep_member_tax_id={rep_member_tax_id}")
                    # print(f"product={product}")
                    # lookup seed sequence
                    is_seed = rep_member_dbref.find('./{*}property[@type="isSeed"]')
                    if(is_seed is not None):  # representative is seed sequence
                        if(product is None or '(fragment)' not in product.lower()):  # skip protein fragments
                            conn.execute('INSERT INTO pscc (uniref50_id, product) VALUES (?,?)', (uniref50_id, product))
                            log_pscc.info('INSERT INTO pscc (uniref50_id, product) VALUES (%s,%s)', uniref50_id, product)
                            pscc_total += 1
                            if(member_count >= PSCC_MIN_MEMBER_COUNT):
                                seq = rep_member.find('./{*}sequence').text.upper()
                                seed_db_type = rep_member_dbref.get('type')
                                seed_db_id = rep_member_dbref.get('id')
                                if(len(seq) >= MAX_SORF_LENGTH):  # normael CDS
                                    fh_fasta_pscc.write(f'>{uniref50_id}\n{seq}\n')
                                    log_pscc.info('write seed: uniref50-id=%s, type=%s, id=%s, length=%s', uniref50_id, seed_db_type, seed_db_id, len(seq))
                                    pscc_seqs += 1
                                else:  # short ORF
                                    fh_fasta_pscc_sorf.write(f'>{uniref50_id}\n{seq}\n')
                                    log_sorf.info('write seed: uniref50-id=%s, type=%s, id=%s, length=%s', uniref50_id, seed_db_type, seed_db_id, len(seq))
                                    pscc_sorf_seqs += 1
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
                                    conn.execute('INSERT INTO pscc (uniref50_id, product) VALUES (?,?)', (uniref50_id, product))
                                    log_pscc.info('INSERT INTO pscc (uniref50_id, product) VALUES (%s,%s)', uniref50_id, product)
                                    pscc_total += 1
                                    if(member_count >= PSCC_MIN_MEMBER_COUNT):
                                        if(uniparc_id is not None):  # use UniParc ID
                                            seed_db_type = 'UniParc ID'
                                            seed_db_id = uniparc_id.attrib['value']
                                        else:  # use DBRef type (either UniParc or UniProtKB)
                                            seed_db_type = member_dbref.get('type')
                                            seed_db_id = member_dbref.get('id')
                                        if(seed_db_type == 'UniParc ID'):
                                            uniref50_uniparc_ids[seed_db_id] = uniref50_id
                                        else:
                                            print(f'detected additional seed type! UniRef50-id={uniref50_id}, seed-type={seed_db_type}, seed-id={seed_db_id}')
                                break
                            member_dbref.clear()
                    i += 1
                    if((i % 1_000_000) == 0):
                        conn.commit()
            except Exception as e:
                try:
                    uniref50_id = elem.attrib['id']
                except:
                    uniref50_id = '-'
                print(f'Error: {uniref50_id}; {elem}')
                print(e)
            elem.clear()  # forstall out of memory errors
            bar()
    conn.commit()
print(f'parsed PSCC: {pscc_total}')
log_pscc.debug('summary: # PSCC=%i', pscc_total)
print('\n')


print(f'UniParc ({len(uniref50_uniparc_ids)})...')
log_pscc.debug('lookup non-representative UniParc seed sequences: %s', len(uniref50_uniparc_ids))
with xopen(str(uniparc_path), mode='rt') as fh_uniparc, pscc_path.open(mode='at') as fh_fasta_psc, pscc_sorf_path.open(mode='at') as fh_fasta_sorf, alive_bar() as bar:
    for record in SeqIO.parse(fh_uniparc, 'fasta'):
        uniref50_id = uniref50_uniparc_ids.get(record.id, None)
        if(uniref50_id):
            seq = str(record.seq).upper()
            if(len(seq) >= MAX_SORF_LENGTH):  # normael CDS
                fh_fasta_psc.write(f'>{uniref50_id}\n{seq}\n')
                log_pscc.info('write seed: uniref50-id=%s, type=UniParc, id=%s, length=%s', uniref50_id, record.id, len(record.seq))
                pscc_seqs += 1
            else:  # short ORF
                fh_fasta_sorf.write(f'>{uniref50_id}\n{seq}\n')
                log_sorf.info('write seed: uniref50-id=%s, type=UniParc, id=%s, length=%s', uniref50_id, record.id, len(record.seq))
                pscc_sorf_seqs += 1
            uniref50_uniparc_ids.pop(record.id)
        bar()
uniparc_total_seqs = pscc_seqs + pscc_sorf_seqs
print(f'written UniParc seed sequences: {uniparc_total_seqs}')
log_pscc.debug('written UniParc seed sequences: %i', uniparc_total_seqs)

print(f'PSCC normal seqs: {pscc_seqs}')
log_pscc.debug('summary: # PSCC normal=%i', pscc_seqs)
print(f'PSCC sORF seqs: {pscc_sorf_seqs}')
log_pscc.debug('summary: # PSCC sORFs=%i', pscc_sorf_seqs)

print("\nsuccessfully initialized PSCC table!")
