
import argparse
import logging
import gzip
import sqlite3
import xml.etree.ElementTree as et
from pathlib import Path

from Bio import SeqIO


parser = argparse.ArgumentParser(
    description='Filter Uniprot\'s UniRef90 XML files to bacterial subsequences and init pc db.'
)
parser.add_argument('--taxonomy', action='store', help='Path to NCBI taxonomy node.dmp file.')
parser.add_argument('--xml', action='store', help='Path to UniRef90 xml file.')
parser.add_argument('--uniprotkb', action='store', help='Path to UniProt KB fasta file.')
parser.add_argument('--uniparc', action='store', help='Path to UniParc fasta file.')
parser.add_argument('--db', action='store', help='Path to Bakta sqlite3 db file.')
parser.add_argument('--fasta', action='store', help='Path to PSC fasta file.')
args = parser.parse_args()


taxonomy_path = Path(args.taxonomy).resolve()
xml_path = Path(args.xml).resolve()
uniprotkb_path = Path(args.uniprotkb).resolve()
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

uniref90_uniprotkb_ids = {}
uniref90_uniparc_ids = {}

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
                    rep_member = elem.find("./{%s}representativeMember" % ns)
                    rep_member_dbref = rep_member.find("./{%s}dbReference" % ns)
                    try:
                        rep_member_tax_id = rep_member_dbref.find("./{%s}property[@type='NCBI taxonomy']" % ns).get('value')
                    except:
                        rep_member_tax_id = 1
                    
                    if(is_taxon_child(common_tax_id, '2', taxonomy) or is_taxon_child(rep_member_tax_id, '2', taxonomy)):
                        uniref90_id = elem.attrib['id'][9:]  # remove 'UniRef90_' prefix
                        try:
                            product = rep_member_dbref.find("./{%s}property[@type='protein name']" % ns).get('value')
                            if(product.lower() == 'hypothetical protein'):
                                product = None
                            elif(product.lower() == 'uncharacterized protein'):
                                product = None
                        except:
                            product = None
                        try:
                            uniref50_id = rep_member_dbref.find("./{%s}property[@type='UniRef50 ID']" % ns).get('value')
                            uniref50_id = uniref50_id[9:]  # remove 'UniRef50_' prefix
                        except:
                            uniref50_id = None
                        
                        is_seed = rep_member_dbref.find("./{%s}property[@type='isSeed']" % ns)
                        if(is_seed is not None):  # representative is seed sequence
                            seq = rep_member.find("./{%s}sequence" % ns).text.upper()
                            fh_fasta.write(">%s\n%s\n" % (uniref90_id, seq))
                            seed_db_type = rep_member_dbref.get('type')
                            seed_db_id = rep_member_dbref.get('id')
                        else:  # search for seed member
                            for member_dbref in elem.findall("./{%s}member/{%s}dbReference" % (ns, ns)):
                                if(member_dbref.find("./{%s}property[@type='isSeed']" % ns) is not None):
                                    seed_db_type = member_dbref.get('type')
                                    seed_db_id = member_dbref.get('id')
                                    if(seed_db_type == 'UniProtKB ID'):
                                        uniref90_uniprotkb_ids[seed_db_id] = uniref90_id
                                    elif(seed_db_type == 'UniParc ID'):
                                        uniref90_uniparc_ids[seed_db_id] = uniref90_id
                                    break
                                member_dbref.clear()
                        psc = (
                            uniref90_id,
                            uniref50_id,  # UniRef50
                            product
                        )
                        log.debug(
                            'uniref90-id=%s, uniref50-id=%s, seed-type=%s, seed-id=%s, product=%s', 
                            uniref90_id, uniref50_id, seed_db_type, seed_db_id, product
                        )
                        conn.execute("INSERT INTO psc (uniref90_id, uniref50_id, product) VALUES (?,?,?)", psc)
                        log.info('INSERT INTO psc (uniref90_id, uniref50_id, product) VALUES (%s,%s,%s)', *psc)
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


print('Lookup non-representative seed sequences in:')
print("\tUniProtKb (%i)..." % len(uniref90_uniprotkb_ids))
i = 0
with gzip.open(str(uniprotkb_path), mode='rt') as fh_uniprotkb, fasta_path.open(mode='at') as fh_fasta:
    for record in SeqIO.parse(fh_uniprotkb, 'fasta'):
        uniprotkb_id = record.id.split('|')[2]  # >tr|H9BV67|H9BV67_9EURY
        uniref90_id = uniref90_uniprotkb_ids.get(uniprotkb_id, None)
        if(uniref90_id):
            log.debug('write seed: UniRef90=%s, UniProtKb=%s', uniref90_id, uniprotkb_id)
            fh_fasta.write(">%s\n%s\n" % (uniref90_id, str(record.seq).upper()))
            uniref90_uniprotkb_ids.pop(uniprotkb_id)
            i += 1
print("\twritten UniProtKB seed sequences: %i" % i)
log.debug('written UniProtKB seed sequences: %i', i)


print("UniParc (%i)..." % len(uniref90_uniparc_ids))
i = 0
with gzip.open(str(uniparc_path), mode='rt') as fh_uniparc, fasta_path.open(mode='at') as fh_fasta:
    for record in SeqIO.parse(fh_uniparc, 'fasta'):
        uniref90_id = uniref90_uniparc_ids.get(record.id, None)
        if(uniref90_id):
            log.debug('write seed: UniRef90=%s, UniParc=%s', uniref90_id, record.id)
            fh_fasta.write(">%s\n%s\n" % (uniref90_id, str(record.seq).upper()))
            uniref90_uniparc_ids.pop(record.id)
            i += 1
print("\twritten UniParc seed sequences: %i" % i)
log.debug('written UniParc seed sequences: %i', i)


print("\nsuccessfully initialized PSC table!")
