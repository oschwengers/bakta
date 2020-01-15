
import argparse
import gzip
import pickle
import xml.etree.ElementTree as et
from pathlib import Path


parser = argparse.ArgumentParser(
    description='Filter Uniprot\'s UniRef90 XML files to bacterial subsequences and init pc db.'
)
parser.add_argument('--taxonomy', action='store', help='Path to NCBI taxonomy node.dmp file.')
parser.add_argument('--xml', action='store', help='Path to UniRef xml file.')
parser.add_argument('--db', action='store', help='Path to Pickle db file.')
args = parser.parse_args()


def is_taxon_child(child, LCA, taxonomy):
    parent = taxonomy.get(child, None)
    while(parent is not None and parent != '1'):
        if(parent == LCA):
            return True
        else:
            parent = taxonomy.get(parent, None)
    return False


taxonomy = {}
taxonomy_path = Path(args.taxonomy).resolve()
with taxonomy_path.open() as fh:
    for line in fh:
        cols = line.split('\t|\t', maxsplit=2)
        taxonomy[cols[0]] = cols[1]


unique_seqs = {}
unique_clusters = {}
ns = 'http://uniprot.org/uniref'
xml_path = Path(args.xml).resolve()
with gzip.open(str(xml_path), mode="rt") as fh:
    for event, elem in et.iterparse(fh):
        # print(elem)
        if(elem.tag == "{%s}entry" % ns):
            # print(elem)
            if('Fragment' in elem.find("./{%s}name" % ns).text):  # skip protein fragments
                continue
            tax_id = elem.find("./{%s}property[@type='common taxon ID']" % ns).attrib['value']
            # print("tax-id=%s" % tax_id)
            if is_taxon_child(tax_id, '2', taxonomy):
                # print("tax-id=%s" % tax_id)
                uniref90_id = elem.attrib['id']
                if(uniref90_id in unique_seqs):
                    raise Exception('duplicated hashes! hash=%s, seq-a-id: %s, seq-b-id: %s' % (hash, unique_seqs[hash].uniref100_id), uniref100_id)
                else:
                    rep_member = elem.find("./{%s}representativeMember/{%s}dbReference" % (ns, ns))
                    try:
                        uniparc_id = rep_member.find("./{%s}property[@type='UniParc ID']" % ns).attrib['value']
                    except Exception:
                        uniparc_id = ''
                    uniref100_id = rep_member.find("./{%s}property[@type='UniRef100 ID']" % ns).attrib['value']
                    prot_name = rep_member.find("./{%s}property[@type='protein name']" % ns).attrib['value']
                    seq = elem.find("./{%s}representativeMember/{%s}sequence" % (ns, ns)).text.upper()
                    # cluster = Cluster(uniref90_id, len(seq), seq, prot_name, uniref100_id, uniparc_id)
                    cluster = {
                        'id': uniref90_id,
                        'length': len(seq),
                        'product': prot_name,
                        'uniref100': uniref100_id,
                        'uniparc': uniparc_id
                    }
                    unique_clusters[uniref90_id] = cluster
                    unique_seqs[uniref90_id] = seq
                    # print(seq)


db_path = Path(args.db)
with db_path.open('wb') as fh:
    pickle.dump(unique_clusters, fh, pickle.HIGHEST_PROTOCOL)


with open('pcs.faa', 'w') as fh_faa:
    for k, v in unique_clusters.items():
        fh_faa.write(">%s\n%s\n" % (k, v))


# with open('pcs.tsv', 'w') as fh_tsv, open('pc.faa', 'w') as fh_faa:
#     for k, v in unique_seqs.items():
#         # print id    length  product
#         fh_tsv.write(
#             "%s\t%s\t%s\n"
#             % (k[9:], v.length, v.product)
#         )
#         fh_faa.write(">%s\n%s\n" % (v.id, v.seq))
