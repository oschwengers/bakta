import argparse
import logging
from pathlib import Path

from Bio import SeqIO


parser = argparse.ArgumentParser(
    description='Setup NCBI BlastRules expert system.'
)
parser.add_argument('--expert-sequence', action='store', dest='expert_sequence', required=True, help='Path to Bakta expert sequence file.')
parser.add_argument('--ncbi-blastrule-tsv', action='store', dest='blastrule_tsv', required=True, help='Path to NCBI blastrule file.')
parser.add_argument('--ncbi-blastrule-proteins', action='store', dest='blastrule_proteins', required=True, help='Path to NCBI blastrule proteins file.')
args = parser.parse_args()


expert_sequence_path = Path(args.expert_sequence).resolve()
blastrule_tsv_path = Path(args.blastrule_tsv).resolve()
blastrule_proteins_path = Path(args.blastrule_proteins).resolve()


logging.basicConfig(
    filename='bakta.db.log',
    filemode='a',
    format='%(name)s - EXPERT-SEQ - %(levelname)s - %(message)s',
    level=logging.DEBUG
)
log = logging.getLogger('BlastRules')


print('import NCBI BlastRules information...')
blast_rules = {}
with blastrule_tsv_path.open(encoding='windows-1252') as fh:
    for line in fh:
        (rule_acc, rule_type, product, gbk_acc, refseq_accs, id_threshold, subject_cov_threshold, query_cov_threshold, gene, pm_acc, ec, about) = line.split('\t')
        if(rule_type in ['BlastRuleEquivalog', 'BlastRuleException', 'BlastRuleIS', 'BlastRuleExact']):
            br = {
                'id': rule_acc,
                'gene': gene,
                'product': product,
                'refseq_accs': refseq_accs.split(','),
                'identity': float(id_threshold),
                'query_cov': float(query_cov_threshold),
                'subject_cov': float(subject_cov_threshold),
                'dbxref': []
            }
            if(ec != ''):
                for ec_entry in ec.split(','):
                    br['dbxref'].append(f"EC:{ec_entry}")
            if(rule_type == 'BlastRuleEquivalog'):
                br['rank'] = 70
            elif(rule_type == 'BlastRuleException'):
                br['rank'] = 80
            elif(rule_type == 'BlastRuleIS'):
                br['rank'] = 90
            elif(rule_type == 'BlastRuleExact'):
                br['rank'] = 99
            else:
                br['rank'] = 1
            for refseq_acc in br['refseq_accs']:
                blast_rules[refseq_acc] = br
print(f'\tstored BlastRule ids: {len(blast_rules)}')

blast_rule_seqs = 0
print('import NCBI BlastRule proteins...')
with blastrule_proteins_path.open() as fh_in, expert_sequence_path.open('w+') as fh_out:
    for record in SeqIO.parse(fh_in, 'fasta'):
        br_id = record.id
        br_seq = str(record.seq).upper()
        refseq_id = None
        for entry in br_id.split('|'):
            if(entry.startswith('WP_')):
                refseq_id = entry.split('.')[0]
        br = blast_rules.get(refseq_id, None)
        if(br is not None):
            fh_out.write(f">{refseq_id} BlastRules~~~{br['rank']}~~~{br['identity']}~~~{br['query_cov']}~~~{br['subject_cov']}~~~{br['gene']}~~~{br['product']}~~~{','.join(br['dbxref'])}\n")
            fh_out.write(f"{br_seq}\n")
            log.info('write seq: RefSeq-id=%s, rank=%i, id=%f, q-cov=%f, s-cov=%f, gene=%s, product=%s, dbxrefs=%s', refseq_id, br['rank'], br['identity'], br['query_cov'], br['subject_cov'], br['gene'], br['product'], ','.join(br['dbxref']) )
            blast_rule_seqs += 1
print(f'\tstored BlastRule sequences: {blast_rule_seqs}')
log.debug('written BlastRule sequences: %i', blast_rule_seqs)


print("\nsuccessfully setup BlastRules expert system!")