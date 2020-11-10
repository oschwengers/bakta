
import logging

import bakta.constants as bc

log = logging.getLogger('TSV')

############################################################################
# Inference terms
# 
# tRNA: 'tRNAscan'
# tmRNA: 'aragorn'
# rRNA: 'Rfam:%s' % subject_id
# ncRNA genes: 'Rfam:%s' % subject_id
# ncRNA regions: 'Rfam:%s' % subject_id
# CDS hyp: 'Prodigal'
# CDS PSC: 'UniProtKB:%s' % psc[DB_PSC_COL_UNIREF90]
# CDS IPS: 'UniProtKB:%s' % psc[DB_PSC_COL_UNIREF100]
# CDS sORF: 'similar to AA sequence:UniProtKB:%s' % psc[DB_PSC_COL_UNIREF100]
############################################################################

def write_tsv(contigs, features_by_contig, tsv_path):
    """Export features in TSV format."""
    with tsv_path.open('w') as fh:
        for contig in contigs:
            for feat in features_by_contig[contig['id']]:
                feat_type = feat['type']
                if(feat['type'] == bc.FEATURE_GAP):
                    feat_type = bc.INSDC_FEATURE_ASSEMBLY_GAP if feat['length'] >= 100 else bc.INSDC_FEATURE_GAP
                fh.write('\t'.join([feat['contig'], feat_type, str(feat['start']), str(feat['stop']), feat['strand'], feat.get('gene', ''), feat.get('product', '')]))
                fh.write('\n')
    return

