
import logging

import bakta
import bakta.config as cfg
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
    log.info('write tsv output: path=%s', tsv_path)
    
    with tsv_path.open('w') as fh:
        fh.write(f'#Annotated with Bakta v{bakta.__version__}, https://github.com/oschwengers/bakta\n')
        fh.write(f"#Database v{cfg.db_info['major']}.{cfg.db_info['minor']}, https://doi.org/10.5281/zenodo.4247252\n")
        fh.write('#Sequence Id\tType\tStart\tStop\tStrand\tGene\tProduct\tDbXrefs\n')
        for contig in contigs:
            for feat in features_by_contig[contig['id']]:
                feat_type = feat['type']
                if(feat['type'] == bc.FEATURE_GAP):
                    feat_type = bc.INSDC_FEATURE_ASSEMBLY_GAP if feat['length'] >= 100 else bc.INSDC_FEATURE_GAP
                
                gene = feat['gene'] if feat.get('gene', None) else ''
                fh.write('\t'.join([feat['contig'], feat_type, str(feat['start']), str(feat['stop']), feat['strand'], gene, feat.get('product', ''), ', '.join(sorted(feat.get('db_xrefs', [])))]))
                fh.write('\n')
    return
