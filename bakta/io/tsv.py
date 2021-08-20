import logging

import bakta
import bakta.config as cfg
import bakta.constants as bc


log = logging.getLogger('TSV')


def write_tsv(contigs, features_by_contig, tsv_path):
    """Export features in TSV format."""
    log.info('write tsv: path=%s', tsv_path)

    with tsv_path.open('wt') as fh:
        fh.write(f'#Annotated with Bakta (v{bakta.__version__}): https://github.com/oschwengers/bakta\n')
        fh.write(f"#Database (v{cfg.db_info['major']}.{cfg.db_info['minor']}): https://doi.org/10.5281/zenodo.4247252\n")
        fh.write('#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tGene\tProduct\tDbXrefs\n')
        for contig in contigs:
            for feat in features_by_contig[contig['id']]:
                feat_type = feat['type']
                if(feat['type'] == bc.FEATURE_GAP):
                    feat_type = bc.INSDC_FEATURE_ASSEMBLY_GAP if feat['length'] >= 100 else bc.INSDC_FEATURE_GAP

                gene = feat['gene'] if feat.get('gene', None) else ''
                fh.write('\t'.join([feat['contig'], feat_type, str(feat['start']), str(feat['stop']), feat['strand'], feat.get('locus', ''), gene, feat.get('product', ''), ', '.join(sorted(feat.get('db_xrefs', [])))]))
                fh.write('\n')
    return


def write_hypothetical_tsv(hypotheticals, tsv_path):
    """Export hypothetical information in TSV format."""
    log.info('write hypothetical tsv: path=%s', tsv_path)

    with tsv_path.open('wt') as fh:
        fh.write(f'#Annotated with Bakta v{bakta.__version__}, https://github.com/oschwengers/bakta\n')
        fh.write(f"#Database v{cfg.db_info['major']}.{cfg.db_info['minor']}, https://doi.org/10.5281/zenodo.4247252\n")
        fh.write('#Sequence Id\tStart\tStop\tStrand\tLocus Tag\tMol Weight [kDa]\tIso El. Point\tPfam hits\tDbxrefs\n')
        for hypo in hypotheticals:
            pfams = [f"{pfam['id']}|{pfam['name']}" for pfam in hypo.get('pfams', [])]
            seq_stats = hypo['seq_stats']
            mol_weight = f"{(seq_stats['molecular_weight']/100):.1f}" if seq_stats['molecular_weight'] else 'NA'
            iso_point = f"{seq_stats['isoelectric_point']:.1f}" if seq_stats['isoelectric_point'] else 'NA'
            fh.write(f"{hypo['contig']}\t{hypo['start']}\t{hypo['stop']}\t{hypo['strand']}\t{hypo.get('locus', '')}\t{mol_weight}\t{iso_point}\t{', '.join(sorted(pfams))}\t{', '.join(sorted(hypo.get('db_xrefs', [])))}\n")
    return
