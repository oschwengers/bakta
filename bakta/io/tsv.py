import logging

from pathlib import Path
from types import LambdaType
from typing import Dict, Sequence

import bakta
import bakta.config as cfg
import bakta.constants as bc


log = logging.getLogger('TSV')


def write_features(contigs: Sequence[dict], features_by_contig: Dict[str, dict], tsv_path: Path):
    """Export features in TSV format."""
    log.info('write feature tsv: path=%s', tsv_path)

    with tsv_path.open('wt') as fh:
        fh.write('# Annotated with Bakta\n')
        fh.write(f'# Software: v{bakta.__version__}\n')
        fh.write(f"# Database: v{cfg.db_info['major']}.{cfg.db_info['minor']}, {cfg.db_info['type']}\n")
        fh.write(f'# DOI: {bc.BAKTA_DOI}\n')
        fh.write(f'# URL: {bc.BAKTA_URL}\n')
        fh.write('#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tGene\tProduct\tDbXrefs\n')

        for contig in contigs:
            for feat in features_by_contig[contig['id']]:
                feat_type = feat['type']
                if(feat_type == bc.FEATURE_GAP):
                    feat_type = bc.INSDC_FEATURE_ASSEMBLY_GAP if feat['length'] >= 100 else bc.INSDC_FEATURE_GAP

                gene = feat['gene'] if feat.get('gene', None) else ''
                product = f"(pseudo) {feat.get('product', '')}" if feat.get('pseudo', False) else feat.get('product', '')
                fh.write('\t'.join(
                    [
                        feat['contig'],
                        feat_type,
                        str(feat['start']),
                        str(feat['stop']),
                        feat['strand'],
                        feat.get('locus', ''),
                        gene,
                        product,
                        ', '.join(sorted(feat.get('db_xrefs', [])))
                    ])
                )
                fh.write('\n')
                if(feat_type == bc.FEATURE_CRISPR):
                    i = 0
                    while i < len(feat['spacers']):
                        repeat = feat['repeats'][i]
                        fh.write('\t'.join([feat['contig'], bc.FEATURE_CRISPR_REPEAT, str(repeat['start']), str(repeat['stop']), repeat['strand'], '', '', f"CRISPR repeat", '']))
                        fh.write('\n')
                        spacer = feat['spacers'][i]
                        fh.write('\t'.join([feat['contig'], bc.FEATURE_CRISPR_SPACER, str(spacer['start']), str(spacer['stop']), spacer['strand'], '', '', f"CRISPR spacer, sequence {spacer['sequence']}", '']))
                        fh.write('\n')
                        i += 1
                    if(len(feat['repeats']) - 1 == i):
                        repeat = feat['repeats'][i]
                        fh.write('\t'.join([feat['contig'], bc.FEATURE_CRISPR_REPEAT, str(repeat['start']), str(repeat['stop']), repeat['strand'], '', '', f"CRISPR repeat", '']))
                        fh.write('\n')
    return


def write_feature_inferences(contigs: Sequence[dict], features_by_contig: Dict[str, dict], tsv_path: Path):
    """Export feature inference statistics in TSV format."""
    log.info('write tsv: path=%s', tsv_path)

    with tsv_path.open('wt') as fh:
        fh.write('# Annotated with Bakta\n')
        fh.write(f'# Software: v{bakta.__version__}\n')
        fh.write(f"# Database: v{cfg.db_info['major']}.{cfg.db_info['minor']}, {cfg.db_info['type']}\n")
        fh.write(f'# DOI: {bc.BAKTA_DOI}\n')
        fh.write(f'# URL: {bc.BAKTA_URL}\n')
        fh.write('#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tscore\tevalue\tquery-cov\tsubject-cov\tid')

        for contig in contigs:
            for feat in features_by_contig[contig['id']]:
                if(feat['type'] in [
                            bc.FEATURE_T_RNA,
                            bc.FEATURE_R_RNA,
                            bc.FEATURE_NC_RNA,
                            bc.FEATURE_NC_RNA_REGION,
                            bc.FEATURE_CDS,
                            bc.FEATURE_SORF
                        ]
                    ):
                    fh.write('\t'.join(
                        [
                            feat['contig'],
                            feat['type'],
                            str(feat['start']),
                            str(feat['stop']),
                            feat['strand'],
                            feat['locus'],
                            f"{feat['score']:0.1f}",
                            f"{feat['evalue']:1.1e}" if 'evalue' in feat else '-',
                            f"{feat['query_cov']:0.3f}" if 'query_cov' in feat else '-',
                            f"{feat['subject_cov']:0.3f}" if 'subject_cov' in feat else '-',
                            f"{feat['identity']:0.3f}" if 'identity' in feat else '-'
                        ])
                    )
                    fh.write('\n')
    return


def write_protein_features(features: Sequence[dict], header_columns: Sequence[str], mapping: LambdaType, tsv_path: Path):
    """Export protein features in TSV format."""
    log.info('write protein feature tsv: path=%s', tsv_path)

    with tsv_path.open('wt') as fh:
        fh.write(f'#Annotated with Bakta (v{bakta.__version__}): https://github.com/oschwengers/bakta\n')
        fh.write(f"#Database (v{cfg.db_info['major']}.{cfg.db_info['minor']}): https://doi.org/10.5281/zenodo.4247252\n")
        fh.write('\t'.join(header_columns))
        fh.write('\n')
        for feat in features:
            columns = mapping(feat)
            fh.write('\t'.join(columns))
            fh.write('\n')
    return


def write_hypotheticals(hypotheticals: Sequence[dict], tsv_path: Path):
    """Export hypothetical information in TSV format."""
    log.info('write hypothetical tsv: path=%s', tsv_path)

    with tsv_path.open('wt') as fh:
        fh.write(f'#Annotated with Bakta v{bakta.__version__}, https://github.com/oschwengers/bakta\n')
        fh.write(f"#Database v{cfg.db_info['major']}.{cfg.db_info['minor']}, https://doi.org/10.5281/zenodo.4247252\n")
        fh.write('#Sequence Id\tStart\tStop\tStrand\tLocus Tag\tMol Weight [kDa]\tIso El. Point\tPfam hits\tDbxrefs\n')
        for hypo in hypotheticals:
            pfams = [f"{pfam['id']}|{pfam['name']}" for pfam in hypo.get('pfams', [])]
            seq_stats = hypo['seq_stats']
            mol_weight = f"{(seq_stats['molecular_weight']/1000):.1f}" if seq_stats['molecular_weight'] else 'NA'
            iso_point = f"{seq_stats['isoelectric_point']:.1f}" if seq_stats['isoelectric_point'] else 'NA'
            fh.write(f"{hypo['contig']}\t{hypo['start']}\t{hypo['stop']}\t{hypo['strand']}\t{hypo.get('locus', '')}\t{mol_weight}\t{iso_point}\t{', '.join(sorted(pfams))}\t{', '.join(sorted(hypo.get('db_xrefs', [])))}\n")
    return
