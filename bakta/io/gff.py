
import logging

from Bio import SeqIO

import bakta.constants as bc
import bakta.io.fasta as fasta

log = logging.getLogger('io:gff')

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

def write_gff3(contigs, features_by_contig, gff3_path):
    """Export features in GFF3 format."""

    with gff3_path.open('w') as fh:
        fh.write('##gff-version 3\n')
        for contig in contigs:  # write features
            fh.write('##sequence-region %s %i %i\n' % (contig['id'], 1, contig['length']))
            for feat in features_by_contig[contig['id']]:
                if(feat['type'] is bc.FEATURE_T_RNA):
                    annotations = {
                        'ID': feat['locus'],
                        'NAME': feat['product'],
                        'locus_tag': feat['locus'],
                        'product': feat['product'],
                        'Dbxref': feat['db_xrefs']
                    }
                    # add gene annotation if available
                    if(feat.get('gene', None)):
                        annotations['gene'] = feat['gene']
                    if(feat.get('pseudo', False)):
                        annotations['pseudo'] = True
                    annotations = encode_annotations(annotations)
                    fh.write('\t'.join([feat['contig'], 'tRNAscan-SE', bc.FEATURE_T_RNA, str(feat['start']), str(feat['stop']), '.', feat['strand'], '.', annotations]))
                    fh.write('\n')
                elif(feat['type'] is bc.FEATURE_TM_RNA):
                    annotations = {
                        'ID': feat['locus'],
                        'NAME': feat['product'],
                        'locus_tag': feat['locus'],
                        'gene': feat['gene'],
                        'product': feat['product'],
                        'Dbxref': feat['db_xrefs']
                    }
                    annotations = encode_annotations(annotations)
                    fh.write('\t'.join([feat['contig'], 'Aragorn', bc.FEATURE_TM_RNA, str(feat['start']), str(feat['stop']), '.', feat['strand'], '.', annotations]))
                    fh.write('\n')
                elif(feat['type'] is bc.FEATURE_R_RNA):
                    annotations = {
                        'ID': feat['locus'],
                        'NAME': feat['product'],
                        'locus_tag': feat['locus'],
                        'gene': feat['gene'],
                        'product': feat['product'],
                        'Dbxref': feat['db_xrefs']
                    }
                    annotations = encode_annotations(annotations)
                    fh.write('\t'.join([feat['contig'], 'Infernal', bc.FEATURE_R_RNA, str(feat['start']), str(feat['stop']), str(feat['evalue']), feat['strand'], '.', annotations]))
                    fh.write('\n')
                elif(feat['type'] is bc.FEATURE_NC_RNA):
                    annotations = {
                        'ID': feat['locus'],
                        'NAME': feat['product'],
                        'locus_tag': feat['locus'],
                        'gene': feat['gene'],
                        'product': feat['product'],
                        'Dbxref': feat['db_xrefs']
                    }
                    annotations = encode_annotations(annotations)
                    fh.write('\t'.join([feat['contig'], 'Infernal', bc.FEATURE_NC_RNA, str(feat['start']), str(feat['stop']), str(feat['evalue']), feat['strand'], '.', annotations]))
                    fh.write('\n')
                elif(feat['type'] is bc.FEATURE_NC_RNA_REGION):
                    annotations = {
                        'ID': feat['locus'],
                        'NAME': feat['product'],
                        'locus_tag': feat['locus'],
                        'product': feat['product'],
                        'Dbxref': feat['db_xrefs']
                    }
                    annotations = encode_annotations(annotations)
                    fh.write('\t'.join([feat['contig'], 'Infernal', 'MISC_RNA', str(feat['start']), str(feat['stop']), str(feat['evalue']), feat['strand'], '.', annotations]))
                    fh.write('\n')
                elif(feat['type'] == bc.FEATURE_CRISPR):
                    annotations = {
                        'NAME': feat['product'],
                        'product': feat['product']
                    }
                    annotations = encode_annotations(annotations)
                    fh.write('\t'.join([feat['contig'], 'PILER-CR', 'CRISPR', str(feat['start']), str(feat['stop']), '.', feat['strand'], '0', annotations]))
                    fh.write('\n')
                elif(feat['type'] is bc.FEATURE_CDS):
                    annotations = {
                        'ID': feat['locus'],
                        'NAME': bc.HYPOTHETICAL_PROTEIN if feat.get('hypothetical', False) else feat['product'],
                        'locus_tag': feat['locus'],
                        'product': bc.HYPOTHETICAL_PROTEIN if feat.get('hypothetical', False) else feat['product']
                    }
                    if('db_xrefs' in feat):
                        annotations['Dbxref'] = feat['db_xrefs']
                    # add gene annotation if available
                    if(feat.get('gene', None)):
                        annotations['gene'] = feat['gene']
                    # add DbXrefs
                    annotations = encode_annotations(annotations)
                    fh.write('\t'.join([feat['contig'], 'Prodigal', 'CDS', str(feat['start']), str(feat['stop']), '.', feat['strand'], '0', annotations]))
                    fh.write('\n')
                elif(feat['type'] is bc.FEATURE_SORF):
                    annotations = {
                        'ID': feat['locus'],
                        'NAME': bc.HYPOTHETICAL_PROTEIN if feat.get('hypothetical', False) else feat['product'],
                        'locus_tag': feat['locus'],
                        'product': bc.HYPOTHETICAL_PROTEIN if feat.get('hypothetical', False) else feat['product']
                    }
                    if('db_xrefs' in feat):
                        annotations['Dbxref'] = feat['db_xrefs']
                    # add gene annotation if available
                    if(feat.get('gene', None)):
                        annotations['gene'] = feat['gene']
                    # add DbXrefs
                    annotations = encode_annotations(annotations)
                    fh.write('\t'.join([feat['contig'], 'Bakta', 'CDS', str(feat['start']), str(feat['stop']), '.', feat['strand'], '0', annotations]))
                    fh.write('\n')
                elif(feat['type'] is bc.FEATURE_GAP):
                    annotations = {
                        'NAME': 'assembly gap [length=%s]' % feat['length'],
                        'product': 'assembly gap [length=%s]' % feat['length']
                    }
                    annotations = encode_annotations(annotations)
                    feat_type = bc.INSDC_FEATURE_ASSEMBLY_GAP if feat['length'] >= 100 else bc.INSDC_FEATURE_GAP
                    fh.write('\t'.join([feat['contig'], 'Bakta', feat_type, str(feat['start']), str(feat['stop']), '.', feat['strand'], '0', annotations]))
                    fh.write('\n')
                elif(feat['type'] == bc.FEATURE_ORIC):
                    annotations = {
                        'NAME': 'oriC',
                        'product': 'oriC'
                    }
                    annotations = encode_annotations(annotations)
                    fh.write('\t'.join([feat['contig'], 'Blast+', bc.FEATURE_ORIC, str(feat['start']), str(feat['stop']), '.', feat['strand'], '0', annotations]))
                    fh.write('\n')
                elif(feat['type'] == bc.FEATURE_ORIT):
                    annotations = {
                        'NAME': 'oriT',
                        'product': 'oriT'
                    }
                    annotations = encode_annotations(annotations)
                    fh.write('\t'.join([feat['contig'], 'Blast+', bc.FEATURE_ORIT, str(feat['start']), str(feat['stop']), '.', feat['strand'], '0', annotations]))
                    fh.write('\n')

        fh.write('##FASTA\n')
        for contig in contigs:  # write sequences
            fh.write(fasta.format_fasta(contig, line_wrapping=True))
    
    return


def encode_annotations(annotations):
    annotation_strings = []
    for key, val in annotations.items():
        if(type(val) is list):
            if(len(val) >= 1):
                annotation = "%s=%s" % (key, ','.join(val) if type(val) is list else val)
                annotation_strings.append(annotation)
        else:
            annotation_strings.append("%s=%s" % (key, val))
    return ';'.join(annotation_strings)
