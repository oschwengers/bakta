
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
                    fh.write('\t'.join([feat['contig'], 'tRNAscan-SE', 'tRNA', str(feat['start']), str(feat['stop']), 'NA', feat['strand'], '', annotations]))
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
                    fh.write('\t'.join([feat['contig'], 'Aragorn', 'tmRNA', str(feat['start']), str(feat['stop']), 'NA', feat['strand'], '', annotations]))
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
                    fh.write('\t'.join([feat['contig'], 'Infernal', 'rRNA', str(feat['start']), str(feat['stop']), 'NA', feat['strand'], '', annotations]))
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
                    fh.write('\t'.join([feat['contig'], 'Infernal', 'ncRNA', str(feat['start']), str(feat['stop']), 'NA', feat['strand'], '', annotations]))
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
                    fh.write('\t'.join([feat['contig'], 'Infernal', 'MISC_RNA', str(feat['start']), str(feat['stop']), 'NA', feat['strand'], '', annotations]))
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
                    fh.write('\t'.join([feat['contig'], 'Prodigal', 'CDS', str(feat['start']), str(feat['stop']), 'NA', feat['strand'], str(feat['frame']), annotations]))
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
                    fh.write('\t'.join([feat['contig'], 'Bakta', 'CDS', str(feat['start']), str(feat['stop']), 'NA', feat['strand'], str(feat['frame']), annotations]))
                    fh.write('\n')
        
        for contig in contigs:  # write sequences
            fh.write('##FASTA\n')
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
