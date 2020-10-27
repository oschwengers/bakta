
import logging

from Bio import SeqIO

import bakta.constants as bc
import bakta.io.fasta as fasta
import bakta.so as so

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
                        'Name': feat['product'],
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
                    fh.write('\t'.join([feat['contig'], 'tRNAscan-SE', so.SO_TRNA.name, str(feat['start']), str(feat['stop']), '.', feat['strand'], '.', annotations]))
                    fh.write('\n')
                elif(feat['type'] is bc.FEATURE_TM_RNA):
                    annotations = {
                        'ID': feat['locus'],
                        'Name': feat['product'],
                        'locus_tag': feat['locus'],
                        'gene': feat['gene'],
                        'product': feat['product'],
                        'Dbxref': feat['db_xrefs']
                    }
                    annotations = encode_annotations(annotations)
                    fh.write('\t'.join([feat['contig'], 'Aragorn', so.SO_TMRNA.name, str(feat['start']), str(feat['stop']), '.', feat['strand'], '.', annotations]))
                    fh.write('\n')
                elif(feat['type'] is bc.FEATURE_R_RNA):
                    annotations = {
                        'ID': feat['locus'],
                        'Name': feat['product'],
                        'locus_tag': feat['locus'],
                        'gene': feat['gene'],
                        'product': feat['product'],
                        'Dbxref': feat['db_xrefs']
                    }
                    annotations = encode_annotations(annotations)
                    fh.write('\t'.join([feat['contig'], 'Infernal', so.SO_RRNA.name, str(feat['start']), str(feat['stop']), str(feat['evalue']), feat['strand'], '.', annotations]))
                    fh.write('\n')
                elif(feat['type'] is bc.FEATURE_NC_RNA):
                    annotations = {
                        'ID': feat['locus'],
                        'Name': feat['product'],
                        'locus_tag': feat['locus'],
                        'gene': feat['gene'],
                        'product': feat['product'],
                        'Dbxref': feat['db_xrefs']
                    }
                    annotations = encode_annotations(annotations)
                    fh.write('\t'.join([feat['contig'], 'Infernal', so.SO_NCRNA_GENE.name, str(feat['start']), str(feat['stop']), str(feat['evalue']), feat['strand'], '.', annotations]))
                    fh.write('\n')
                elif(feat['type'] is bc.FEATURE_NC_RNA_REGION):
                    annotations = {
                        'ID': feat['locus'],
                        'Name': feat['product'],
                        'locus_tag': feat['locus'],
                        'product': feat['product'],
                        'Dbxref': feat['db_xrefs']
                    }
                    annotations = encode_annotations(annotations)
                    fh.write('\t'.join([feat['contig'], 'Infernal', so.SO_REGULATORY_REGION.name, str(feat['start']), str(feat['stop']), str(feat['evalue']), feat['strand'], '.', annotations]))
                    fh.write('\n')
                elif(feat['type'] == bc.FEATURE_CRISPR):
                    annotations = {
                        'Name': feat['product'],
                        'product': feat['product']
                    }
                    annotations = encode_annotations(annotations)
                    fh.write('\t'.join([feat['contig'], 'PILER-CR', so.SO_CRISPR.name, str(feat['start']), str(feat['stop']), '.', feat['strand'], '0', annotations]))
                    fh.write('\n')
                elif(feat['type'] is bc.FEATURE_CDS):
                    annotations = {
                        'ID': feat['locus'],
                        'Name': bc.HYPOTHETICAL_PROTEIN if feat.get('hypothetical', False) else feat['product'],
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
                    fh.write('\t'.join([feat['contig'], 'Prodigal', so.SO_CDS.name, str(feat['start']), str(feat['stop']), '.', feat['strand'], '0', annotations]))
                    fh.write('\n')
                elif(feat['type'] is bc.FEATURE_SORF):
                    annotations = {
                        'ID': feat['locus'],
                        'Name': bc.HYPOTHETICAL_PROTEIN if feat.get('hypothetical', False) else feat['product'],
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
                    fh.write('\t'.join([feat['contig'], 'Bakta', so.SO_CDS.name, str(feat['start']), str(feat['stop']), '.', feat['strand'], '0', annotations]))
                    fh.write('\n')
                elif(feat['type'] is bc.FEATURE_GAP):
                    annotations = {
                        'Name': 'assembly gap [length=%s]' % feat['length'],
                        'product': 'assembly gap [length=%s]' % feat['length']
                    }
                    annotations = encode_annotations(annotations)
                    fh.write('\t'.join([feat['contig'], 'Bakta', so.SO_GAP.name, str(feat['start']), str(feat['stop']), '.', feat['strand'], '0', annotations]))
                    fh.write('\n')
                elif(feat['type'] == bc.FEATURE_ORIC):
                    annotations = {
                        'Name': 'oriC',
                        'product': 'oriC'
                    }
                    annotations = encode_annotations(annotations)
                    fh.write('\t'.join([feat['contig'], 'Blast+', so.SO_ORIC.name, str(feat['start']), str(feat['stop']), '.', feat['strand'], '0', annotations]))
                    fh.write('\n')
                elif(feat['type'] == bc.FEATURE_ORIV):
                    annotations = {
                        'Name': 'oriV',
                        'product': 'oriV'
                    }
                    annotations = encode_annotations(annotations)
                    fh.write('\t'.join([feat['contig'], 'Blast+', so.SO_ORIV.name, str(feat['start']), str(feat['stop']), '.', feat['strand'], '0', annotations]))
                    fh.write('\n')
                elif(feat['type'] == bc.FEATURE_ORIT):
                    annotations = {
                        'Name': 'oriT',
                        'product': 'oriT'
                    }
                    annotations = encode_annotations(annotations)
                    fh.write('\t'.join([feat['contig'], 'Blast+', so.SO_ORIT.name, str(feat['start']), str(feat['stop']), '.', feat['strand'], '0', annotations]))
                    fh.write('\n')

        fh.write('##FASTA\n')
        for contig in contigs:  # write sequences
            fh.write('>%s\n' % contig['id'])
            fh.write(fasta.wrap_sequence(contig['sequence']))
    
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
