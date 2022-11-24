import logging

from pathlib import Path
from typing import Dict, Sequence, Union

import bakta
import bakta.config as cfg
import bakta.constants as bc
import bakta.io.fasta as fasta
import bakta.io.insdc as insdc
import bakta.so as so


log = logging.getLogger('GFF')


def write_gff3(genome: dict, features_by_contig: Dict[str, dict], gff3_path: Path):
    """Export features in GFF3 format."""
    log.info('write GFF3: path=%s', gff3_path)

    with gff3_path.open('wt') as fh:
        fh.write('##gff-version 3\n')  # GFF version
        fh.write('##feature-ontology https://github.com/The-Sequence-Ontology/SO-Ontologies/blob/v3.1/so.obo\n')  # SO feature version

        if(genome['taxon']):  # write organism info
            fh.write(f"# organism {genome['taxon']}\n")

        fh.write('# Annotated with Bakta\n')
        fh.write(f'# Software: v{bakta.__version__}\n')
        fh.write(f"# Database: v{cfg.db_info['major']}.{cfg.db_info['minor']}\n")
        fh.write(f'# DOI: {bc.BAKTA_DOI}\n')
        fh.write(f'# URL: {bc.BAKTA_URL}\n')

        for contig in genome['contigs']:  # write features
            fh.write(f"##sequence-region {contig['id']} 1 {contig['length']}\n")  # sequence region

            # write landmark region
            annotations = {
                'ID': contig['id'],
                'Name': contig['id']
            }
            if(contig['topology'] == bc.TOPOLOGY_CIRCULAR):
                annotations['Is_circular'] = 'true'
            annotations = encode_annotations(annotations)
            fh.write(f"{contig['id']}\tBakta\tregion\t1\t{str(contig['length'])}\t.\t+\t.\t{annotations}\n")

            for feat in features_by_contig[contig['id']]:
                start = feat['start']
                stop = feat['stop']
                if('edge' in feat):
                    stop += contig['length']

                if(feat['type'] is bc.FEATURE_T_RNA):
                    annotations = {
                        'ID': feat['locus'],
                        'Name': feat['product'],
                        'locus_tag': feat['locus'],
                        'product': feat['product'],
                        'Dbxref': feat['db_xrefs']
                    }
                    if(feat.get('gene', None)):  # add gene annotation if available
                        annotations['gene'] = feat['gene']
                    if(feat.get('pseudo', False)):
                        annotations['pseudo'] = True
                    if(feat.get('anti_codon', False)):
                        annotations['anti_codon'] = feat['anti_codon']
                    if(feat.get('amino_acid', False)):
                        annotations['amino_acid'] = feat['amino_acid']
                    if(cfg.compliant):
                        gene_id = f"{feat['locus']}_gene"
                        annotations['Parent'] = gene_id
                        annotations['inference'] = 'profile:tRNAscan:2.0'
                        annotations['Dbxref'], annotations['Note'] = insdc.revise_dbxref_insdc(feat['db_xrefs'])  # remove INSDC invalid DbXrefs
                        gene_annotations = {
                            'ID': gene_id,
                            'locus_tag': feat['locus']
                        }
                        if(feat.get('gene', None)):
                            gene_annotations['gene'] = feat['gene']
                        if(feat.get('pseudo', None)):
                            annotations[bc.INSDC_FEATURE_PSEUDOGENE] = 'true'
                            gene_annotations[bc.INSDC_FEATURE_PSEUDOGENE] = 'true'
                        gene_annotations = encode_annotations(gene_annotations)
                        fh.write(f"{feat['contig']}\ttRNAscan-SE\tgene\t{start}\t{stop}\t.\t{feat['strand']}\t.\t{gene_annotations}\n")
                    annotations = encode_annotations(annotations)
                    fh.write(f"{feat['contig']}\ttRNAscan-SE\t{so.SO_TRNA.name}\t{start}\t{stop}\t.\t{feat['strand']}\t.\t{annotations}\n")
                elif(feat['type'] is bc.FEATURE_TM_RNA):
                    annotations = {
                        'ID': feat['locus'],
                        'Name': feat['product'],
                        'locus_tag': feat['locus'],
                        'gene': feat['gene'],
                        'product': feat['product'],
                        'Dbxref': feat['db_xrefs']
                    }
                    if(cfg.compliant):
                        gene_id = f"{feat['locus']}_gene"
                        annotations['Parent'] = gene_id
                        annotations['inference'] = 'profile:aragorn:1.2'
                        annotations['Dbxref'], annotations['Note'] = insdc.revise_dbxref_insdc(feat['db_xrefs'])  # remove INSDC invalid DbXrefs
                        gene_annotations = {
                            'ID': gene_id,
                            'locus_tag': feat['locus'],
                            'gene': feat['gene']
                        }
                        gene_annotations = encode_annotations(gene_annotations)
                        fh.write(f"{feat['contig']}\tAragorn\tgene\t{start}\t{stop}\t.\t{feat['strand']}\t.\t{gene_annotations}\n")
                    annotations = encode_annotations(annotations)
                    fh.write(f"{feat['contig']}\tAragorn\t{so.SO_TMRNA.name}\t{start}\t{stop}\t.\t{feat['strand']}\t.\t{annotations}\n")
                elif(feat['type'] is bc.FEATURE_R_RNA):
                    annotations = {
                        'ID': feat['locus'],
                        'Name': feat['product'],
                        'locus_tag': feat['locus'],
                        'gene': feat['gene'],
                        'product': feat['product'],
                        'Dbxref': feat['db_xrefs']
                    }
                    if(cfg.compliant):
                        gene_id = f"{feat['locus']}_gene"
                        annotations['Parent'] = gene_id
                        for rfam_id in [dbxref.split(':')[1] for dbxref in feat['db_xrefs'] if dbxref.split(':')[0] == bc.DB_XREF_RFAM]:
                            annotations['inference'] = f'profile:Rfam:{rfam_id}'
                        annotations['Dbxref'], annotations['Note'] = insdc.revise_dbxref_insdc(feat['db_xrefs'])  # remove INSDC invalid DbXrefs
                        gene_annotations = {
                            'ID': gene_id,
                            'locus_tag': feat['locus'],
                            'gene': feat['gene']
                        }
                        gene_annotations = encode_annotations(gene_annotations)
                        fh.write(f"{feat['contig']}\tInfernal\tgene\t{start}\t{stop}\t.\t{feat['strand']}\t.\t{gene_annotations}\n")
                    annotations = encode_annotations(annotations)
                    fh.write(f"{feat['contig']}\tInfernal\t{so.SO_RRNA.name}\t{start}\t{stop}\t{feat['evalue']}\t{feat['strand']}\t.\t{annotations}\n")
                elif(feat['type'] is bc.FEATURE_NC_RNA):
                    annotations = {
                        'ID': feat['locus'],
                        'Name': feat['product'],
                        'locus_tag': feat['locus'],
                        'gene': feat['gene'],
                        'product': feat['product'],
                        'Dbxref': feat['db_xrefs']
                    }
                    if(cfg.compliant):
                        gene_id = f"{feat['locus']}_gene"
                        gene_annotations = {
                            'ID': gene_id,
                            'locus_tag': feat['locus'],
                            'gene': feat['gene']
                        }
                        annotations['Parent'] = gene_id
                        for rfam_id in [dbxref.split(':')[1] for dbxref in feat['db_xrefs'] if dbxref.split(':')[0] == bc.DB_XREF_RFAM]:
                            annotations['inference'] = f'profile:Rfam:{rfam_id}'
                        annotations['Dbxref'], annotations['Note'] = insdc.revise_dbxref_insdc(feat['db_xrefs'])  # remove INSDC invalid DbXrefs
                        annotations[bc.INSDC_FEATURE_NC_RNA_CLASS] = insdc.select_ncrna_class(feat)
                        gene_annotations = encode_annotations(gene_annotations)
                        fh.write(f"{feat['contig']}\tInfernal\tgene\t{start}\t{stop}\t.\t{feat['strand']}\t.\t{gene_annotations}\n")
                    annotations = encode_annotations(annotations)
                    fh.write(f"{feat['contig']}\tInfernal\t{so.SO_NCRNA_GENE.name}\t{start}\t{stop}\t{feat['evalue']}\t{feat['strand']}\t.\t{annotations}\n")
                elif(feat['type'] is bc.FEATURE_NC_RNA_REGION):
                    annotations = {
                        'ID': feat['id'],
                        'Name': feat['product'],
                        'product': feat['product'],
                        'Dbxref': feat['db_xrefs']
                    }
                    if(cfg.compliant):
                        for rfam_id in [dbxref.split(':')[1] for dbxref in feat['db_xrefs'] if dbxref.split(':')[0] == bc.DB_XREF_RFAM]:
                            annotations['inference'] = f'profile:Rfam:{rfam_id}'
                        annotations['Dbxref'], annotations['Note'] = insdc.revise_dbxref_insdc(feat['db_xrefs'])  # remove INSDC invalid DbXrefs
                        annotations[bc.INSDC_FEATURE_REGULATORY_CLASS] = insdc.select_regulatory_class(feat)
                    annotations = encode_annotations(annotations)
                    fh.write(f"{feat['contig']}\tInfernal\t{so.SO_REGULATORY_REGION.name}\t{start}\t{stop}\t{feat['evalue']}\t{feat['strand']}\t.\t{annotations}\n")
                elif(feat['type'] == bc.FEATURE_CRISPR):
                    annotations = {
                        'ID': feat['id'],
                        'Name': feat['product'],
                        'product': feat['product']
                    }
                    feat_type = so.SO_CRISPR.name
                    if(cfg.compliant):
                        feat_type = bc.INSDC_FEATURE_REPEAT_REGION
                        annotations['inference'] = 'COORDINATES:alignment:pilercr:1.02'
                        annotations['Dbxref'], annotations['Note'] = insdc.revise_dbxref_insdc(feat['db_xrefs'])  # remove INSDC invalid DbXrefs
                        annotations[bc.INSDC_FEATURE_REPEAT_FAMILY] = 'CRISPR'
                        annotations[bc.INSDC_FEATURE_REPEAT_TYPE] = 'direct'
                        annotations[bc.INSDC_FEATURE_REPEAT_UNIT_SEQ] = feat['repeat_consensus']
                    annotations = encode_annotations(annotations)
                    fh.write(f"{feat['contig']}\tPILER-CR\t{feat_type}\t{start}\t{stop}\t.\t{feat['strand']}\t.\t{annotations}\n")
                elif(feat['type'] is bc.FEATURE_CDS):
                    annotations = {
                        'ID': feat['locus'],
                        'Name': feat['product'],
                        'locus_tag': feat['locus'],
                        'product': feat['product'],
                        'Dbxref': feat['db_xrefs']
                    }
                    if(feat.get('pseudo', False)):
                        annotations['pseudo'] = True
                    if(feat.get('gene', None)):  # add gene annotation if available
                        annotations['gene'] = feat['gene']
                    if(cfg.compliant):
                        gene_id = f"{feat['locus']}_gene"
                        gene_annotations = {
                            'ID': gene_id,
                            'locus_tag': feat['locus']
                        }
                        if(feat.get('gene', None)):
                            gene_annotations['gene'] = feat['gene']
                        annotations['Parent'] = gene_id
                        annotations['inference'] = 'ab initio prediction:Prodigal:2.6'
                        annotations['Dbxref'], annotations['Note'] = insdc.revise_dbxref_insdc(feat['db_xrefs'])  # remove INSDC invalid DbXrefs
                        annotations['Note'], ec_number = insdc.extract_ec_from_notes_insdc(annotations, 'Note')
                        if(feat.get('pseudo', False)):
                            annotations[bc.INSDC_FEATURE_PSEUDOGENE] = bc.INSDC_FEATURE_PSEUDOGENE_TYPE_UNPROCESSED if feat[bc.PSEUDOGENE]['paralog'] else bc.INSDC_FEATURE_PSEUDOGENE_TYPE_UNITARY
                            gene_annotations[bc.INSDC_FEATURE_PSEUDOGENE] = bc.INSDC_FEATURE_PSEUDOGENE_TYPE_UNPROCESSED if feat[bc.PSEUDOGENE]['paralog'] else bc.INSDC_FEATURE_PSEUDOGENE_TYPE_UNITARY
                        if(ec_number is not None):
                            annotations['ec_number'] = ec_number
                        gene_annotations = encode_annotations(gene_annotations)
                        fh.write(f"{feat['contig']}\tProdigal\tgene\t{start}\t{stop}\t.\t{feat['strand']}\t.\t{gene_annotations}\n")
                    if('exception' in feat):
                        ex = feat['exception']
                        pos = f"{ex['start']}..{ex['stop']}"
                        if(feat['strand'] == bc.STRAND_REVERSE):
                            pos = f"complement({pos})"
                        annotations['transl_except']=f"(pos:{pos},aa:{ex['aa']})"
                        notes = annotations.get('Note', [])
                        notes.append(f"codon on position {ex['codon_position']} is a {ex['type']} codon")
                        if('Notes' not in annotations):
                            annotations['Note'] = notes
                    annotations = encode_annotations(annotations)
                    fh.write(f"{feat['contig']}\tProdigal\t{so.SO_CDS.name}\t{start}\t{stop}\t.\t{feat['strand']}\t0\t{annotations}\n")
                    if(bc.FEATURE_SIGNAL_PEPTIDE in feat):
                        write_signal_peptide(fh, feat)
                elif(feat['type'] is bc.FEATURE_SORF):
                    annotations = {
                        'ID': feat['locus'],
                        'Name': feat['product'],
                        'locus_tag': feat['locus'],
                        'product': feat['product'],
                        'Dbxref': feat['db_xrefs']
                    }
                    if(feat.get('gene', None)):  # add gene annotation if available
                        annotations['gene'] = feat['gene']
                    if(cfg.compliant):
                        gene_id = f"{feat['locus']}_gene"
                        gene_annotations = {
                            'ID': gene_id,
                            'locus_tag': feat['locus']
                        }
                        if(feat.get('gene', None)):
                            gene_annotations['gene'] = feat['gene']
                        annotations['Parent'] = gene_id
                        annotations['inference'] = 'ab initio prediction:Bakta'
                        annotations['Dbxref'], annotations['Note'] = insdc.revise_dbxref_insdc(feat['db_xrefs'])  # remove INSDC invalid DbXrefs
                        annotations['Note'], ec_number = insdc.extract_ec_from_notes_insdc(annotations, 'Note')
                        if(ec_number is not None):
                            annotations['ec_number'] = ec_number
                        gene_annotations = encode_annotations(gene_annotations)
                        fh.write(f"{feat['contig']}\tBakta\tgene\t{start}\t{stop}\t.\t{feat['strand']}\t.\t{gene_annotations}\n")
                    annotations = encode_annotations(annotations)
                    fh.write(f"{feat['contig']}\tBakta\t{so.SO_CDS.name}\t{start}\t{stop}\t.\t{feat['strand']}\t0\t{annotations}\n")
                    if(bc.FEATURE_SIGNAL_PEPTIDE in feat):
                        write_signal_peptide(fh, feat)
                elif(feat['type'] is bc.FEATURE_GAP):
                    annotations = {
                        'ID': feat['id'],
                        'Name': f"gap ({feat['length']} bp)",
                        'product': f"gap ({feat['length']} bp)"
                    }
                    annotations = encode_annotations(annotations)
                    fh.write(f"{feat['contig']}\tBakta\t{so.SO_GAP.name}\t{start}\t{stop}\t.\t{feat['strand']}\t.\t{annotations}\n")
                elif(feat['type'] == bc.FEATURE_ORIC):
                    annotations = {
                        'ID': feat['id'],
                        'Name': feat['product']
                    }
                    if(cfg.compliant):
                        annotations['Note'] = feat['product']
                    else:
                        annotations['product'] = feat['product']
                        annotations['inference'] = 'similar to DNA sequence'
                    annotations = encode_annotations(annotations)
                    feat_type = bc.INSDC_FEATURE_ORIGIN_REPLICATION if cfg.compliant else so.SO_ORIC.name
                    fh.write(f"{feat['contig']}\tBLAST+\t{feat_type}\t{start}\t{stop}\t.\t{feat['strand']}\t.\t{annotations}\n")
                elif(feat['type'] == bc.FEATURE_ORIV):
                    annotations = {
                        'ID': feat['id'],
                        'Name': feat['product']
                    }
                    if(cfg.compliant):
                        annotations['Note'] = feat['product']
                    else:
                        annotations['product'] = feat['product']
                        annotations['inference'] = 'similar to DNA sequence'
                    annotations = encode_annotations(annotations)
                    feat_type = bc.INSDC_FEATURE_ORIGIN_REPLICATION if cfg.compliant else so.SO_ORIC.name
                    fh.write(f"{feat['contig']}\tBLAST+\t{feat_type}\t{start}\t{stop}\t.\t{feat['strand']}\t.\t{annotations}\n")
                elif(feat['type'] == bc.FEATURE_ORIT):
                    annotations = {
                        'ID': feat['id'],
                        'Name': feat['product']
                    }
                    if(cfg.compliant):
                        annotations['Note'] = feat['product']
                    else:
                        annotations['product'] = feat['product']
                        annotations['inference'] = 'similar to DNA sequence'
                    annotations = encode_annotations(annotations)
                    feat_type = bc.INSDC_FEATURE_ORIGIN_TRANSFER if cfg.compliant else so.SO_ORIT.name
                    fh.write(f"{feat['contig']}\tBLAST+\t{feat_type}\t{start}\t{stop}\t.\t{feat['strand']}\t.\t{annotations}\n")

        if(not cfg.compliant):
            fh.write('##FASTA\n')
            for contig in genome['contigs']:  # write sequences
                fh.write(f">{contig['id']}\n")
                fh.write(fasta.wrap_sequence(contig['sequence']))
    return


def encode_attribute(product: str) -> str:
    """Replace special characters forbidden in column 9 of the GFF3 format: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md"""
    product = str(product)
    product = product.replace('%', '%25')
    product = product.replace(';', '%3B')
    product = product.replace('=', '%3D')
    product = product.replace('&', '%26')
    product = product.replace(',', '%2C')
    return product


def encode_annotations(annotations: Dict[str, Union[str, Sequence[str]]]) -> str:
    annotation_strings = []
    for key, val in annotations.items():
        if(type(val) is list):
            if(len(val) >= 1):
                val = [encode_attribute(k) for k in val]
                annotation = f"{key}={','.join(val)}"
                annotation_strings.append(annotation)
        else:
            annotation_strings.append(f'{key}={encode_attribute(val)}')
    return ';'.join(annotation_strings)


def write_signal_peptide(fh, feat: dict):
    sig_peptide = feat[bc.FEATURE_SIGNAL_PEPTIDE]
    annotations = {
        'ID': f"{feat['locus']}_sigpep",
        'Name': 'signal peptide',
        'product': 'signal peptide',
        'score': sig_peptide['score'],
        'Parent': feat['locus']
    }
    annotations = encode_annotations(annotations)
    fh.write(f"{feat['contig']}\tDeepSig\t{so.SO_SIGNAL_PEPTIDE.name}\t{sig_peptide['start']}\t{sig_peptide['stop']}\t{sig_peptide['score']:.2f}\t{feat['strand']}\t.\t{annotations}\n")
