import logging

import bakta
import bakta.config as cfg
import bakta.constants as bc
import bakta.io.fasta as fasta
import bakta.io.insdc as insdc
import bakta.so as so


log = logging.getLogger('GFF')


def write_gff3(genome, features_by_contig, gff3_path):
    """Export features in GFF3 format."""
    log.info('write GFF3: path=%s', gff3_path)

    with gff3_path.open('wt') as fh:
        fh.write('##gff-version 3\n')  # GFF version
        fh.write('##feature-ontology https://github.com/The-Sequence-Ontology/SO-Ontologies/blob/v3.1/so.obo\n')  # SO feature version

        if(genome['taxon']):  # write organism info
            fh.write(f"# organism {genome['taxon']}\n")

        fh.write(f'# annotated with Bakta (v{bakta.__version__}): https://github.com/oschwengers/bakta\n')
        fh.write(f"# database (v{cfg.db_info['major']}.{cfg.db_info['minor']}): https://doi.org/10.5281/zenodo.4247252\n")

        feature_id_counter = 1
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
                            annotations['pseudo'] = 'true'
                            gene_annotations['pseudo'] = 'true'
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
                        for dbxref in feat['db_xrefs']:
                            if(dbxref.split(':')[0] == 'RFAM'):
                                rfam_id = dbxref.split(':')[1]
                                annotations['inference'] = f'profile:Rfam:{rfam_id}'
                        annotations['Dbxref'], annotations['Note'] = insdc.revise_dbxref_insdc(feat['db_xrefs'])  # remove INSDC invalid DbXrefs
                        gene_annotations = {
                            'ID': gene_id,
                            'locus_tag': feat['locus'],
                            'gene': feat['gene']
                        }
                        if(feat.get('pseudo', None)):
                            annotations['pseudo'] = 'true'
                            gene_annotations['pseudo'] = 'true'
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
                        for dbxref in feat['db_xrefs']:
                            if(dbxref.split(':')[0] == 'RFAM'):
                                rfam_id = dbxref.split(':')[1]
                                annotations['inference'] = f'profile:Rfam:{rfam_id}'
                        annotations['Dbxref'], annotations['Note'] = insdc.revise_dbxref_insdc(feat['db_xrefs'])  # remove INSDC invalid DbXrefs
                        annotations[bc.INSDC_FEATURE_NC_RNA_CLASS] = insdc.select_ncrna_class(feat)
                        gene_annotations = encode_annotations(gene_annotations)
                        fh.write(f"{feat['contig']}\tInfernal\tgene\t{start}\t{stop}\t.\t{feat['strand']}\t.\t{gene_annotations}\n")
                    annotations = encode_annotations(annotations)
                    fh.write(f"{feat['contig']}\tInfernal\t{so.SO_NCRNA_GENE.name}\t{start}\t{stop}\t{feat['evalue']}\t{feat['strand']}\t.\t{annotations}\n")
                elif(feat['type'] is bc.FEATURE_NC_RNA_REGION):
                    annotations = {
                        'ID': feature_id_counter,
                        'Name': feat['product'],
                        'product': feat['product'],
                        'Dbxref': feat['db_xrefs']
                    }
                    if(cfg.compliant):
                        annotations['Dbxref'], annotations['Note'] = insdc.revise_dbxref_insdc(feat['db_xrefs'])  # remove INSDC invalid DbXrefs
                        annotations[bc.INSDC_FEATURE_REGULATORY_CLASS] = insdc.select_regulatory_class(feat)
                    feature_id_counter += 1
                    annotations = encode_annotations(annotations)
                    fh.write(f"{feat['contig']}\tInfernal\t{so.SO_REGULATORY_REGION.name}\t{start}\t{stop}\t{feat['evalue']}\t{feat['strand']}\t.\t{annotations}\n")
                elif(feat['type'] == bc.FEATURE_CRISPR):
                    annotations = {
                        'ID': feat['locus'],
                        'Name': feat['product'],
                        'locus_tag': feat['locus'],
                        'product': feat['product']
                    }
                    feat_type = so.SO_CRISPR.name
                    if(cfg.compliant):
                        feat_type = bc.INSDC_FEATURE_REPEAT_REGION
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
                        for note in annotations['Note']:
                            if(bc.DB_XREF_EC in note):
                                annotations['ec_number'] = note.replace('EC:', '')
                        annotations['Note'] = [note for note in annotations['Note'] if bc.DB_XREF_EC not in note]
                        gene_annotations = encode_annotations(gene_annotations)
                        fh.write(f"{feat['contig']}\tProdigal\tgene\t{start}\t{stop}\t.\t{feat['strand']}\t.\t{gene_annotations}\n")
                    annotations = encode_annotations(annotations)
                    fh.write(f"{feat['contig']}\tProdigal\t{so.SO_CDS.name}\t{start}\t{stop}\t.\t{feat['strand']}\t0\t{annotations}\n")
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
                        for note in annotations['Note']:
                            if(bc.DB_XREF_EC in note):
                                annotations['ec_number'] = note.replace('EC:', '')
                        annotations['Note'] = [note for note in annotations['Note'] if bc.DB_XREF_EC not in note]
                        gene_annotations = encode_annotations(gene_annotations)
                        fh.write(f"{feat['contig']}\tBakta\tgene\t{start}\t{stop}\t.\t{feat['strand']}\t.\t{gene_annotations}\n")
                    annotations = encode_annotations(annotations)
                    fh.write(f"{feat['contig']}\tBakta\t{so.SO_CDS.name}\t{start}\t{stop}\t.\t{feat['strand']}\t0\t{annotations}\n")
                elif(feat['type'] is bc.FEATURE_GAP):
                    annotations = {
                        'ID': feature_id_counter,
                        'Name': f"gap ({feat['length']} bp)",
                        'product': f"gap ({feat['length']} bp)"
                    }
                    feature_id_counter += 1
                    annotations = encode_annotations(annotations)
                    fh.write(f"{feat['contig']}\tBakta\t{so.SO_GAP.name}\t{start}\t{stop}\t.\t{feat['strand']}\t.\t{annotations}\n")
                elif(feat['type'] == bc.FEATURE_ORIC):
                    annotations = {
                        'ID': feat['locus'],
                        'Name': 'oriC',
                        'locus_tag': feat['locus'],
                        'product': 'oriC'
                    }
                    annotations = encode_annotations(annotations)
                    feat_type = bc.INSDC_FEATURE_ORIGIN_REPLICATION if cfg.compliant else so.SO_ORIC.name
                    fh.write(f"{feat['contig']}\tBLAST+\t{feat_type}\t{start}\t{stop}\t.\t{feat['strand']}\t.\t{annotations}\n")
                elif(feat['type'] == bc.FEATURE_ORIV):
                    annotations = {
                        'ID': feat['locus'],
                        'Name': 'oriV',
                        'locus_tag': feat['locus'],
                        'product': 'oriV'
                    }
                    annotations = encode_annotations(annotations)
                    feat_type = bc.INSDC_FEATURE_ORIGIN_REPLICATION if cfg.compliant else so.SO_ORIC.name
                    fh.write(f"{feat['contig']}\tBLAST+\t{feat_type}\t{start}\t{stop}\t.\t{feat['strand']}\t.\t{annotations}\n")
                elif(feat['type'] == bc.FEATURE_ORIT):
                    annotations = {
                        'ID': feat['locus'],
                        'Name': 'oriT',
                        'locus_tag': feat['locus'],
                        'product': 'oriT'
                    }
                    annotations = encode_annotations(annotations)
                    feat_type = bc.INSDC_FEATURE_ORIGIN_TRANSFER if cfg.compliant else so.SO_ORIT.name
                    fh.write(f"{feat['contig']}\tBLAST+\t{feat_type}\t{start}\t{stop}\t.\t{feat['strand']}\t.\t{annotations}\n")

        if(not cfg.compliant):
            fh.write('##FASTA\n')
            for contig in genome['contigs']:  # write sequences
                fh.write(f">{contig['id']}\n")
                fh.write(fasta.wrap_sequence(contig['sequence']))

    return


def encode_attribute(product):
    """Replace special characters forbidden in column 9 of the GFF3 format: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md"""
    product = str(product)
    product = product.replace('%', '%25')
    product = product.replace(';', '%3B')
    product = product.replace('=', '%3D')
    product = product.replace('&', '%26')
    product = product.replace(',', '%2C')
    return product


def encode_annotations(annotations):
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
