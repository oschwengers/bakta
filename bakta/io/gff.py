
import logging

import bakta
import bakta.config as cfg
import bakta.constants as bc
import bakta.io.fasta as fasta
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
                    # add gene annotation if available
                    if(feat.get('gene', None)):
                        annotations['gene'] = feat['gene']
                    if(feat.get('pseudo', False)):
                        annotations['pseudo'] = True
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
                    annotations = encode_annotations(annotations)
                    fh.write(f"{feat['contig']}\tInfernal\t{so.SO_NCRNA_GENE.name}\t{start}\t{stop}\t{feat['evalue']}\t{feat['strand']}\t.\t{annotations}\n")
                elif(feat['type'] is bc.FEATURE_NC_RNA_REGION):
                    annotations = {
                        'Name': feat['product'],
                        'product': feat['product'],
                        'Dbxref': feat['db_xrefs']
                    }
                    annotations = encode_annotations(annotations)
                    fh.write(f"{feat['contig']}\tInfernal\t{so.SO_REGULATORY_REGION.name}\t{start}\t{stop}\t{feat['evalue']}\t{feat['strand']}\t.\t{annotations}\n")
                elif(feat['type'] == bc.FEATURE_CRISPR):
                    annotations = {
                        'ID': feat['locus'],
                        'Name': feat['product'],
                        'locus_tag': feat['locus'],
                        'product': feat['product']
                    }
                    annotations = encode_annotations(annotations)
                    fh.write(f"{feat['contig']}\tPILER-CR\t{so.SO_CRISPR.name}\t{start}\t{stop}\t.\t{feat['strand']}\t.\t{annotations}\n")
                elif(feat['type'] is bc.FEATURE_CDS):
                    annotations = {
                        'ID': feat['locus'],
                        'Name': feat['product'],
                        'locus_tag': feat['locus'],
                        'product': feat['product']
                    }
                    if('db_xrefs' in feat):
                        annotations['Dbxref'] = feat['db_xrefs']
                    # add gene annotation if available
                    if(feat.get('gene', None)):
                        annotations['gene'] = feat['gene']
                    # add DbXrefs
                    annotations = encode_annotations(annotations)
                    fh.write(f"{feat['contig']}\tProdigal\t{so.SO_CDS.name}\t{start}\t{stop}\t.\t{feat['strand']}\t0\t{annotations}\n")
                elif(feat['type'] is bc.FEATURE_SORF):
                    annotations = {
                        'ID': feat['locus'],
                        'Name': feat['product'],
                        'locus_tag': feat['locus'],
                        'product': feat['product']
                    }
                    if('db_xrefs' in feat):
                        annotations['Dbxref'] = feat['db_xrefs']
                    # add gene annotation if available
                    if(feat.get('gene', None)):
                        annotations['gene'] = feat['gene']
                    # add DbXrefs
                    annotations = encode_annotations(annotations)
                    fh.write(f"{feat['contig']}\tBakta\t{so.SO_CDS.name}\t{start}\t{stop}\t.\t{feat['strand']}\t0\t{annotations}\n")
                elif(feat['type'] is bc.FEATURE_GAP):
                    annotations = {
                        'Name': f"gap ({feat['length']} bp)",
                        'product': f"gap ({feat['length']} bp)"
                    }
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
                    fh.write(f"{feat['contig']}\tBlast+\t{so.SO_ORIC.name}\t{start}\t{stop}\t.\t{feat['strand']}\t.\t{annotations}\n")
                elif(feat['type'] == bc.FEATURE_ORIV):
                    annotations = {
                        'ID': feat['locus'],
                        'Name': 'oriV',
                        'locus_tag': feat['locus'],
                        'product': 'oriV'
                    }
                    annotations = encode_annotations(annotations)
                    fh.write(f"{feat['contig']}\tBlast+\t{so.SO_ORIV.name}\t{start}\t{stop}\t.\t{feat['strand']}\t.\t{annotations}\n")
                elif(feat['type'] == bc.FEATURE_ORIT):
                    annotations = {
                        'ID': feat['locus'],
                        'Name': 'oriT',
                        'locus_tag': feat['locus'],
                        'product': 'oriT'
                    }
                    annotations = encode_annotations(annotations)
                    fh.write(f"{feat['contig']}\tBlast+\t{so.SO_ORIT.name}\t{start}\t{stop}\t.\t{feat['strand']}\t.\t{annotations}\n")

        fh.write('##FASTA\n')
        for contig in genome['contigs']:  # write sequences
            fh.write(f">{contig['id']}\n")
            fh.write(fasta.wrap_sequence(contig['sequence']))
    
    return


def encode_annotations(annotations):
    annotation_strings = []
    for key, val in annotations.items():
        if(type(val) is list):
            if(len(val) >= 1):
                annotation = f"{key}={','.join(val)}" if type(val) is list else val
                annotation_strings.append(annotation)
        else:
            annotation_strings.append(f'{key}={val}')
    return ';'.join(annotation_strings)
