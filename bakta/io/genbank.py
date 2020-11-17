import logging
import json

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation

from datetime import date

import bakta.constants as bc
import bakta.config as cfg
import bakta.psc as psc

log = logging.getLogger('io:genbank')

def write_genbank(genome, features, genbank_path):
    contig_list = []
    time = date.today()
    for contig in genome['contigs']:
        contig_annotations = {
            'molecule_type': 'DNA',
            'source': genome['taxon'],
            'organism': genome['taxon'],
            'date': time.strftime('%d-%b-%Y').upper(),
            'topology': contig['topology'],
            'data_file_division': 'HGT' if contig['type'] == bc.REPLICON_CONTIG else 'BCT'
            # TODO: taxonomy
        }
        if(contig.get('name', '')):
            if(contig['type'] == bc.REPLICON_PLASMID):
                description = f"{genome['taxon']} plasmid {contig['name']}"
            else:
                description = f"{genome['taxon']} {contig['name']}"
        else:
            description = f"{genome['taxon']} {contig['type']}"
        if(contig['complete']):
                description = f"{description}, complete sequece"

        contig_rec = SeqIO.SeqRecord(id='.', name=contig['id'], description=description, annotations=contig_annotations, seq=Seq(contig['sequence']))

        seq_feature_list = []
        source_qualifiers = {
            'organism': genome['taxon'],
            'mol_type': 'DNA',
            'strain': cfg.strain
        }
        source = SeqFeature(FeatureLocation(0, contig['length'], strand=+1), type='source', qualifiers=source_qualifiers)
        seq_feature_list.append(source)

        for feature in features:
            if(feature['contig'] == contig['id']):
                insdc_feature_type = None
                qualifiers = {}
                if('db_xrefs' in feature):
                    qualifiers['db_xref'] = feature['db_xrefs']
                if('product' in feature):
                    qualifiers['product'] = feature['product']
                if('locus' in feature):
                    qualifiers['locus_tag'] = feature['locus']
                
                if(feature['type'] == bc.FEATURE_GAP):
                    insdc_feature_type = bc.INSDC_FEATURE_GAP
                    qualifiers['estimated_length'] = feature['length']
                elif(feature['type'] == bc.FEATURE_ORIC or feature['type'] == bc.FEATURE_ORIV):
                    #TODO: Add fuzzy positions for oriC/oriV
                    insdc_feature_type = bc.INSDC_FEATURE_ORIGIN_REPLICTION
                elif(feature['type'] == feature['type'] == bc.FEATURE_ORIT):
                    #TODO: Add fuzzy positions for oriT
                    insdc_feature_type = bc.INSDC_FEATURE_ORIGIN_TRANSFER
                elif(feature['type'] == bc.FEATURE_CDS) or (feature['type'] == bc.FEATURE_SORF):
                    qualifiers['translation'] = feature['sequence']
                    qualifiers['codon_start'] = 1
                    insdc_feature_type = bc.INSDC_FEATURE_CDS
                    inference = []
                    inference.append('ab initio prediction:Prodigal:2.6' if feature['type'] == bc.FEATURE_CDS else 'ab initio prediction:Bakta')
                    if('hypothetical' not in feature):
                        if('ups' in feature):
                            if('ncbi_nrp_id' in feature['ups']):
                                qualifiers['protein_id'] = feature['ups']['ncbi_nrp_id']
                        if('ips' in feature):
                            if('uniref100_id' in feature['ips']):
                                ips_subject_id = feature['ips']['uniref100_id']
                                inference.append(f'similar to AA sequence:UniProtKB:{ips_subject_id}')
                        if('psc' in feature):
                            if('uniref90_id' in feature['psc']):
                                psc_subject_id = feature['psc']['uniref90_id']
                                inference.append(f'similar to AA sequence:UniProtKB:{psc_subject_id}')
                    qualifiers['inference'] = inference
                elif(feature['type'] == bc.FEATURE_T_RNA):
                    # TODO: Position anticodon
                    if('notes' in feature):
                        if('anti_codon' in feature):
                            qualifiers['note'] = feature['notes']
                            i = feature['product'].find('-')
                            t_rna_type = feature['product'][i+1:].lower()
                            qualifiers['anticodon'] = f"(aa:{t_rna_type},seq:{feature['anti_codon'].lower()})"
                    qualifiers['inference'] = 'profile:tRNAscan:2.0'
                    insdc_feature_type = bc.INSDC_FEATURE_T_RNA
                elif(feature['type'] == bc.FEATURE_TM_RNA):
                    qualifiers['inference'] = 'profile:aragorn:1.2'
                    insdc_feature_type = bc.INSDC_FEATURE_TM_RNA
                elif(feature['type'] == bc.FEATURE_R_RNA):
                    for reference in feature['db_xrefs']:
                        if(reference.split(':')[0] == 'RFAM'):
                            r_subject_id = reference.split(':')[1]
                    qualifiers['inference'] = f'profile:Rfam:{r_subject_id}'
                    insdc_feature_type = bc.INSDC_FEATURE_R_RNA
                elif(feature['type'] == bc.FEATURE_NC_RNA):
                    # TODO: ncRNA_class
                    for reference in feature['db_xrefs']:
                        if(reference.split(':')[0] == 'RFAM'):
                            nc_subject_id = reference.split(':')[1]
                    qualifiers['ncRNA_class'] = 'other'
                    qualifiers['inference'] = f'profile:Rfam:{nc_subject_id}'
                    insdc_feature_type = bc.INSDC_FEATURE_NC_RNA
                elif(feature['type']==bc.FEATURE_NC_RNA_REGION):
                    for reference in feature['db_xrefs']:
                        if(reference.split(':')[0] == 'RFAM'):
                            nc_subject_id = reference.split(':')[1]
                    qualifiers['ncRNA_class'] = 'other'
                    qualifiers['inference'] = f'profile:Rfam:{nc_subject_id}'
                    insdc_feature_type = bc.INSDC_FEATURE_REGULATORY
                elif(feature['type']==bc.FEATURE_CRISPR):
                    qualifiers['repeats'] = feature['repeats']
                    qualifiers['repeat_consensus'] = feature['repeat_consensus']
                    qualifiers['repeat_length'] = feature['repeat_length']
                    qualifiers['spacer_length'] = feature['spacer_length']
                    feature['type'] = 'misc_feature'
                    insdc_feature_type = bc.INSDC_FEATURE_MISC_FEATURE

                strand_dir = -1 if feature['strand'] == bc.STRAND_REVERSE else +1

                if(feature.get('gene', '') and feature['type'] != bc.FEATURE_NC_RNA_REGION):
                    qualifiers['gene'] = feature['gene']
                    gene_qualifier = {
                        'gene': feature['gene'],
                        'locus_tag': feature['locus']
                    }
                    gen_seqfeat = SeqFeature(FeatureLocation(feature['start'] - 1, feature['stop'], strand=strand_dir), type='gene', qualifiers=gene_qualifier)
                    seq_feature_list.append(gen_seqfeat)
                feat_seqfeat = SeqFeature(FeatureLocation(feature['start'] - 1, feature['stop'], strand=strand_dir), type=insdc_feature_type, qualifiers=qualifiers)
                seq_feature_list.append(feat_seqfeat)
        contig_rec.features = seq_feature_list
        contig_list.append(contig_rec)

    with genbank_path.open('w', encoding='utf-8') as fh:
        SeqIO.write(contig_list, fh, format='genbank')

    return
