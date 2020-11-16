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
    with open(genbank_path, 'w', encoding='utf-8') as genbank_file:
        contig_list = []
        time = date.today()
        for contig in genome['contigs']:
            if contig['type'] == bc.REPLICON_CONTIG:
                division = "HGT"
            else:
                division = "BCT"
            contig_annotations = {"molecule_type": "DNA",
                                  "source": genome['taxon'],
                                  "organism": genome['taxon'],
                                  "date": time.strftime("%d-%b-%Y").upper(),
                                  "topology": contig['topology'],
                                  "data_file_division": division}
                                  #"taxonomy":['Bacteria', cfg.genus]}  TODO: taxonomy
            if 'name' not in contig:
                if contig['complete'] == True:
                    description = f'{genome["taxon"]} {contig["type"]}, complete sequence'
                else:
                    description = f'{genome["taxon"]} {contig["type"]}'
            else:
                if contig['type'] == bc.REPLICON_PLASMID and contig['complete'] == True:
                    description = f'{genome["taxon"]} plasmid {contig["name"]}, complete sequence'
                elif contig['type'] == bc.REPLICON_PLASMID and contig['complete'] == False:
                    description = f'{genome["taxon"]} plasmid {contig["name"]}'
                elif contig['type'] != bc.REPLICON_PLASMID and contig['complete'] == True:
                    description = f'{genome["taxon"]} {contig["name"]}, complete sequence'
                else:
                    description = f'{genome["taxon"]} {contig["name"]}'

            contig_rec = SeqIO.SeqRecord(id='.', name=contig['id'], description=description, annotations=contig_annotations, seq=Seq(contig['sequence']))

            seq_feature_list = []
            source_qualifiers = {'organism': genome['taxon'],
                                 'molecule_type': 'DNA',
                                 'strain': cfg.strain}
            source = SeqFeature(FeatureLocation(0, contig['length'], strand=+1), type='source', qualifiers=source_qualifiers)
            seq_feature_list.append(source)

            for feature in features:
                if feature['contig'] == contig['id']:
                    feature_type = feature['type']
                    insdc_feature_type = None
                    qualifier = {}
                    if (feature_type == bc.FEATURE_ORIC) or (feature_type == bc.FEATURE_ORIT) or (feature_type == bc.FEATURE_ORIV):
                        #TODO: Add fuzzy positions for oriCTV
                        qualifier['strand'] = feature['strand']
                        if (feature_type == bc.FEATURE_ORIT):
                            insdc_feature_type = bc.INSDC_FEATURE_ORIGIN_TRANSFER
                        else:
                            insdc_feature_type = bc.INSDC_FEATURE_ORIGIN_REPLICTION
                    else:
                        qualifier['db_xref'] = feature['db_xrefs']
                        qualifier['product'] = feature['product']
                        if (feature['type'] == bc.FEATURE_CDS) or (feature_type == bc.FEATURE_SORF):
                            qualifier['locus_tag'] = feature['locus']
                            qualifier['translation'] = feature['sequence']
                            qualifier['codon_start'] = 1
                            insdc_feature_type = bc.INSDC_FEATURE_CDS
                            inference = []
                            if (feature['type'] == bc.FEATURE_CDS):
                                inference.append('ab initio prediction:Prodigal:2.6')
                            else:
                                inference.append('ab initio prediction:bakta')
                            if 'hypothetical' not in feature:
                                if 'ups' in feature:
                                    if 'ncbi_nrp_id' in feature['ups']:
                                        protein_id = feature['ups']['ncbi_nrp_id']
                                        qualifier['protein_id'] = protein_id
                                if 'psc' in feature:
                                    if 'uniref90_id' in feature['psc']:
                                        psc_subject_id = feature['psc']['uniref90_id']
                                        inference.append('similar to AA sequence:UniProtKB:%s' % psc_subject_id)
                                if 'ips' in feature:
                                    if 'uniref100_id' in feature['ips']:
                                        ips_subject_id = feature['ips']['uniref100_id']
                                        if ips_subject_id != psc_subject_id:
                                            inference.append('similar to AA sequence:UniProtKB:%s' % ips_subject_id)
                            qualifier['inference'] = inference
                        elif (feature_type == bc.FEATURE_T_RNA):
                            # TODO: Position anticodon
                            qualifier['locus_tag'] = feature['locus']
                            if 'notes' in feature:
                                if 'anti_codon' in feature:
                                    qualifier['note'] = feature['notes']
                                    qualifier['anticodon'] = "(aa:" + str(feature['notes'])[-10:-7] + ",seq:" + str(feature['anti_codon']).lower() + ")"
                                else:
                                    qualifier['anticodon'] = "(seq:" + str(feature['anti_codon']).lower() + ")"
                            qualifier['inference'] = 'profile:tRNAscan:2.0'
                            insdc_feature_type = bc.INSDC_FEATURE_T_RNA
                        elif (feature_type == bc.FEATURE_TM_RNA):
                            qualifier['locus_tag'] = feature['locus']
                            qualifier['inference'] = 'profile:aragorn:1.2'
                            insdc_feature_type = bc.INSDC_FEATURE_TM_RNA
                        elif (feature_type == bc.FEATURE_R_RNA):
                            for reference in feature['db_xrefs']:
                                if reference.split(":")[0]=="RFAM":
                                    r_subject_id = reference.split(":")[1]
                            qualifier['locus_tag'] = feature['locus']
                            qualifier['inference'] = 'profile:Rfam:%s' % r_subject_id
                            insdc_feature_type = bc.INSDC_FEATURE_R_RNA
                        elif (feature_type == bc.FEATURE_NC_RNA):
                            # TODO: ncRNA_class
                            for reference in feature['db_xrefs']:
                                if reference.split(":")[0]=="RFAM":
                                    nc_subject_id = reference.split(":")[1]
                            qualifier['locus_tag'] = feature['locus']
                            qualifier['ncRNA_class'] = 'other'
                            qualifier['inference'] = 'profile:Rfam:%s' % nc_subject_id
                            insdc_feature_type = bc.INSDC_FEATURE_NC_RNA
                        elif (feature_type==bc.FEATURE_NC_RNA_REGION):
                            for reference in feature['db_xrefs']:
                                if reference.split(":")[0]=="RFAM":
                                    nc_subject_id = reference.split(":")[1]
                            qualifier['locus_tag'] = feature['locus']
                            qualifier['ncRNA_class'] = 'other'
                            qualifier['inference'] = 'profile:Rfam:%s' % nc_subject_id
                            insdc_feature_type = bc.INSDC_FEATURE_REGULATORY
                        elif (feature_type==bc.FEATURE_CRISPR):
                            qualifier['repeats'] = feature['repeats']
                            qualifier['repeat_consensus'] = feature['repeat_consensus']
                            qualifier['repeat_length'] = feature['repeat_length']
                            qualifier['spacer_length'] = feature['spacer_length']
                            feature['type'] = 'misc_feature'
                            insdc_feature_type = bc.INSDC_FEATURE_MISC_FEATURE

                    if feature['strand'] == bc.STRAND_REVERSE:
                        strand_dir = -1
                    else:
                        strand_dir = +1

                    if 'gene' in feature and feature['type'] != bc.FEATURE_NC_RNA_REGION:
                        qualifier['gene']=feature['gene']
                        gene_qualifier = {'gene': feature['gene'], 'locus_tag': feature['locus']}
                        gen_seqfeat = SeqFeature(FeatureLocation(feature['start']-1, feature['stop'], strand=strand_dir), type='gene', qualifiers=gene_qualifier)
                        seq_feature_list.append(gen_seqfeat)
                    feat_seqfeat = SeqFeature(FeatureLocation(feature['start']-1, feature['stop'], strand=strand_dir), type=insdc_feature_type, qualifiers=qualifier)

                    seq_feature_list.append(feat_seqfeat)
            contig_rec.features = seq_feature_list
            contig_list.append(contig_rec)

        SeqIO.write(contig_list, genbank_file, format="genbank")

        return
