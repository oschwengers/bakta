import logging
import re

from datetime import date, datetime
from pathlib import Path
from typing import Sequence, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation, AfterPosition, BeforePosition

import bakta
import bakta.config as cfg
import bakta.constants as bc
import bakta.so as so


log = logging.getLogger('INSDC')


def write_insdc(genome: dict, features:Sequence[dict], genbank_output_path: Path, embl_output_path: Path):
    log.debug('prepare: genbank=%s, embl=%s', genbank_output_path, embl_output_path)

    contig_list = []
    for contig in genome['contigs']:
        contig_features = [feat for feat in features if feat['contig'] == contig['id']]
        comment = (
            'Annotated with Bakta',
            f"Software: v{bakta.__version__}\n",
            f"Database: v{cfg.db_info['major']}.{cfg.db_info['minor']}\n",
            f'DOI: {bc.BAKTA_DOI}\n',
            f'URL: {bc.BAKTA_URL}\n',
            '\n',
            '##Genome Annotation Summary:##\n',
            f"{'Annotation Date':<30} :: {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}\n",
            f"{'Annotation Pipeline':<30} :: Bakta\n",
            f"{'Annotation Software version':<30} ::  v{bakta.__version__}\n",
            f"{'Annotation Database version':<30} ::  v{cfg.db_info['major']}.{cfg.db_info['minor']}\n",
            f"{'CDSs':<30} :: {len([feat for feat in contig_features if feat['type'] == bc.FEATURE_CDS or feat['type'] == bc.FEATURE_SORF]):5,}\n",
            f"{'tRNAs':<30} :: {len([feat for feat in contig_features if feat['type'] == bc.FEATURE_T_RNA]):5,}\n",
            f"{'tmRNAs':<30} :: {len([feat for feat in contig_features if feat['type'] == bc.FEATURE_TM_RNA]):5,}\n",
            f"{'rRNAs':<30} :: {len([feat for feat in contig_features if feat['type'] == bc.FEATURE_R_RNA]):5,}\n",
            f"{'ncRNAs':<30} :: {len([feat for feat in contig_features if feat['type'] == bc.FEATURE_NC_RNA]):5,}\n",
            f"{'regulatory ncRNAs':<30} :: {len([feat for feat in contig_features if feat['type'] == bc.FEATURE_NC_RNA_REGION]):5,}\n",
            f"{'CRISPR Arrays':<30} :: {len([feat for feat in contig_features if feat['type'] == bc.FEATURE_CRISPR]):5,}",
            f"{'oriCs/oriVs':<30} :: {len([feat for feat in contig_features if feat['type'] == bc.FEATURE_ORIC or feat['type'] == bc.FEATURE_ORIV]):5,}",
            f"{'oriTs':<30} :: {len([feat for feat in contig_features if feat['type'] == bc.FEATURE_ORIT]):5,}",
            f"{'gaps':<30} :: {len([feat for feat in contig_features if feat['type'] == bc.FEATURE_GAP]):5,}",
        )
        contig_annotations = {
            'molecule_type': 'DNA',
            'source': genome['taxon'],
            'date': date.today().strftime('%d-%b-%Y').upper(),
            'topology': contig['topology'],
            'data_file_division': 'HGT' if contig['type'] == bc.REPLICON_CONTIG else 'BCT',
            # 'accession': '*',  # hold back unil EMBL output bug is fixed in BioPython (https://github.com/biopython/biopython/pull/3572)
            'comment': comment
            # TODO: taxonomy
        }
        source_qualifiers = {
            'mol_type': 'genomic DNA'
            # 'molecule_type': 'DNA' #  might be necessary in BioPython > 1.78 along with removal of Seq(..., generic_dna)
        }

        description = ''
        if(genome['taxon']):
            contig_annotations['organism'] = genome['taxon']
            source_qualifiers['organism'] = genome['taxon']
            description = genome['taxon']
        if(genome['strain']):
            source_qualifiers['strain'] = genome['strain']

        if(contig['type'] == bc.REPLICON_PLASMID):
            source_qualifiers['plasmid'] = contig['name'] if contig.get('name', None) else 'unnamed'
            description = f"{description} plasmid {contig.get('name', 'unnamed')}"
            description += ', complete sequence' if contig['complete'] else ', whole genome shotgun sequence'
        elif(contig['type'] == bc.REPLICON_CHROMOSOME):
            if contig.get('name', None):
                source_qualifiers['chromosome'] = contig['name']
            description = f'{description} chromosome, complete genome' if contig['complete'] else f"{description} chromosome {contig['id']}, whole genome shotgun sequence"
        else:
            description += f" {contig['id']}, whole genome shotgun sequence"

        if(len(description) > 0 and description[0] == ' '):  # discard potential leading whitespace
            description = description[1:]

        contig_rec = SeqIO.SeqRecord(id=contig['id'], name=contig['id'], description=description, annotations=contig_annotations, seq=Seq(contig['sequence']))

        source = SeqFeature(FeatureLocation(0, contig['length'], strand=+1), type='source', qualifiers=source_qualifiers)
        seq_feature_list = [source]

        for feature in contig_features:
            insdc_feature_type = None
            qualifiers = {
                'note': []
            }
            if('db_xrefs' in feature):
                if(cfg.compliant):
                    qualifiers['db_xref'], qualifiers['note'] = revise_dbxref_insdc(feature['db_xrefs'])
                else:
                    qualifiers['db_xref'] = feature['db_xrefs']
            if('product' in feature):
                qualifiers['product'] = feature['product']
            if('locus' in feature):
                qualifiers['locus_tag'] = feature['locus']

            accompanying_features=[]
            if(feature['type'] == bc.FEATURE_GAP):
                insdc_feature_type = bc.INSDC_FEATURE_GAP
                qualifiers['estimated_length'] = feature['length']
            elif(feature['type'] == bc.FEATURE_ORIC or feature['type'] == bc.FEATURE_ORIV):
                # TODO: Add fuzzy positions for oriC/oriV
                insdc_feature_type = bc.INSDC_FEATURE_ORIGIN_REPLICATION
                qualifiers['inference'] = 'similar to DNA sequence'
                qualifiers['note'] = feature['product']
                if('product' in qualifiers):
                    qualifiers['note'] = feature['product']
                    del qualifiers['product']
            elif(feature['type'] == bc.FEATURE_ORIT):
                # TODO: Add fuzzy positions for oriT
                insdc_feature_type = bc.INSDC_FEATURE_ORIGIN_TRANSFER
                qualifiers['inference'] = 'similar to DNA sequence'
                qualifiers['note'] = feature['product']
                if('product' in qualifiers):
                    del qualifiers['product']
            elif(feature['type'] == bc.FEATURE_CDS) or (feature['type'] == bc.FEATURE_SORF):
                qualifiers['translation'] = feature['aa']
                qualifiers['codon_start'] = 1
                qualifiers['transl_table'] = cfg.translation_table
                insdc_feature_type = bc.INSDC_FEATURE_CDS
                inference = []
                inference.append('ab initio prediction:Prodigal:2.6' if feature['type'] == bc.FEATURE_CDS else f"ab initio prediction:Bakta:{'.'.join(bakta.__version__.split('.')[0:2])}")
                qualifiers['protein_id'] = f"gnl|Bakta|{feature['locus']}"
                if('ncbi_nrp_id' in feature.get('ups', {})):
                    nrp_id = feature['ups']['ncbi_nrp_id']
                    inference.append(f'similar to AA sequence:{bc.DB_XREF_REFSEQ_NRP}:{nrp_id}')
                elif('uniref100_id' in feature.get('ips', {})):
                    ips_subject_id = feature['ips']['uniref100_id']
                    inference.append(f'similar to AA sequence:{bc.DB_XREF_UNIREF}:{ips_subject_id}')
                elif('uniparc_id' in feature.get('ups', {})):
                    uniparc_id = feature['ups']['uniparc_id']
                    inference.append(f'similar to AA sequence:{bc.DB_XREF_UNIPARC}:{uniparc_id}')
                elif('uniref90_id' in feature.get('psc', {}) and feature.get('psc', {}).get('valid', False)):
                    psc_subject_id = feature['psc']['uniref90_id']
                    inference.append(f'similar to AA sequence:{bc.DB_XREF_UNIREF}:{psc_subject_id}')
                elif('uniref50_id' in feature.get('pscc', {})):
                    pscc_subject_id = feature['psc']['uniref50_id']
                    inference.append(f'similar to AA sequence:{bc.DB_XREF_UNIREF}:{pscc_subject_id}')
                qualifiers['inference'] = inference
                if(cfg.compliant):
                    qualifiers['note'], ec_number = extract_ec_from_notes_insdc(qualifiers, 'note')
                    if(ec_number is not None):
                            qualifiers['EC_number'] = ec_number
                if('exception' in feature):
                    ex = feature['exception']
                    pos = f"{ex['start']}..{ex['stop']}"
                    if(feature['strand'] == bc.STRAND_REVERSE):
                        pos = f"complement({pos})"
                    qualifiers['transl_except']=f"(pos:{pos},aa:{ex['aa']})"
                    qualifiers['note'].append(f"codon on position {ex['codon_position']} is a {ex['type']} codon")
                if(bc.FEATURE_SIGNAL_PEPTIDE in feature):
                    sigpep_qualifiers = {}
                    sigpep_qualifiers['locus_tag'] = feature['locus']
                    sigpep_qualifiers['inference'] = 'ab initio prediction:DeepSig:1.2'
                    sigpep = feature[bc.FEATURE_SIGNAL_PEPTIDE]
                    sigpep_strand = 1 if feature['strand'] == bc.STRAND_FORWARD else -1 if feature['strand'] == bc.STRAND_REVERSE else 0
                    sigpep_location = FeatureLocation(sigpep['start'] - 1, sigpep['stop'], strand = sigpep_strand)
                    sigpep_feature = SeqFeature(sigpep_location, type=bc.INSDC_FEATURE_SIGNAL_PEPTIDE, qualifiers=sigpep_qualifiers)
                    accompanying_features.append(sigpep_feature)
            elif(feature['type'] == bc.FEATURE_T_RNA):
                if('amino_acid' in feature and 'anti_codon' in feature):
                    if('anti_codon_pos' in feature):
                        anti_codon_pos = feature['anti_codon_pos']
                        qualifiers['anticodon'] = f"(pos:{anti_codon_pos[0]}..{anti_codon_pos[1]},aa:{feature['amino_acid']},seq:{feature['anti_codon']})"
                    else:
                        qualifiers['note'] = f"tRNA-{feature['amino_acid']} ({feature['anti_codon']})"
                qualifiers['inference'] = 'profile:tRNAscan:2.0'
                insdc_feature_type = bc.INSDC_FEATURE_T_RNA
                if('pseudo' in feature):
                    qualifiers['pseudo'] = None
            elif(feature['type'] == bc.FEATURE_TM_RNA):
                qualifiers['inference'] = 'profile:aragorn:1.2'
                insdc_feature_type = bc.INSDC_FEATURE_TM_RNA
            elif(feature['type'] == bc.FEATURE_R_RNA):
                for rfam_id in [dbxref.split(':')[1] for dbxref in feature['db_xrefs'] if dbxref.split(':')[0] == bc.DB_XREF_RFAM]:
                    qualifiers['inference'] = f'profile:Rfam:{rfam_id}'
                insdc_feature_type = bc.INSDC_FEATURE_R_RNA
            elif(feature['type'] == bc.FEATURE_NC_RNA):
                for rfam_id in [dbxref.split(':')[1] for dbxref in feature['db_xrefs'] if dbxref.split(':')[0] == bc.DB_XREF_RFAM]:
                    qualifiers['inference'] = f'profile:Rfam:{rfam_id}'
                qualifiers[bc.INSDC_FEATURE_NC_RNA_CLASS] = select_ncrna_class(feature)
                insdc_feature_type = bc.INSDC_FEATURE_NC_RNA
            elif(feature['type'] == bc.FEATURE_NC_RNA_REGION):
                for rfam_id in [dbxref.split(':')[1] for dbxref in feature['db_xrefs'] if dbxref.split(':')[0] == bc.DB_XREF_RFAM]:
                    qualifiers['inference'] = f'profile:Rfam:{rfam_id}'
                qualifiers[bc.INSDC_FEATURE_REGULATORY_CLASS] = select_regulatory_class(feature)
                insdc_feature_type = bc.INSDC_FEATURE_REGULATORY
                qualifiers['note'] = feature['product']
                qualifiers.pop('product', None)
            elif(feature['type'] == bc.FEATURE_CRISPR):
                qualifiers[bc.INSDC_FEATURE_REPEAT_FAMILY] = 'CRISPR'
                qualifiers[bc.INSDC_FEATURE_REPEAT_TYPE] = 'direct'
                qualifiers[bc.INSDC_FEATURE_REPEAT_UNIT_SEQ] = feature['repeat_consensus']
                qualifiers['inference'] = 'COORDINATES:alignment:pilercr:1.02'
                insdc_feature_type = bc.INSDC_FEATURE_REPEAT_REGION
                qualifiers['note'] = feature['product']
                qualifiers.pop('product', None)

            strand = None
            if(feature['strand'] == bc.STRAND_FORWARD):
                strand = 1
            elif(feature['strand'] == bc.STRAND_REVERSE):
                strand = -1
            elif(feature['strand'] == bc.STRAND_UNKNOWN):
                strand = 0

            start = feature['start'] - 1
            stop = feature['stop']
            if('edge' in feature):
                fl_1 = FeatureLocation(start, contig['length'], strand=strand)
                fl_2 = FeatureLocation(0, stop, strand=strand)
                if(feature['strand'] == bc.STRAND_REVERSE):
                    feature_location = CompoundLocation([fl_2, fl_1])
                else:
                    feature_location = CompoundLocation([fl_1, fl_2])
            else:
                if('truncated' in feature):
                    if(feature['truncated'] == bc.FEATURE_END_5_PRIME):
                        if(feature['strand'] == bc.STRAND_FORWARD):
                            start = BeforePosition(start)
                        else:
                            stop = AfterPosition(stop)
                    elif(feature['truncated'] == bc.FEATURE_END_3_PRIME):
                        if(feature['strand'] == bc.STRAND_FORWARD):
                            stop = AfterPosition(stop)
                        else:
                            start = BeforePosition(start)
                    else:
                        start = BeforePosition(start)
                        stop = AfterPosition(stop)
                feature_location = FeatureLocation(start, stop, strand=strand)
            if(feature.get('locus', None)):
                gene_qualifier = {
                    'locus_tag': feature['locus']
                }
                if(feature.get('gene', None)):
                    qualifiers['gene'] = feature['gene']
                    gene_qualifier['gene'] = feature['gene']
                gen_seqfeat = SeqFeature(feature_location, type='gene', qualifiers=gene_qualifier)
                seq_feature_list.append(gen_seqfeat)
            feat_seqfeat = SeqFeature(feature_location, type=insdc_feature_type, qualifiers=qualifiers)
            seq_feature_list.append(feat_seqfeat)
            for acc_feature in accompanying_features:  # add accompanying features, e.g. signal peptides
                seq_feature_list.append(acc_feature)
        contig_rec.features = seq_feature_list
        contig_list.append(contig_rec)

    with genbank_output_path.open('wt', encoding='utf-8') as fh:
        log.info('write GenBank: path=%s', genbank_output_path)
        SeqIO.write(contig_list, fh, format='genbank')

    with embl_output_path.open('wt', encoding='utf-8') as fh:
        log.info('write EMBL: path=%s', embl_output_path)
        SeqIO.write(contig_list, fh, format='embl')


def select_ncrna_class(feature: dict) -> str:
    if(feature['class'] is None):
        return bc.INSDC_FEATURE_NC_RNA_CLASS_OTHER
    elif(feature['class'].id == so.SO_NCRNA_GENE_ANTISENSE.id):
        return bc.INSDC_FEATURE_NC_RNA_CLASS_ANTISENSE
    elif(feature['class'].id == so.SO_NCRNA_GENE_RIBOZYME.id):
        return bc.INSDC_FEATURE_NC_RNA_CLASS_RIBOZYME
    elif(feature['class'].id == so.SO_NCRNA_GENE_RNASEP.id):
        return bc.INSDC_FEATURE_NC_RNA_CLASS_RNASEP
    else:
        return bc.INSDC_FEATURE_NC_RNA_CLASS_OTHER


def select_regulatory_class(feature: dict) -> str:
    if(feature['class'] is None):
        return bc.INSDC_FEATURE_REGULATORY_CLASS_OTHER
    elif(feature['class'].id == so.SO_CIS_REG_ATTENUATOR.id):
        return bc.INSDC_FEATURE_REGULATORY_CLASS_ATTENUATOR
    elif(feature['class'].id == so.SO_CIS_REG_RIBOSWITCH.id):
        return bc.INSDC_FEATURE_REGULATORY_CLASS_RIBOSWITCH
    elif(feature['class'].id == so.SO_CIS_REG_THERMOMETER.id):
        return bc.INSDC_FEATURE_REGULATORY_CLASS_RESPONSE_ELEMENT
    elif(feature['class'].id == so.SO_CIS_REG_RECODING_STIMULATION_REGION.id or feature['class'].id == so.SO_CIS_REG_FRAMESHIFT.id):
        return bc.INSDC_FEATURE_REGULATORY_CLASS_RECODING_STIMULATORY_REGION
    elif(feature['class'].id == so.SO_CIS_REG_RIBOSOME_BINDING_SITE.id):
        return bc.INSDC_FEATURE_REGULATORY_CLASS_RIBOSOME_BINDING_SITE
    else:
        return bc.INSDC_FEATURE_REGULATORY_CLASS_OTHER


def revise_product_insdc(feature: dict):
    """Revise product name for INSDC compliant submissions"""
    product = feature['product']

    old_product = product
    if(re.search(r'(uncharacteri[sz]ed)', product, flags=re.IGNORECASE)):  # replace putative synonyms)
        product = re.sub(r'(uncharacteri[sz]ed)', 'putative', product, flags=re.IGNORECASE)
        log.info('fix product: replace putative synonyms. new=%s, old=%s', product, old_product)

    old_product = product
    if(product.count('(') != product.count(')')):  # remove unbalanced parentheses
        product = product.replace('(', '').replace(')', '')  # ToDo: find and replace only legend parentheses
        log.info('fix product: remove unbalanced parantheses. new=%s, old=%s', product, old_product)

    old_product = product
    if(product.count('[') != product.count(']')):  # remove unbalanced brackets
        product = product.replace('[', '').replace(']', '')  # ToDo: find and replace only legend bracket
        log.info('fix product: remove unbalanced brackets. new=%s, old=%s', product, old_product)

    feature['product'] = product


def revise_dbxref_insdc(dbxrefs: Sequence[str]) -> Tuple[Sequence[str], Sequence[str]]:
    """Remove INSDC non-compliant DbXrefs."""
    insdc_valid_dbxrefs = [bc.DB_XREF_UNIPROTKB, bc.DB_XREF_GO, bc.DB_XREF_PFAM, bc.DB_XREF_RFAM]
    valid_dbxrefs = []
    invalid_dbxrefs = []
    for dbxref in dbxrefs:
        if(dbxref.split(':')[0] in insdc_valid_dbxrefs):
            valid_dbxrefs.append(dbxref)
        else:
            invalid_dbxrefs.append(dbxref)
    return valid_dbxrefs, invalid_dbxrefs


def extract_ec_from_notes_insdc(qualifiers: dict, note_key: str):
    note_tmp = []
    ec_number = None
    for note in qualifiers[note_key]:  # move EC numbers from note to EC_number qualifier
        if(note.split(':')[0] == bc.DB_XREF_EC):
            ec_number = note.replace('EC:', '')
        else:
            note_tmp.append(note)
    return note_tmp, ec_number
