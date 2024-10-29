import concurrent.futures as cf
import copy
import logging
import subprocess as sp
import sys
import xml.etree.ElementTree as ET

from collections import OrderedDict
from typing import Dict, Sequence, Set, Tuple, Union
from pathlib import Path

import pyrodigal
import pyhmmer

from Bio import SeqIO
from Bio import SeqFeature
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from xopen import xopen

import bakta.config as cfg
import bakta.constants as bc
import bakta.features.orf as orf
import bakta.utils as bu
import bakta.so as so

from bakta.psc import DB_PSC_COL_UNIREF90


log = logging.getLogger('CDS')


def predict(data: dict):
    """Predict open reading frames with Pyrodigal."""
    # create Pyrodigal trainining file if not provided by the user
    prodigal_tf_path = cfg.prodigal_tf
    trainings_info = None
    prodigal_metamode = cfg.meta  or  data['stats']['size'] < pyrodigal.MIN_SINGLE_GENOME  # 20_000 bp
    log.debug('prodigal mode: meta=%s', prodigal_metamode)
    if(prodigal_tf_path is None):
        closed = not data['genome']['complete']
        if(not prodigal_metamode):
            log.info('create prodigal training info object: meta=%s, closed=%s', prodigal_metamode, closed)
            gene_finder = pyrodigal.GeneFinder(meta=prodigal_metamode, closed=closed)
            seqs = [seq['nt'] for seq in data['sequences']]
            trainings_info = gene_finder.train(*seqs, translation_table=cfg.translation_table)
        else:
            log.info('skip creation of prodigal training info object: meta=%s, closed=%s', prodigal_metamode, closed)
    else:
        try:
            with prodigal_tf_path.open('rb') as fh_tf:
                trainings_info = pyrodigal.TrainingInfo.load(fh_tf)
                prodigal_metamode = False
                log.info('use provided prodigal training info: file=%s, meta mode=%s', prodigal_tf_path, prodigal_metamode)
        except:
            log.error('Cannot use provided prodigal training file: file=%s', prodigal_tf_path)
            sys.exit(f'Error! Cannot use provided prodigal training file: file={prodigal_tf_path}')

    cdss = []
    # predict genes on linear sequences
    linear_sequences = [seq for seq in data['sequences'] if seq['topology'] == bc.TOPOLOGY_LINEAR]
    if(len(linear_sequences) > 0):
        if prodigal_metamode:
            gene_finder = pyrodigal.GeneFinder(meta=True, metagenomic_bins=None, closed=True, mask=True)
        else:
            gene_finder = pyrodigal.GeneFinder(trainings_info, meta=False, closed=True, mask=True)
        sequences = [seq['nt'] for seq in linear_sequences]
        with cf.ThreadPoolExecutor(max_workers=cfg.threads) as tpe:
            for seq, genes in zip(linear_sequences, tpe.map(gene_finder.find_genes, sequences)):
                cdss_per_sequence = create_cdss(genes, seq)
                cdss.extend(cdss_per_sequence)

    # predict genes on circular replicons (chromosomes/plasmids)
    circular_sequences = [seq for seq in data['sequences'] if seq['topology'] == bc.TOPOLOGY_CIRCULAR]
    if(len(circular_sequences) > 0):
        if prodigal_metamode:
            gene_finder = pyrodigal.GeneFinder(meta=True, metagenomic_bins=None, closed=False, mask=True)
        else:
            gene_finder = pyrodigal.GeneFinder(trainings_info, meta=False, closed=False, mask=True)
        sequences = [seq['nt'] for seq in circular_sequences]
        with cf.ThreadPoolExecutor(max_workers=cfg.threads) as tpe:
            for seq, genes in zip(circular_sequences, tpe.map(gene_finder.find_genes, sequences)):
                cdss_per_sequence = create_cdss(genes, seq)
                cdss.extend(cdss_per_sequence)

    log.info('predicted=%i', len(cdss))
    return cdss


def create_cds(sequence: dict, start: int, stop: int, strand: str, edge:bool, nt: str, aa: str):
    cds = OrderedDict()
    cds['type'] = bc.FEATURE_CDS
    cds['sequence'] = sequence['id']
    cds['start'] = start
    cds['stop'] = stop
    cds['strand'] = strand
    cds['frame'] = (start - 1) % 3 + 1 if strand == bc.STRAND_FORWARD else (sequence['length'] - stop) % 3 + 1
    cds['gene'] = None
    cds['product'] = None
    cds['db_xrefs'] = [so.SO_CDS.id]
    cds['nt'] = nt
    cds['aa'] = aa
    cds['aa_digest'], cds['aa_hexdigest'] = bu.calc_aa_hash(aa)
    if(edge):
        cds['edge'] = True
    return cds


def create_cdss(genes, sequence):
    partial_cdss_per_sequence = []
    cdss_per_sequence = []
    for gene in genes:
        strand = bc.STRAND_FORWARD if gene.strand == 1 else bc.STRAND_REVERSE
        cds = create_cds(sequence, gene.begin, gene.end, strand, False, '', '')
        cds['start_type'] = gene.start_type
        cds['rbs_motif'] = gene.rbs_motif
        if gene.partial_begin  and  gene.partial_end:
            cds['truncated'] = bc.FEATURE_END_BOTH
            partial_cdss_per_sequence.append(cds)
        elif gene.partial_begin:
            cds['truncated'] = bc.FEATURE_END_5_PRIME if cds['strand'] == bc.STRAND_FORWARD else bc.FEATURE_END_3_PRIME
            partial_cdss_per_sequence.append(cds)
        elif gene.partial_end:
            cds['truncated'] = bc.FEATURE_END_3_PRIME if cds['strand'] == bc.STRAND_FORWARD else bc.FEATURE_END_5_PRIME
            partial_cdss_per_sequence.append(cds)
        else:
            cdss_per_sequence.append(cds)

        # In meta mode, let pyrodigal use the translation table the gene was predicted with
        if cfg.meta:
            aa = gene.translate(translation_table=None).upper()
        else:
            aa = gene.translate(translation_table=cfg.translation_table).upper()
        
        if(aa[-1] == '*'):
            aa = aa[:-1]  # discard trailing asterisk
        cds['aa'] = aa
        cds['aa_digest'], cds['aa_hexdigest'] = bu.calc_aa_hash(aa)
        
        log.info(
            'seq=%s, start=%i, stop=%i, strand=%s, frame=%s, truncated=%s, start-type=%s, RBS-motif=%s',
            cds['sequence'], cds['start'], cds['stop'], cds['strand'], cds['frame'], cds.get('truncated', 'no'), cds['start_type'], cds['rbs_motif']
        )
    if(sequence['topology'] == bc.TOPOLOGY_CIRCULAR and len(partial_cdss_per_sequence) >= 2):
        first_partial_cds = partial_cdss_per_sequence[0]  # first partial CDS per sequence
        last_partial_cds = partial_cdss_per_sequence[-1]  # last partial CDS per sequence
        # check if partial CDSs are on same strand and have opposite truncated edges
        # and first starts at 1 and last ends at sequence end (length)
        if(first_partial_cds['strand'] == last_partial_cds['strand']
            and first_partial_cds['truncated'] != last_partial_cds['truncated']
            and first_partial_cds['start'] == 1
            and last_partial_cds['stop'] == sequence['length']):
            cds = last_partial_cds
            cds['stop'] = first_partial_cds['stop']
            if(last_partial_cds['truncated'] == bc.FEATURE_END_3_PRIME):
                aa = last_partial_cds['aa'] + first_partial_cds['aa']  # merge sequence
            else:
                aa = first_partial_cds['aa'] + last_partial_cds['aa']  # merge sequence
                cds['start_type'] = first_partial_cds['start_type']
                cds['rbs_motif'] = first_partial_cds['rbs_motif']
            log.debug(f'trunc seq: seq-start={aa[:10]}, seq-end={aa[-10:]}')
            cds['edge'] = True  # mark CDS as edge feature
            cds.pop('truncated')
            cds['aa'] = aa
            cds['aa_digest'], cds['aa_hexdigest'] = bu.calc_aa_hash(aa)
            cdss_per_sequence.append(cds)
            log.info(
                'edge CDS: seq=%s, start=%i, stop=%i, strand=%s, frame=%s, start-type=%s, RBS-motif=%s, aa-hexdigest=%s, aa=[%s..%s]',
                cds['sequence'], cds['start'], cds['stop'], cds['strand'], cds['frame'], cds['start_type'], cds['rbs_motif'], cds['aa_hexdigest'], aa[:10], aa[-10:]
            )
            partial_cdss_per_sequence = partial_cdss_per_sequence[1:-1]  # remove first/last partial CDS
    for partial_cds in partial_cdss_per_sequence:
        cdss_per_sequence.append(partial_cds)
        log.info(
            'truncated CDS: seq=%s, start=%i, stop=%i, strand=%s, frame=%s, truncated=%s, start-type=%s, RBS-motif=%s, aa-hexdigest=%s, aa=[%s..%s]',
            partial_cds['sequence'], partial_cds['start'], partial_cds['stop'], partial_cds['strand'], partial_cds['frame'], partial_cds['truncated'], partial_cds['start_type'], partial_cds['rbs_motif'], partial_cds['aa_hexdigest'], partial_cds['aa'][:10], partial_cds['aa'][-10:]
        )
    for cds in cdss_per_sequence:  # extract nt sequences
        nt = bu.extract_feature_sequence(cds, sequence)
        cds['nt'] = nt
        log.info(
            'seq=%s, start=%i, stop=%i, strand=%s, nt=[%s..%s]',
            cds['sequence'], cds['start'], cds['stop'], cds['strand'], nt[:10], nt[-10:]
        )
    return cdss_per_sequence


def import_user_cdss(data: dict, import_path: Path):
    """Import user-provided CDS regions.
    Only CDS region information are imported skipping any existing functional annotations.
    
    Parameters
    ----------
    data : dict
        data dictionary holding sequence information
    import_path : Path
        Path to GFF3 or Genbank file with regions or features.

    Returns
    -------
    list
        a list of CDS features - without functional annotations.
    """
    user_cdss = []
    if(cfg.keep_sequence_headers):
        sequences_by_id = {seq['id']: seq for seq in data['sequences']}  # use ID as it's not altered -> no 'orig_id' field
    else:
        sequences_by_id = {seq['orig_id']: seq for seq in data['sequences']}  # use 'orig_id' instead of autogenerated new 'id'
    file_suffix = import_path.suffix.lower()
    if(file_suffix in ['.gff', '.gff3']):  # parse GFF3 format
        try:
            with xopen(str(import_path), threads=0) as fh_in:
                skip_lines = False
                for line in  fh_in:
                    line = line.strip()
                    if(line == '##FASTA'):
                        skip_lines = True
                    elif(skip_lines  or  line[0] == '#'):
                        continue
                    else:
                        sequence_id, tool, feature_type, start, stop, score, strand, phase, attributes = line.split('\t')
                        if(feature_type.lower() == 'cds'):
                            attributes = attributes.lower().split(';')
                            seq = sequences_by_id.get(sequence_id, None)
                            if(seq is None):
                                log.error('user-provided CDS: No seq found for id=%s', sequence_id)
                                raise Exception(f'user-provided CDS: No sequence found for id={sequence_id}')
                            edge = False
                            start = int(start)
                            stop = int(stop)
                            if(stop > seq['length']):  # check for features spanning sequence edges
                                stop = stop - seq['length']
                                edge = True
                                
                            user_cds = create_cds(seq, start, stop, strand, edge, '', '')
                            user_cds['source'] = bc.CDS_SOURCE_USER
                            if('pseudo=' in attributes  or  bc.INSDC_FEATURE_PSEUDOGENE in attributes):  # skip pseudo genes
                                log.debug(
                                    'skip user-provided CDS: reason=pseudogene seq=%s, start=%i, stop=%i, strand=%s',
                                    user_cds['sequence'], user_cds['start'], user_cds['stop'], user_cds['strand']
                                )
                                continue
                            try:
                                nt = bu.extract_feature_sequence(user_cds, seq)
                                user_cds['nt'] = nt
                            except:
                                log.error('user-provided CDS out of range! seq=%s, start=%i, stop=%i', user_cds['sequence'], user_cds['start'], user_cds['stop'])
                                raise ValueError(f"User-provided CDS out of range! sequence={user_cds['sequence']}, start={user_cds['start']}, stop={user_cds['stop']}")
                            try:
                                aa = str(Seq(nt).translate(table=cfg.translation_table, cds=True))
                                user_cds['aa'] = aa
                                user_cds['aa_digest'], user_cds['aa_hexdigest'] = bu.calc_aa_hash(aa)
                            except:
                                log.error('user-provided CDS could not be translated into a valid amino acid sequence! seq=%s, start=%i, stop=%i, cds=%s', user_cds['sequence'], user_cds['start'], user_cds['stop'], nt)
                                raise ValueError(f"User-provided CDS could not be translated into a valid amino acid sequence! sequence={user_cds['sequence']}, start={user_cds['start']}, stop={user_cds['stop']}, cds={nt}")
                            
                            log.info(
                                'user-provided CDS: seq=%s, start=%i, stop=%i, strand=%s, nt=[%s..%s], aa=[%s..%s]',
                                user_cds['sequence'], user_cds['start'], user_cds['stop'], user_cds['strand'], nt[:10], nt[-10:], aa[:10], aa[-10:]
                            )
                            user_cdss.append(user_cds)
        except Exception as e:
            log.error('user-provided regions/features file GFF3 format not valid!', exc_info=True)
            sys.exit(f'ERROR: User-provided regions/features file GFF3 format not valid!')
    elif(file_suffix in ['.gb', '.genbank', '.gbk', '.gbff']):  # parse GenBank format
        try:
            with xopen(str(import_path), threads=0) as fh_in:
                for record in SeqIO.parse(fh_in, 'genbank'):
                    for feature in record.features:
                        if(feature.type.lower() == 'cds'):
                            seq = sequences_by_id.get(record.id, None)
                            if(seq is None):
                                log.error('user-provided CDS: No seq found for id=%s', record.id)
                                raise Exception(f'user-provided CDS: No sequence found for id={record.id}')
                            if(feature.location.strand is None):  # weird mixed-stranded compound locations
                                strand = bc.STRAND_UNKNOWN
                            else:
                                strand = bc.STRAND_FORWARD if feature.location.strand == +1 else bc.STRAND_REVERSE
                            start = feature.location.start + 1
                            end = feature.location.end
                            edge = False
                            if('<' in str(feature.location.start)  or  '>' in str(feature.location.end)):
                                log.debug(
                                    'skip user-provided CDS: reason=partial, seq=%s, start=%s, stop=%s, strand=%s',
                                    seq['id'], feature.location.start, feature.location.end, strand
                                )
                                continue
                            elif(bc.INSDC_FEATURE_PSEUDO in feature.qualifiers  or  bc.INSDC_FEATURE_PSEUDOGENE in feature.qualifiers):
                                log.debug(
                                    'skip user-provided CDS: reason=pseudo, seq=%s, start=%i, stop=%i, strand=%s',
                                    seq['id'], feature.location.start, feature.location.end, strand
                                )
                                continue
                            elif('ribosomal_slippage' in feature.qualifiers):
                                log.debug(
                                    'skip user-provided CDS: reason=ribosomal slippage, seq=%s, start=%i, stop=%i, strand=%s',
                                    seq['id'], feature.location.start, feature.location.end, strand
                                )
                                continue
                            elif(isinstance(feature.location, SeqFeature.CompoundLocation)  and  len(feature.location.parts) == 2):
                                strand = feature.location.strand
                                if(strand != bc.STRAND_UNKNOWN):  # only accept equal strands -> edge feature 
                                    strand = bc.STRAND_FORWARD if feature.location.strand == +1 else bc.STRAND_REVERSE
                                    edge = True
                                    edge_left, edge_right = feature.location.parts
                                    if(strand == bc.STRAND_FORWARD):
                                        start = edge_left.start + 1
                                        end = edge_right.end
                                    else:
                                        start = edge_right.start + 1
                                        end = edge_left.end

                            user_cds = create_cds(seq, start, end, strand, edge, '', '')
                            user_cds['source'] = bc.CDS_SOURCE_USER
                            try:
                                nt = bu.extract_feature_sequence(user_cds, seq)
                                user_cds['nt'] = nt
                            except:
                                log.error('user-provided CDS: CDS out of range! seq=%s, start=%i, stop=%i', user_cds['sequence'], user_cds['start'], user_cds['stop'])
                                raise ValueError(f"User-provided CDS out of range! sequence={user_cds['sequence']}, start={user_cds['start']}, stop={user_cds['stop']}")
                            try:
                                aa = str(Seq(nt).translate(table=cfg.translation_table, cds=True))
                                user_cds['aa'] = aa
                                user_cds['aa_digest'], user_cds['aa_hexdigest'] = bu.calc_aa_hash(aa)
                            except:
                                log.error('user-provided CDS: CDS could not be translated into a valid amino acid sequence! seq=%s, start=%i, stop=%i, cds=%s', user_cds['sequence'], user_cds['start'], user_cds['stop'], nt)
                                raise ValueError(f"User-provided CDS could not be translated into a valid amino acid sequence! sequence={user_cds['sequence']}, start={user_cds['start']}, stop={user_cds['stop']}, cds={nt}")
                            
                            log.info(
                                'user-provided CDS: seq=%s, start=%i, stop=%i, strand=%s, nt=[%s..%s], aa=[%s..%s]',
                                user_cds['sequence'], user_cds['start'], user_cds['stop'], user_cds['strand'], nt[:10], nt[-10:], aa[:10], aa[-10:]
                            )
                            user_cdss.append(user_cds)
        except Exception as e:
            log.error('user-provided CDS: regions/features file GenBank format not valid!', exc_info=True)
            sys.exit(f'ERROR: User-provided regions/features file GenBank format not valid!')
    else:
        log.warn('user-provided regions/features file suffix not detected! suffix=%s, path=%s', file_suffix, str(import_path))
    
    return user_cdss


def predict_pfam(cdss: Sequence[dict]) -> Sequence[dict]:
    """Detect Pfam-A entries"""
    pfam_hits = []
    cds_with_pfams_hits = {}
    orf_by_aa_digest = orf.get_orf_dictionary(cdss)
    alphabet: pyhmmer.easel.Alphabet = pyhmmer.easel.Alphabet.amino()
    proteins: list[pyhmmer.easel.DigitalSequence] = [ pyhmmer.easel.TextSequence(sequence=cds['aa'], name=bytes(orf.get_orf_key(cds), 'UTF-8')).digitize(alphabet) for cds in cdss ]
    with pyhmmer.plan7.HMMFile(cfg.db_path.joinpath('pfam')) as hmm:
        for top_hits in pyhmmer.hmmsearch(hmm, proteins, bit_cutoffs='gathering', cpus=cfg.threads):
            for hit in top_hits:
                aa_identifier = hit.name.decode()
                cds = orf_by_aa_digest[aa_identifier]
                domain_cov = (hit.best_domain.alignment.hmm_to - hit.best_domain.alignment.hmm_from + 1) / len(hit.best_domain.alignment.hmm_sequence)
                aa_cov = (hit.best_domain.alignment.target_to - hit.best_domain.alignment.target_from + 1) / len(cds['aa'])

                pfam = OrderedDict()
                pfam['id'] = hit.best_domain.alignment.hmm_accession.decode()
                pfam['name'] = hit.best_domain.alignment.hmm_name.decode()
                pfam['length'] = len(hit.best_domain.alignment.hmm_sequence)
                pfam['aa_cov'] = aa_cov
                pfam['hmm_cov'] = domain_cov
                pfam['evalue'] = hit.evalue
                pfam['score'] = hit.score
                pfam['start'] = hit.best_domain.alignment.target_from
                pfam['stop'] = hit.best_domain.alignment.target_to

                cds.setdefault('pfams', [])
                cds['pfams'].append(pfam)
                cds.setdefault('db_xrefs', [])
                cds['db_xrefs'].append(f"PFAM:{pfam['id']}")
                pfam_hits.append(cds)
                cds_with_pfams_hits[aa_identifier] = cds
                log.info(
                    'pfam detected: seq=%s, start=%i, stop=%i, strand=%s, pfam-id=%s, length=%i, aa-start=%i, aa-stop=%i, aa-cov=%1.1f, hmm-cov=%1.1f, evalue=%1.1e, bitscore=%1.1f, name=%s',
                    cds['sequence'], cds['start'], cds['stop'], cds['strand'], pfam['id'], pfam['length'], pfam['start'],
                    pfam['stop'], pfam['aa_cov'], pfam['hmm_cov'], pfam['evalue'], pfam['score'], pfam['name']
                )
    log.info('predicted-pfams=%i, CDS-w/-pfams=%i', len(pfam_hits), len(cds_with_pfams_hits))
    return cds_with_pfams_hits.values()


def analyze_proteins(cdss: Sequence[dict]):
    for cds in cdss:
        seq = ProteinAnalysis(cds['aa'])
        seq_stats = OrderedDict()
        try:
            seq_stats['molecular_weight'] = seq.molecular_weight()
        except:
            log.warning(
                'could not calc molecular weight! seq=%s, start=%i, stop=%i, strand=%s, frame=%s',
                cds['sequence'], cds['start'], cds['stop'], cds['strand'], cds['frame']
            )
            seq_stats['molecular_weight'] = None
        try:
            seq_stats['isoelectric_point'] = seq.isoelectric_point()
        except:
            log.warning(
                'could not calc isoelectric point! seq=%s, start=%i, stop=%i, strand=%s, frame=%s',
                cds['sequence'], cds['start'], cds['stop'], cds['strand'], cds['frame']
            )
            seq_stats['isoelectric_point'] = None
        cds['seq_stats'] = seq_stats


def revise_translational_exceptions(data: dict, cdss: Sequence[dict]):
    """
    Revise translational exceptions as for istance selenocystein proteins.
    """
    no_revised = 0
    if(len([feat for feat in data['features'] if feat['type'] == bc.FEATURE_NC_RNA_REGION]) == 0):  # check if ncRNA regions have been detected, otherwise skip analysis and return
        return no_revised

    sequences = {seq['id']: seq for seq in data['sequences']}
    # detect splitted orphan ORFs of selenocystein proteins that are subject to stop codon recoding.
    cdss_per_sequences = {seq['id']: [] for seq in data['sequences']}  # get CDS per sequence
    for cds in cdss:
        cdss_per_sequence = cdss_per_sequences[cds['sequence']]
        if('truncated' not in cds):  # exclude truncated CDS for now
            cdss_per_sequence.append(cds)
    cds_pairs_per_sequence = {seq['id']: [] for seq in data['sequences']}  # extract inframe primate CDS neighbouring pairs
    for id, cdss_per_sequence in cdss_per_sequences.items():
        cdss_per_sequence = sorted(cdss_per_sequence, key=lambda k: k['start'])
        for i in range(1, len(cdss_per_sequence)):
            cds_a = cdss_per_sequence[i-1]
            cds_b = cdss_per_sequence[i]
            strand = cds_a['strand']
            upstream_stop_codon = cds_a['nt'][-3:] if strand == bc.STRAND_FORWARD else cds_b['nt'][-3:]
            if(
                cds_a['strand'] == cds_b['strand'] and  # up- and downstream ORFs on the same strand
                cds_a['frame'] == cds_b['frame'] and  # up- and downstream ORFs on the same frame
                upstream_stop_codon == 'TGA' and  # tRNAScan-SE 2.0 only predicts tRNA-Sec with UCA anticodons, therefore we can only detect TGA stop codons
                (cds_b['start'] - cds_a['stop']) < 100):  # up- and downstream ORFs in close proximity
                cds_pairs = cds_pairs_per_sequence[cds_a['sequence']]
                cds_pairs.append((cds_a, cds_b))

    recoding_regions = [ncrna_region for ncrna_region in data['features'] if ncrna_region['type'] == bc.FEATURE_NC_RNA_REGION  and  ncrna_region['class'] == so.SO_CIS_REG_RECODING_STIMULATION_REGION]  #  Selenocysteine insertion sequences
    for recoding_region in recoding_regions:
        if('selenocysteine' in recoding_region.get('product', '').lower()):
            cds_pairs = cds_pairs_per_sequence[recoding_region['sequence']]
            for cds_a, cds_b in cds_pairs:  # find CDS pair around recoding region
                strand = cds_a['strand']
                if(
                    strand == recoding_region['strand'] and  # everything is on the same strand
                    cds_a['start'] < recoding_region['start'] and recoding_region['stop'] < cds_b['stop']):  # recoding region lies between up- and downstream ORFs
                    log.debug(
                        'selenocysteine recoding ncRNA/CDS pair detected: seq=%s, strand=%s, CDS-A=[%i...%i] (%s..%s), recoding-ie=[%i..%i], CDS-B=[%i...%i] (%s..%s)',
                        recoding_region['sequence'], recoding_region['strand'], cds_a['start'], cds_a['stop'], cds_a['nt'][:10], cds_a['nt'][-10:], recoding_region['start'], recoding_region['stop'], cds_b['start'], cds_b['stop'], cds_b['nt'][:10], cds_b['nt'][-10:]
                    )
                    seleno_cds = copy.deepcopy(cds_a)
                    seleno_cds['stop'] = cds_b['stop']
                    seleno_cds['rbs_motif'] = cds_a['rbs_motif'] if strand == bc.STRAND_FORWARD else cds_b['rbs_motif']
                    seq = sequences[seleno_cds['sequence']]
                    nt = bu.extract_feature_sequence(seleno_cds, seq)
                    seleno_cds['nt'] = nt
                    aa = str(Seq(nt).translate(table=cfg.translation_table, stop_symbol='*', to_stop=False, cds=False))
                    if(
                        aa[0] == 'M' and  # starts with M
                        aa[-1] == '*' and  # ends with stop *
                        aa[1:-1].count('*') == 1  # contains exactly 1 additional stop (*) somewhere in between
                        ):
                        aa = aa.replace('*', 'U', 1)  # replace internal stop codon by U -> selenocysteine
                        aa = aa[:-1]  # remove stop asterisk
                        seleno_cds['aa'] = aa
                        seleno_cds['aa_digest'], seleno_cds['aa_hexdigest'] = bu.calc_aa_hash(aa)
                        seleno_cds['exception'] = {
                            'type': 'selenocysteine',
                            'aa': 'Sec',
                            'start': cds_a['stop'] - 2 if strand == bc.STRAND_FORWARD else cds_b['start'],
                            'stop': cds_a['stop'] if strand == bc.STRAND_FORWARD else cds_b['start'] + 2,
                            'codon_position': aa.find('U') + 1
                        }                    
                        cdss.append(seleno_cds)
                        log.info(
                            'selenocysteine CDS detected: seq=%s, start=%i, stop=%i, strand=%s, frame=%i, exception=[%i..%i], nt=[%s..%s], aa=[%s..%s], aa-hexdigest=%s',
                            seleno_cds['sequence'], seleno_cds['start'], seleno_cds['stop'], seleno_cds['strand'], seleno_cds['frame'], seleno_cds['exception']['start'], seleno_cds['exception']['stop'], nt[:10], nt[-10:], aa[:10], aa[-10:], seleno_cds['aa_hexdigest']
                        )
                        discard = {  # mark CDS a/b as discarded
                            'type': bc.DISCARD_TYPE_RECODING,
                            'description': f"selenocysteine recoding at ({seleno_cds['exception']['start']}..{seleno_cds['exception']['stop']})"
                        }
                        cds_a['discarded'] = discard
                        cds_b['discarded'] = discard
                        no_revised += 1
                    else:
                        log.warning(
                            'spurious selenocysteine CDS detected: seq=%s, start=%i, stop=%i, strand=%s, frame=%i, nt=[%s], aa=[%s]',
                            seleno_cds['sequence'], seleno_cds['start'], seleno_cds['stop'], seleno_cds['strand'], seleno_cds['frame'], nt, aa
                        )
    return no_revised


def revise_special_cases_annotated(data: dict, cdss: Sequence[dict]):
    """
    Revise rare but known special cases as for istance supposedly truncated dnaA genes on rotated chromosome starts
    which often appear on re-annotated genomes.
    """
    
    sequences = {seq['id']: seq for seq in data['sequences']}
    # look for supposedly truncated dnaA genes on rotated chromosome starts: start=1, strand=+
    dnaA = None
    for cds in cdss:
        seq = sequences[cds['sequence']]
        if(
            seq['complete'] and
            cds['start'] == 1 and 
            cds['strand'] == bc.STRAND_FORWARD and 
            cds['start_type'] == 'Edge' and 
            cds['rbs_motif'] is None and
            ('dnaa' in cds['product'].lower().split() or cds['gene'] == 'dnaA')):
            dnaA = cds
            break
    if(dnaA is not None and 'truncated' in dnaA):
        dnaA.pop('truncated')
        gene = dnaA.get('gene', '-')
        log.info(
            'revise supposedly truncated dnaA gene on rotated chromosome start: seq=%s, start=%i, stop=%i, strand=%s, gene=%s, product=%s, nt=[%s..%s], aa=[%s..%s]',
            dnaA['sequence'], dnaA['start'], dnaA['stop'], dnaA['strand'], gene, dnaA['product'], dnaA['nt'][:10], dnaA['nt'][-10:], dnaA['aa'][:10], dnaA['aa'][-10:]
        )
    
    # look for supposedly truncated repA genes on rotated plasmid starts: start=1, strand=+
    repAs = []
    for cds in cdss:
        seq = sequences[cds['sequence']]
        if(
            seq['complete'] and
            cds['start'] == 1 and 
            cds['strand'] == bc.STRAND_FORWARD and 
            cds['start_type'] == 'Edge' and 
            cds['rbs_motif'] is None and
            ('repa' in cds['product'].lower().split() or cds['gene'] == 'repA')):
            repAs.append(cds)
    for repA in repAs:
        if('truncated' in repA):
            repA.pop('truncated')
            gene = repA.get('gene', '-')
            log.info(
                'revise supposedly truncated repA gene on rotated plasmid start: seq=%s, start=%i, stop=%i, strand=%s, gene=%s, product=%s, nt=[%s..%s], aa=[%s..%s]',
                repA['sequence'], repA['start'], repA['stop'], repA['strand'], gene, repA['product'], repA['nt'][:10], repA['nt'][-10:], repA['aa'][:10], repA['aa'][-10:]
            )


def predict_pseudo_candidates(hypotheticals: Sequence[dict]) -> Sequence[dict]:
    """
    Conduct homology search of hypothetical CDSs against the PSC db to find pseudogene candidates.
    """
    diamond_db_path = cfg.db_path.joinpath('psc.dmnd')
    diamond_output_path = cfg.tmp_path.joinpath('cds.pseudo.candidates.diamond.tsv')
    cds_hypotheticals_faa_path = cfg.tmp_path.joinpath('cds.pseudo.candidates.faa')
    orf.write_internal_faa(hypotheticals, cds_hypotheticals_faa_path)
    # TODO allow multiple hits
    cmd = [
        'diamond',
        'blastp',
        '--db', str(diamond_db_path),
        '--query', str(cds_hypotheticals_faa_path),
        '--out', str(diamond_output_path),
        '--id', str(int(bc.MIN_PSEUDOGENE_IDENTITY * 100)),                     # '80'
        '--query-cover', str(int(bc.MIN_PSEUDOGENE_QUERY_COVERAGE * 100)),      # '80'
        '--subject-cover', str(int(bc.MIN_PSEUDOGENE_SUBJECT_COVERAGE * 100)),  # '40'
        '--max-target-seqs', '1',  # single best output
        '--outfmt', '6', 'qseqid', 'sseqid', 'qlen', 'slen', 'length', 'pident', 'evalue', 'bitscore', 'qstart', 'qend', 'sstart', 'send', 'full_sseq',
        '--threads', str(cfg.threads),
        '--tmpdir', str(cfg.tmp_path),
        '--block-size', '3',  # slightly increase block size for faster executions
        '--fast'
    ]
    log.debug('cmd=%s', cmd)
    proc = sp.run(
        cmd,
        cwd=str(cfg.tmp_path),
        env=cfg.env,
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        universal_newlines=True
    )
    if proc.returncode != 0:  # exit-code 132 CORE DUMP -> CPU misses AVX extension
        log.debug('stdout=\'%s\', stderr=\'%s\'', proc.stdout, proc.stderr)
        log.warning('Diamond failed! diamond-error-code=%d', proc.returncode)
        raise Exception(f'diamond error! error code: {proc.returncode}')

    pseudo_candidates = []
    cds_by_hexdigest = orf.get_orf_dictionary(hypotheticals)
    with diamond_output_path.open() as fh:
        for line in fh:
            (aa_identifier, cluster_id, query_length, subject_length, alignment_length, identity, evalue, bitscore, query_start, query_end, subject_start, subject_end, subject_sequence) = line.rstrip('\n').split('\t')
            cds = cds_by_hexdigest[aa_identifier]
            query_cov = int(alignment_length) / len(cds['aa'])
            subject_cov = int(alignment_length) / int(subject_length)
            identity = float(identity) / 100
            bitscore = float(bitscore)
            evalue = float(evalue)
            if(query_cov >= bc.MIN_PSEUDOGENE_QUERY_COVERAGE and bc.MIN_PSEUDOGENE_SUBJECT_COVERAGE <= subject_cov < bc.MIN_PSC_COVERAGE and identity >= bc.MIN_PSEUDOGENE_IDENTITY):
                cds['pseudo-inference'] = {
                    DB_PSC_COL_UNIREF90: cluster_id,
                    'query_cov': query_cov,
                    'subject_cov': subject_cov,
                    'identity': identity,
                    'score': bitscore,
                    'evalue': evalue,
                    'gene_start': int(query_start),
                    'gene_end': int(query_end),
                    'reference_start': int(subject_start),
                    'reference_end': int(subject_end),
                    'reference_sequence': subject_sequence
                }
                pseudo_candidates.append(cds)
                log.debug(
                    'pseudogene-candidate: seq=%s, start=%i, stop=%i, strand=%s, aa-length=%i, query-cov=%0.3f, subject-cov=%0.3f, identity=%0.3f, score=%0.1f, evalue=%1.1e, UniRef90=%s',
                    cds['sequence'], cds['start'], cds['stop'], cds['strand'], len(cds['aa']), query_cov, subject_cov, identity, bitscore, evalue, cluster_id
                )
    log.info('found: pseudogene-candidates=%i', len(pseudo_candidates))
    return pseudo_candidates


def detect_pseudogenes(candidates: Sequence[dict], cdss: Sequence[dict], data: dict) -> Sequence[dict]:
    """
    Conduct a BLASTX search of 5'/3'-extended sequences of pseudogene candidates against matching PSCs.
    Search for and determine possible pseudogenization causes in the resulting alignments.
    """
    psc_references_faa_path = cfg.tmp_path.joinpath('cds.pseudo.references.faa')
    psc_references_dmnd_path = cfg.tmp_path.joinpath('cds.pseudo.references.dmnd')
    candidates_elongated_sequences_path = cfg.tmp_path.joinpath('cds.pseudo.elongated-sequences.fna')
    candidates_blastx_output_path = cfg.tmp_path.joinpath('cds.pseudo.blastx.xml')

    # Write unique PSC cluster sequences to a FAA file
    with psc_references_faa_path.open(mode='w') as fh:
        for cluster_id, faa_seq in {cds['pseudo-inference'][DB_PSC_COL_UNIREF90]: cds['pseudo-inference']['reference_sequence'] for cds in candidates}.items():
            fh.write(f">{cluster_id}\n{faa_seq}\n")

    # Get extended cds sequences
    sequences = {seq['id']: seq for seq in data['sequences']}
    candidates_extended_positions = {}
    with candidates_elongated_sequences_path.open(mode='w') as fh:
        for cds in candidates:
            seq = sequences[cds['sequence']]
            cds_elongated = get_elongated_cds(cds, seq)
            seq = bu.extract_feature_sequence(cds_elongated, seq)
            orf_key = orf.get_orf_key(cds)
            fh.write(f">{orf_key}\n{seq}\n")
            candidates_extended_positions[orf_key] = cds_elongated

    commands = [
        [
            'diamond',
            'makedb',
            '--in', str(psc_references_faa_path),
            '--db', str(psc_references_dmnd_path)
        ],
        [
            'diamond',
            'blastx',
            '--db', str(psc_references_dmnd_path),  # PSC pseudogene candidates
            '--query', str(candidates_elongated_sequences_path),  # nucleotide sequences of hypotheticals
            '--out', str(candidates_blastx_output_path),
            '--outfmt', '5',
            '--threads', str(cfg.threads),
            '--tmpdir', str(cfg.tmp_path),  # use tmp folder
            '--block-size', '3',  # slightly increase block size for faster executions
            '--query-gencode', str(cfg.translation_table),
            '--strand', 'plus',
            '--frameshift', '15',
            '--ultra-sensitive'
        ]
    ]

    for cmd in commands:
        log.debug('cmd=%s', cmd)
        proc = sp.run(
            cmd,
            cwd=str(cfg.tmp_path),
            env=cfg.env,
            stdout=sp.PIPE,
            stderr=sp.PIPE,
            universal_newlines=True
        )
        if proc.returncode != 0:
            log.debug('stdout=\'%s\', stderr=\'%s\'', proc.stdout, proc.stderr)
            log.warning('PSEUDO failed! diamond-error-code=%d', proc.returncode)
            raise Exception(f'diamond error! error code: {proc.returncode}\n{proc.stdout}')

    pseudogenes = []
    cds_by_hexdigest = orf.get_orf_dictionary(candidates)
    uniref90_by_hexdigest = {aa_identifier: cds['psc']['uniref90_id'] for aa_identifier, cds in orf.get_orf_dictionary(cdss).items() if 'psc' in cds and 'uniref90_id' in cds['psc']}
    with candidates_blastx_output_path.open() as fh:
        root = ET.parse(fh).getroot()
        for query in root.findall('./BlastOutput_iterations/Iteration'):
            aa_identifier = query.find('Iteration_query-def').text
            cds = cds_by_hexdigest[aa_identifier]
            extended_positions = candidates_extended_positions[aa_identifier]
            for hit in query.findall('./Iteration_hits/Hit'):
                cluster_id = hit.find('Hit_id').text
                if cluster_id == cds['pseudo-inference'][DB_PSC_COL_UNIREF90]:
                    query_alignment = hit.find('Hit_hsps/Hsp/Hsp_qseq').text
                    ref_alignment = hit.find('Hit_hsps/Hsp/Hsp_hseq').text
                    query_alignment_start = int(hit.find('Hit_hsps/Hsp/Hsp_query-from').text)
                    query_alignment_stop = int(hit.find('Hit_hsps/Hsp/Hsp_query-to').text)
                    alignment_length = int(hit.find('Hit_hsps/Hsp/Hsp_align-len').text)
                    identity = float(hit.find('Hit_hsps/Hsp/Hsp_identity').text) / alignment_length
                    bitscore = float(hit.find('Hit_hsps/Hsp/Hsp_bit-score').text)
                    evalue = float(hit.find('Hit_hsps/Hsp/Hsp_evalue').text)

                    if alignment_length == len(cds['aa']):  # skip non-extended genes (full match)
                        log.debug(
                            'no pseudogene (full match): seq=%s, start=%i, stop=%i, strand=%s',
                            cds['sequence'], cds['start'], cds['stop'], cds['strand']
                        )
                        continue

                    observations, positions = detect_pseudogenization_observations(
                        query_alignment,
                        ref_alignment,
                        query_alignment_start,
                        query_alignment_stop,
                        extended_positions,
                        cds
                    )

                    directions = observations.get('directions', [])
                    if bc.FEATURE_END_5_PRIME in directions or bc.FEATURE_END_3_PRIME in directions:
                        pseudogene = {
                            'start': positions['start'],
                            'stop': positions['stop'],
                            'observations': clean_observations(observations),
                            'inference': cds['pseudo-inference'],
                            'paralog': is_paralog(uniref90_by_hexdigest, aa_identifier, cluster_id),
                            'identity': identity,
                            'score': bitscore,
                            'evalue': evalue
                        }

                        effects = []
                        if len(observations.get(bc.PSEUDOGENE_EFFECT_START, [])) > 0:
                            start_codon = ', '.join(map(str, observations[bc.PSEUDOGENE_EFFECT_START]))
                            effects.append(f'Internal start codon at {start_codon}')
                        if len(observations.get(bc.PSEUDOGENE_EFFECT_STOP, [])) > 0:
                            stop_codon = ', '.join(map(str, observations[bc.PSEUDOGENE_EFFECT_STOP]))
                            effects.append(f'Internal stop codon at {stop_codon}')
                        effects = '; '.join(effects)

                        causes = []
                        if len(observations.get(bc.PSEUDOGENE_CAUSE_INSERTION, [])) > 0:
                            insertions = ', '.join(map(str, observations[bc.PSEUDOGENE_CAUSE_INSERTION]))
                            causes.append(f"Frameshift due to insertion around {insertions}.")
                        if len(observations.get(bc.PSEUDOGENE_CAUSE_DELETION, [])) > 0:
                            deletions = ', '.join(map(str, observations[bc.PSEUDOGENE_CAUSE_DELETION]))
                            causes.append(f"Frameshift due to deletion around {deletions}.")
                        if len(observations.get(bc.PSEUDOGENE_CAUSE_MUTATION, [])) > 0:
                            mutations = ', '.join(map(str, observations[bc.PSEUDOGENE_CAUSE_MUTATION]))
                            causes.append(f"Nonsense mutation around {mutations}.")
                        # if observations.get(bc.PSEUDOGENE_EXCEPTION_SELENOCYSTEINE, None):  # only for pseudogenes with translation exception + other cause
                        #     causes.append('Translation exception: Selenocysteine around ' + ', '.join(map(str, observations[bc.PSEUDOGENE_EXCEPTION_SELENOCYSTEINE])) + '.')
                        # if observations.get(bc.PSEUDOGENE_EXCEPTION_PYROLYSINE, None):  # only for pseudogenes with translation exception + other cause
                        #     causes.append('Translation exception: Pyrolysin around ' + ', '.join(map(str, observations[bc.PSEUDOGENE_EXCEPTION_PYROLYSINE])) + '.')
                        causes = ' '.join(causes)  # pseudogene cause
                        pseudogene['description'] = f"{effects}. {causes}" if effects != '' else causes

                        if bc.FEATURE_END_5_PRIME in directions and bc.FEATURE_END_3_PRIME in directions:
                            cds['truncated'] = bc.FEATURE_END_BOTH
                        elif bc.FEATURE_END_5_PRIME in directions:
                            cds['truncated'] = bc.FEATURE_END_5_PRIME if cds['strand'] == bc.STRAND_FORWARD else bc.FEATURE_END_3_PRIME
                        elif bc.FEATURE_END_3_PRIME in directions:
                            cds['truncated'] = bc.FEATURE_END_3_PRIME if cds['strand'] == bc.STRAND_FORWARD else bc.FEATURE_END_5_PRIME
                        cds[bc.PSEUDOGENE] = pseudogene
                        cds.pop('hypothetical')
                        pseudogenes.append(cds)
                        log.info(
                            'pseudogene: seq=%s, start=%i, stop=%i, strand=%s, insertions=%s, deletions=%s, mutations=%s, effect=%s',
                            cds['sequence'], cds['start'], cds['stop'], cds['strand'], observations.get(bc.PSEUDOGENE_CAUSE_INSERTION, []), observations.get(bc.PSEUDOGENE_CAUSE_DELETION, []), observations.get(bc.PSEUDOGENE_CAUSE_MUTATION, []), effects
                        )

                    elif observations[bc.PSEUDOGENE_EXCEPTION_SELENOCYSTEINE] or observations[bc.PSEUDOGENE_EXCEPTION_PYROLYSINE]:
                        # TODO handle translation exceptions, correct annotation
                        pass

    for cds in candidates:
        cds.pop('pseudo-inference')
    log.info('found: pseudogenes=%i', len(pseudogenes))
    return pseudogenes


def get_elongated_cds(cds: dict, sequence: dict, offset: int = bc.PSEUDOGENE_OFFSET) -> Dict[str, Union[int, str, bool]]:
    """
    Elongate the given CDS sequence with the offset in upstream and downstream direction, if possible.
    """
    elongated_cds = {
        'start': cds['start'],
        'stop': cds['stop'],
        'strand': cds['strand'],
        'edge': cds.get('edge', False),
        'elongation_upstream': offset,
        'elongation_downstream': offset
    }

    sequence_length = len(sequence['nt'])
    if sequence['topology'] == 'circular' and elongated_cds['start'] - offset < 0:
        elongated_cds['start'] = sequence_length + elongated_cds['start'] - offset
        elongated_cds['edge'] = True
    elif elongated_cds['start'] - offset < 0:
        elongated_cds['start'] = 1
        elongated_cds['elongation_upstream'] = cds['start']
    else:
        elongated_cds['start'] = elongated_cds['start'] - offset

    if sequence['topology'] == 'circular' and elongated_cds['stop'] + offset > sequence_length:
        elongated_cds['stop'] = elongated_cds['stop'] + offset - sequence_length
        elongated_cds['edge'] = True
    elif elongated_cds['stop'] + offset > sequence_length:
        elongated_cds['stop'] = sequence_length
        elongated_cds['elongation_downstream'] = sequence_length - cds['stop']
    else:
        elongated_cds['stop'] = elongated_cds['stop'] + offset

    return elongated_cds


def is_paralog(uniref90_by_hexdigest: Dict[str, str], aa_identifier: str, cluster: str) -> bool:
    """
    Check if the pseudogene is 'unprocessed' (the pseudogene has arisen from a copy of a parent gene by duplication
    followed by accumulation of random mutation).
    """
    excluded_uniref90 = set()
    for hexdigest, uniref90 in uniref90_by_hexdigest.items():
        if hexdigest != aa_identifier:
            excluded_uniref90.add(uniref90)
    return cluster in excluded_uniref90


def detect_pseudogenization_observations(alignment: str, ref_alignment: str, qstart: int, qstop: int, extended_positions: dict,
                            cds: dict) -> Tuple[Dict[str, Union[Set[int], bool]], Dict]:
    """
    Search for pseudogenization observations in the given alignments.
    """
    observations = {
        bc.PSEUDOGENE_CAUSE_INSERTION: set(),
        bc.PSEUDOGENE_CAUSE_DELETION: set(),
        bc.PSEUDOGENE_CAUSE_MUTATION: set(),
        bc.PSEUDOGENE_EFFECT_START: set(),
        bc.PSEUDOGENE_EFFECT_STOP: set(),
        bc.PSEUDOGENE_EXCEPTION_SELENOCYSTEINE: set(),
        bc.PSEUDOGENE_EXCEPTION_PYROLYSINE: set(),
        'directions': set()
    }

    elongated_edge = False
    if extended_positions.get('edge', False):
        # TODO implement edge case
        return observations, dict()

    if cds['strand'] == bc.STRAND_FORWARD:
        positions = {
            'upstream': -1 * (extended_positions['elongation_upstream'] - qstart + 1),  # index
            'downstream': (extended_positions['start'] + qstop) - cds['stop'] - 1 + 3  # index, stop codon
        }
        positions['start'] = get_abs_position(cds, cds['start'], positions['upstream'], elongated_edge)
        positions['stop'] = get_abs_position(cds, cds['stop'], positions['downstream'], elongated_edge)
    else:
        positions = {
            'upstream': -1 * (extended_positions['elongation_downstream'] - qstart + 1),   # index
            'downstream': qstop - (extended_positions['stop'] - cds['start']) - 1 + 3  # index, stop codon
        }
        positions['start'] = get_abs_position(cds, cds['start'], positions['downstream'], elongated_edge)
        positions['stop'] = get_abs_position(cds, cds['stop'], positions['upstream'], elongated_edge)

    # Skip pseudogenes with wrong blastx hsp hits
    # either cds['start'] or cds['stop'] are not within the pseudogene positions and less than 80% of the cds are covered
    if (not (positions['start'] <= cds['start'] + 3 <= positions['stop']) or not (positions['start'] <= cds['stop'] - 3 <= positions['stop'])) and \
            len(range(max(cds['start'], positions['start']), min(cds['stop'], positions['stop']) + 1)) < int(len(cds['nt']) * 0.8):
        return observations, positions

    compare_alignments(observations, alignment, ref_alignment, cds, positions, elongated_edge)

    return observations, positions


def compare_alignments(observations: dict, alignment: str, ref_alignment: str, cds: dict, positions: dict, edge: bool):
    """
    Compare the alignment and reference alignment to find the causes of pseudogenization.
    """
    alignment_position = positions['upstream']
    start = cds['start'] if cds['strand'] == bc.STRAND_FORWARD else cds['stop']
    insertions = 0
    deletions = 0

    # Check for anomalous start codons
    if positions['upstream'] < 0:
        up_length = (-1 * positions['upstream']) // 3
        if alignment[up_length-1] == 'M' and ref_alignment[up_length-1] != 'M':
            if cds['rbs_motif'] is None:  # point mutation -> internal start codon
                genome_position = get_abs_position(cds, start, alignment_position + (up_length * 3), edge)
                observations[bc.PSEUDOGENE_EFFECT_START].add(genome_position)
                observations['directions'].add(bc.FEATURE_END_3_PRIME)
                log.info(
                    'pseudogene observation: seq=%s, start=%i, stop=%i, strand=%s, original start=%i',
                    cds['sequence'], cds['start'], cds['stop'], cds['strand'], cds['start'] + genome_position
                )
            else:  # RBS was predicted (protein iso-form) -> skip
                pass

        # TODO implement case: point mutation -> loss of original start codon
        if (alignment[0] == 'M' and ref_alignment[0] == 'M' and up_length != 0) and \
                (alignment[up_length + 1] == 'M' and ref_alignment[up_length + 1] == 'M') and cds['rbs_motif'] is None:
            # TODO correct structural annotation
            pass

    # Check alignment
    for char, ref_char in zip(alignment, ref_alignment):
        if char == '-':
            continue
        elif char == '\\':  # insertion
            insertions += 1
            genome_position = get_abs_position(cds, start, alignment_position, edge)
            observations[bc.PSEUDOGENE_CAUSE_INSERTION].add(genome_position)
            observations['directions'].add(get_direction(alignment_position, edge))
            log.info(
                'pseudogene observation: seq=%s, start=%i, stop=%i, strand=%s, cause=insertion, position=%i',
                cds['sequence'], cds['start'], cds['stop'], cds['strand'], genome_position
            )
            alignment_position += 1
        elif char == '/':  # deletion
            deletions += 1
            genome_position = get_abs_position(cds, start, alignment_position, edge)
            observations[bc.PSEUDOGENE_CAUSE_DELETION].add(genome_position)
            observations['directions'].add(get_direction(alignment_position, edge))
            log.info(
                'pseudogene observation: seq=%s, start=%i, stop=%i, strand=%s, cause=deletion, position=%i',
                cds['sequence'], cds['start'], cds['stop'], cds['strand'], genome_position
            )
        elif char == '*':  # stop codon, selenocysteine, pyrolysine
            if ref_char == 'U':  # selenocysteine
                genome_position = get_abs_position(cds, start, alignment_position, edge)
                observations[bc.PSEUDOGENE_EXCEPTION_SELENOCYSTEINE].add(genome_position)
                log.info(
                    'pseudogene observation: seq=%s, start=%i, stop=%i, strand=%s, exception=selenocysteine, position=%i',
                    cds['sequence'], cds['start'], cds['stop'], cds['strand'], genome_position
                )
            elif ref_char == 'O':  # pyrolysine
                genome_position = get_abs_position(cds, start, alignment_position, edge)
                observations[bc.PSEUDOGENE_EXCEPTION_PYROLYSINE].add(genome_position)
                log.info(
                    'pseudogene observation: seq=%s, start=%i, stop=%i, strand=%s, exception=pyrolysin, position=%i',
                    cds['sequence'], cds['start'], cds['stop'], cds['strand'], genome_position
                )
            else:  # stop codon
                mutation = ''
                genome_position = get_abs_position(cds, start, alignment_position, edge)
                if abs(insertions - deletions) % 3 == 0:
                    observations[bc.PSEUDOGENE_CAUSE_MUTATION].add(genome_position)
                    mutation = ', cause=mutation'
                observations[bc.PSEUDOGENE_EFFECT_STOP].add(genome_position)
                observations['directions'].add(get_direction(alignment_position, edge))
                log.info(
                    'pseudogene observation: seq=%s, start=%i, stop=%i, strand=%s, effect=stop%s, position=%i',
                    cds['sequence'], cds['start'], cds['stop'], cds['strand'], mutation, genome_position
                )
            alignment_position += 3
        else:
            alignment_position += 3


def get_abs_position(cds, initial_position, movement, edge):
    """
    Return the absolute cds start position.
    """
    if edge:
        pass
    else:
        if cds['strand'] == bc.STRAND_FORWARD:
            return initial_position + movement
        else:
            return initial_position - movement


def get_direction(position: int, edge: bool) -> str:
    if edge:
        pass
    else:
        if position < 0:
            return bc.FEATURE_END_5_PRIME
        else:
            return bc.FEATURE_END_3_PRIME


def clean_observations(observations: dict) -> dict:
    """
    Cleanup pseudogenization observations.
    """
    drop_keys = []
    for key, value in observations.items():
        if len(value) == 0:
            drop_keys.append(key)
        else:
            observations[key] = sorted(list(value))
    for key in drop_keys:
        observations.pop(key)
    return observations
