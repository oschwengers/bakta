import copy
import logging
import subprocess as sp
import xml.etree.ElementTree as ET

from collections import OrderedDict
from math import ceil
from typing import Dict, Sequence, Set, Tuple, Union
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis

import bakta.config as cfg
import bakta.constants as bc
import bakta.features.orf as orf
import bakta.io.fasta as fasta
import bakta.utils as bu
import bakta.so as so

from bakta.psc import DB_PSC_COL_UNIREF90


log = logging.getLogger('CDS')


def predict(genome: dict, sequences_path: Path):
    """Predict open reading frames with Prodigal."""
    # create prodigal trainining file if not provided by the user
    prodigal_tf_path = cfg.prodigal_tf
    if(prodigal_tf_path is None):
        if(genome['size'] >= 20000):
            prodigal_tf_path = cfg.tmp_path.joinpath('prodigal.tf')
            log.info('create prodigal training file: file=%s', prodigal_tf_path)
            execute_prodigal(genome, sequences_path, prodigal_tf_path, train=True, complete=genome['complete'])
        else:
            log.info('skip creation of prodigal training file: genome-size=%i', genome['size'])
    else:
        log.info('use provided prodigal training file: file=%s', prodigal_tf_path)

    sequences = {k['id']: k for k in genome['contigs']}
    cdss = []

    # execute prodigal for non-complete sequences (contigs)
    contigs_path = cfg.tmp_path.joinpath('contigs.fasta')
    contigs = [c for c in genome['contigs'] if not c['complete']]
    if(len(contigs) > 0):
        fasta.export_contigs(contigs, contigs_path)
        log.debug('export contigs: # contigs=%i, path=%s', len(contigs), contigs_path)
        proteins_contigs_path = cfg.tmp_path.joinpath('prodigal.contigs.faa')
        gff_contigs_path = cfg.tmp_path.joinpath('prodigal.contigs.gff')
        log.info('run prodigal: type=contigs, # sequences=%i', len(contigs))
        execute_prodigal(genome, contigs_path, prodigal_tf_path, proteins_path=proteins_contigs_path, gff_path=gff_contigs_path, complete=False)
        cds = parse_prodigal_output(genome, sequences, gff_contigs_path, proteins_contigs_path)
        log.info('contig cds: predicted=%i', len(cds))
        cdss.extend(cds)

    # execute prodigal for complete replicons (chromosomes/plasmids)
    replicons_path = cfg.tmp_path.joinpath('replicons.fasta')
    replicons = [c for c in genome['contigs'] if c['complete']]
    if(len(replicons) > 0):
        fasta.export_contigs(replicons, replicons_path)
        log.debug('export replicons: # sequences=%i, path=%s', len(replicons), replicons_path)
        proteins_replicons_path = cfg.tmp_path.joinpath('prodigal.replicons.faa')
        gff_replicons_path = cfg.tmp_path.joinpath('prodigal.replicons.gff')
        log.info('run prodigal: type=replicons, # sequences=%i', len(replicons))
        execute_prodigal(genome, replicons_path, prodigal_tf_path, proteins_path=proteins_replicons_path, gff_path=gff_replicons_path, complete=True)
        cds = parse_prodigal_output(genome, sequences, gff_replicons_path, proteins_replicons_path)
        log.info('replicon cds: predicted=%i', len(cds))
        cdss.extend(cds)

    log.info('predicted=%i', len(cdss))
    return cdss


def execute_prodigal(genome: dict, contigs_path: Path, traininng_file_path: Path, proteins_path: Path=None, gff_path: Path=None, train: bool=False, complete: bool=False):
    log.debug('execute-prodigal: contigs-path=%s, traininng-file-path=%s, proteins-path=%s, gff-path=%s, train=%s, complete=%s', contigs_path, traininng_file_path, proteins_path, gff_path, train, complete)
    cmd = [
        'prodigal',
        '-i', str(contigs_path),
        '-g', str(cfg.translation_table)  # set translation table
    ]
    if(not complete):
        cmd.append('-c')  # closed ends
    if(train):
        cmd.append('-t')
        cmd.append(str(traininng_file_path))
    else:
        if(genome['size'] < 20000):  # no trainings file provided and not enough sequence information available
            cmd.append('-p')  # run prodigal in meta mode
            cmd.append('meta')
        elif(traininng_file_path):
            cmd.append('-t')
            cmd.append(str(traininng_file_path))
        if(proteins_path):  # aa fasta output
            cmd.append('-a')
            cmd.append(str(proteins_path))
        if(gff_path):  # GFF output
            cmd.append('-f')
            cmd.append('gff')
            cmd.append('-o')
            cmd.append(str(gff_path))

    log.debug('cmd=%s', cmd)
    proc = sp.run(
        cmd,
        cwd=str(cfg.tmp_path),
        env=cfg.env,
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        universal_newlines=True
    )
    if(proc.returncode != 0):
        log.debug('stdout=\'%s\', stderr=\'%s\'', proc.stdout, proc.stderr)
        log.warning('ORFs failed! prodigal-error-code=%d', proc.returncode)
        raise Exception(f'prodigal error! error code: {proc.returncode}')


def parse_prodigal_output(genome: dict, sequences, gff_path: Path, proteins_path: Path):
    log.debug('parse-prodigal-output: gff-path=%s, proteins-path=%s', gff_path, proteins_path)
    cdss = {}
    partial_cdss_per_record = {}
    partial_cdss_per_contig = {k['id']: [] for k in genome['contigs']}
    with gff_path.open() as fh:
        for line in fh:
            if(line[0] != '#'):
                (contig_id, inference, _, start, stop, score, strand, _, annotations_raw) = line.strip().split('\t')
                gff_annotations = split_gff_annotation(annotations_raw)
                contig_orf_id = gff_annotations['ID'].split('_')[1]

                cds = OrderedDict()
                cds['type'] = bc.FEATURE_CDS
                cds['contig'] = contig_id
                cds['start'] = int(start)
                cds['stop'] = int(stop)
                cds['strand'] = bc.STRAND_FORWARD if strand == '+' else bc.STRAND_REVERSE
                cds['gene'] = None
                cds['product'] = None
                cds['start_type'] = gff_annotations['start_type']
                cds['rbs_motif'] = gff_annotations['rbs_motif'] if gff_annotations['rbs_motif'] != 'None' else None
                cds['db_xrefs'] = [so.SO_CDS.id]

                if(cds['strand'] == bc.STRAND_FORWARD):
                    cds['frame'] = (cds['start'] - 1) % 3 + 1
                else:
                    cds['frame'] = (sequences[cds['contig']]['length'] - cds['stop']) % 3 + 1

                if(gff_annotations['partial'] == '10'):
                    cds['truncated'] = bc.FEATURE_END_5_PRIME if cds['strand'] == bc.STRAND_FORWARD else bc.FEATURE_END_3_PRIME
                    partial_cdss_per_record[f"{cds['contig']}_{contig_orf_id}"] = cds
                    partial_cdss_per_contig[cds['contig']].append(cds)
                elif(gff_annotations['partial'] == '01'):
                    cds['truncated'] = bc.FEATURE_END_3_PRIME if cds['strand'] == bc.STRAND_FORWARD else bc.FEATURE_END_5_PRIME
                    partial_cdss_per_record[f"{cds['contig']}_{contig_orf_id}"] = cds
                    partial_cdss_per_contig[cds['contig']].append(cds)
                else:
                    cdss[f"{cds['contig']}_{contig_orf_id}"] = cds

                log.info(
                    'contig=%s, start=%i, stop=%i, strand=%s, frame=%s, truncated=%s, start-type=%s, RBS-motif=%s',
                    cds['contig'], cds['start'], cds['stop'], cds['strand'], cds['frame'], cds.get('truncated', 'no'), cds['start_type'], cds['rbs_motif']
                )

    # extract translated orf sequences
    with proteins_path.open() as fh:
        for record in SeqIO.parse(fh, 'fasta'):
            cds = cdss.get(record.id, None)
            if(cds):
                aa = str(record.seq)[:-1]  # discard trailing asterisk
                cds['aa'] = aa
                cds['aa_digest'], cds['aa_hexdigest'] = bu.calc_aa_hash(aa)
            else:
                partial_cds = partial_cdss_per_record.get(record.id, None)
                if(partial_cds):
                    aa = str(record.seq)
                    if(partial_cds['truncated'] == bc.FEATURE_END_5_PRIME):
                        aa = aa[:-1]  # discard trailing asterisk
                    partial_cds['aa'] = aa
                    partial_cds['aa_digest'], partial_cds['aa_hexdigest'] = bu.calc_aa_hash(aa)
                else:
                    log.warning('unknown sequence detected! id=%s, sequence=%s', record.id, record.seq)
    cdss = list(cdss.values())

    # merge matching partial features on sequence edges
    for partial_cdss in partial_cdss_per_contig.values():
        if(len(partial_cdss) >= 2):  # skip unpaired edge features
            first_partial_cds = partial_cdss[0]  # first partial CDS per contig
            last_partial_cds = partial_cdss[-1]  # last partial CDS per contig
            # check if partial CDSs are on same strand and have opposite truncated edges
            # and firtst starts at 1 and last ends at contig end (length)
            if(first_partial_cds['strand'] == last_partial_cds['strand']
                and first_partial_cds['truncated'] != last_partial_cds['truncated']
                and first_partial_cds['start'] == 1
                and last_partial_cds['stop'] == sequences[last_partial_cds['contig']]['length']
                and sequences[first_partial_cds['contig']]['topology'] == bc.TOPOLOGY_CIRCULAR):
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
                cdss.append(cds)
                log.info(
                    'edge CDS: contig=%s, start=%i, stop=%i, strand=%s, frame=%s, start-type=%s, RBS-motif=%s, aa-hexdigest=%s, aa=[%s..%s]',
                    cds['contig'], cds['start'], cds['stop'], cds['strand'], cds['frame'], cds['start_type'], cds['rbs_motif'], cds['aa_hexdigest'], aa[:10], aa[-10:]
                )
                partial_cdss = partial_cdss[1:-1]
        for partial_cds in partial_cdss:
            cdss.append(partial_cds)
            log.info(
                'truncated CDS: contig=%s, start=%i, stop=%i, strand=%s, frame=%s, truncated=%s, start-type=%s, RBS-motif=%s, aa-hexdigest=%s, aa=[%s..%s]',
                partial_cds['contig'], partial_cds['start'], partial_cds['stop'], partial_cds['strand'], partial_cds['frame'], partial_cds['truncated'], partial_cds['start_type'], partial_cds['rbs_motif'], partial_cds['aa_hexdigest'], partial_cds['aa'][:10], partial_cds['aa'][-10:]
            )

    contigs = {c['id']: c for c in genome['contigs']}
    for cds in cdss:  # extract nt sequences
        contig = contigs[cds['contig']]
        nt = bu.extract_feature_sequence(cds, contig)
        cds['nt'] = nt
        log.info(
            'contig=%s, start=%i, stop=%i, strand=%s, nt=[%s..%s]',
            cds['contig'], cds['start'], cds['stop'], cds['strand'], nt[:10], nt[-10:]
        )
    return cdss


def split_gff_annotation(annotation_string: str) -> Dict[str, str]:
    annotations = {}
    for expr in annotation_string.split(';'):
        if(expr != ''):
            try:
                key, value = expr.split('=')
                annotations[key] = value
            except:
                log.error('expr=%s', expr)
    return annotations


def predict_pfam(cdss: Sequence[dict]) -> Sequence[dict]:
    """Detect Pfam-A entries"""
    fasta_path = cfg.tmp_path.joinpath('hypotheticals.faa')
    orf.write_internal_faa(cdss, fasta_path)
    output_path = cfg.tmp_path.joinpath('cds.hypotheticals.pfam.tsv')
    cmd = [
        'hmmsearch',
        '--noali',
        '--cut_ga',  # use gathering cutoff
        '--domtblout', str(output_path),
        '--cpu', str(cfg.threads if cfg.threads <= 4 else 4),
        str(cfg.db_path.joinpath('pfam')),
        str(fasta_path)
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
    if(proc.returncode != 0):
        log.debug('stdout=\'%s\', stderr=\'%s\'', proc.stdout, proc.stderr)
        log.warning('pfam detection failed! hmmsearch-error-code=%d', proc.returncode)
        raise Exception(f'hmmsearch error! error code: {proc.returncode}')

    pfam_hits = []
    cds_pfams_hits = {}
    orf_by_aa_digest = orf.get_orf_dictionary(cdss)
    with output_path.open() as fh:
        for line in fh:
            if(line[0] != '#'):
                cols = bc.RE_MULTIWHITESPACE.split(line.strip())
                aa_identifier = cols[0]
                cds = orf_by_aa_digest[aa_identifier]

                domain_length = int(cols[5])
                domain_start = int(cols[15])
                domain_stop = int(cols[16])
                domain_cov = (domain_stop - domain_start + 1) / domain_length
                aa_start = int(cols[19])
                aa_stop = int(cols[20])
                aa_cov = (aa_stop - aa_start + 1) / len(cds['aa'])

                pfam = OrderedDict()
                pfam['id'] = cols[4]
                pfam['name'] = cols[3]
                pfam['length'] = domain_length
                pfam['aa_cov'] = aa_cov
                pfam['hmm_cov'] = domain_cov
                pfam['evalue'] = float(cols[6])
                pfam['score'] = float(cols[7])
                pfam['start'] = aa_start
                pfam['stop'] = aa_stop

                if('pfams' not in cds):
                    cds['pfams'] = []
                cds['pfams'].append(pfam)
                if('db_xrefs' not in cds):
                    cds['db_xrefs'] = []
                cds['db_xrefs'].append(f"PFAM:{pfam['id']}")
                pfam_hits.append(cds)
                cds_pfams_hits[aa_identifier] = cds
                log.info(
                    'pfam detected: contig=%s, start=%i, stop=%i, strand=%s, pfam-id=%s, length=%i, aa-start=%i, aa-stop=%i, aa-cov=%1.1f, hmm-cov=%1.1f, evalue=%1.1e, bitscore=%1.1f, name=%s',
                    cds['contig'], cds['start'], cds['stop'], cds['strand'], pfam['id'], pfam['length'], pfam['start'], pfam['stop'], pfam['aa_cov'], pfam['hmm_cov'], pfam['evalue'], pfam['score'], pfam['name']
                )
    log.info('predicted-pfams=%i, CDS-w/-pfams=%i', len(pfam_hits), len(cds_pfams_hits))
    return cds_pfams_hits.values()


def analyze_proteins(cdss: Sequence[dict]):
    for cds in cdss:
        seq = ProteinAnalysis(cds['aa'])
        seq_stats = OrderedDict()
        try:
            seq_stats['molecular_weight'] = seq.molecular_weight()
        except:
            log.warning(
                'could not calc molecular weight! contig=%s, start=%i, stop=%i, strand=%s, frame=%s',
                cds['contig'], cds['start'], cds['stop'], cds['strand'], cds['frame']
            )
            seq_stats['molecular_weight'] = None
        try:
            seq_stats['isoelectric_point'] = seq.isoelectric_point()
        except:
            log.warning(
                'could not calc isoelectric point! contig=%s, start=%i, stop=%i, strand=%s, frame=%s',
                cds['contig'], cds['start'], cds['stop'], cds['strand'], cds['frame']
            )
            seq_stats['isoelectric_point'] = None
        cds['seq_stats'] = seq_stats


def revise_translational_exceptions(genome: dict, cdss: Sequence[dict]):
    """
    Revise translational exceptions as for istance selenocystein proteins.
    """
    no_revised = 0
    if(bc.FEATURE_NC_RNA_REGION not in genome['features']):  # check if ncRNA regions have been detected, otherwise skip analysis and return
        return no_revised

    contigs = {c['id']: c for c in genome['contigs']}
    # detect splitted orphan ORFs of selenocystein proteins that are subject to stop codon recoding.
    cdss_per_contigs = {k['id']: [] for k in genome['contigs']}  # get CDS per contig
    for cds in cdss:
        cdss_per_contig = cdss_per_contigs[cds['contig']]
        if('truncated' not in cds):  # exclude truncated CDS for now
            cdss_per_contig.append(cds)
    cds_pairs_per_contig = {k['id']: [] for k in genome['contigs']}  # extract inframe primate CDS neighbouring pairs
    for id, cdss_per_contig in cdss_per_contigs.items():
        cdss_per_contig = sorted(cdss_per_contig, key=lambda k: k['start'])
        for i in range(1, len(cdss_per_contig)):
            cds_a = cdss_per_contig[i-1]
            cds_b = cdss_per_contig[i]
            strand = cds_a['strand']
            upstream_stop_codon = cds_a['nt'][-3:] if strand == bc.STRAND_FORWARD else cds_b['nt'][-3:]
            if(
                cds_a['strand'] == cds_b['strand'] and  # up- and downstream ORFs on the same strand
                cds_a['frame'] == cds_b['frame'] and  # up- and downstream ORFs on the same frame
                upstream_stop_codon == 'TGA' and  # tRNAScan-SE 2.0 only predicts tRNA-Sec with UCA anticodons, therefore we can only detect TGA stop codons
                (cds_b['start'] - cds_a['stop']) < 100):  # up- and downstream ORFs in close proximity
                cds_pairs = cds_pairs_per_contig[cds_a['contig']]
                cds_pairs.append((cds_a, cds_b))

    recoding_regions = [ncrna_region for ncrna_region in genome['features'][bc.FEATURE_NC_RNA_REGION] if ncrna_region['class'] == so.SO_CIS_REG_RECODING_STIMULATION_REGION]  #  Selenocysteine insertion sequences
    for recoding_region in recoding_regions:
        if('selenocysteine' in recoding_region.get('product', '').lower()):
            cds_pairs = cds_pairs_per_contig[recoding_region['contig']]
            for cds_a, cds_b in cds_pairs:  # find CDS pair around recoding region
                strand = cds_a['strand']
                if(
                    strand == recoding_region['strand'] and  # everything is on the same strand
                    cds_a['start'] < recoding_region['start'] and recoding_region['stop'] < cds_b['stop']):  # recoding region lies between up- and downstream ORFs
                    log.debug(
                        'selenocysteine recoding ncRNA/CDS pair detected: contig=%s, strand=%s, CDS-A=[%i...%i] (%s..%s), recoding-ie=[%i..%i], CDS-B=[%i...%i] (%s..%s)',
                        recoding_region['contig'], recoding_region['strand'], cds_a['start'], cds_a['stop'], cds_a['nt'][:10], cds_a['nt'][-10:], recoding_region['start'], recoding_region['stop'], cds_b['start'], cds_b['stop'], cds_b['nt'][:10], cds_b['nt'][-10:]
                    )
                    seleno_cds = copy.deepcopy(cds_a)
                    seleno_cds['stop'] = cds_b['stop']
                    seleno_cds['rbs_motif'] = cds_a['rbs_motif'] if strand == bc.STRAND_FORWARD else cds_b['rbs_motif']
                    contig = contigs[seleno_cds['contig']]
                    nt = bu.extract_feature_sequence(seleno_cds, contig)
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
                            'selenocysteine CDS detected: contig=%s, start=%i, stop=%i, strand=%s, frame=%i, exception=[%i..%i], nt=[%s..%s], aa=[%s..%s], aa-hexdigest=%s',
                            seleno_cds['contig'], seleno_cds['start'], seleno_cds['stop'], seleno_cds['strand'], seleno_cds['frame'], seleno_cds['exception']['start'], seleno_cds['exception']['stop'], nt[:10], nt[-10:], aa[:10], aa[-10:], seleno_cds['aa_hexdigest']
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
                            'spurious selenocysteine CDS detected: contig=%s, start=%i, stop=%i, strand=%s, frame=%i, nt=[%s], aa=[%s]',
                            seleno_cds['contig'], seleno_cds['start'], seleno_cds['stop'], seleno_cds['strand'], seleno_cds['frame'], nt, aa
                        )
    return no_revised


def revise_special_cases_annotated(genome: dict, cdss: Sequence[dict]):
    """
    Revise rare but known special cases as for istance supposedly truncated dnaA genes on rotated chromosome starts
    which often appear on re-annotated genomes.
    """
    
    contigs = {c['id']: c for c in genome['contigs']}
    # look for supposedly truncated dnaA genes on rotated chromosome starts: start=1, strand=+
    dnaA = None
    for cds in cdss:
        contig = contigs[cds['contig']]
        if(
            contig['complete'] and
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
            'revise supposedly truncated dnaA gene on rotated chromosome start: contig=%s, start=%i, stop=%i, strand=%s, gene=%s, product=%s, nt=[%s..%s], aa=[%s..%s]',
            dnaA['contig'], dnaA['start'], dnaA['stop'], dnaA['strand'], gene, dnaA['product'], dnaA['nt'][:10], dnaA['nt'][-10:], dnaA['aa'][:10], dnaA['aa'][-10:]
        )
    
    # look for supposedly truncated repA genes on rotated plasmid starts: start=1, strand=+
    repAs = []
    for cds in cdss:
        contig = contigs[cds['contig']]
        if(
            contig['complete'] and
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
                'revise supposedly truncated repA gene on rotated plasmid start: contig=%s, start=%i, stop=%i, strand=%s, gene=%s, product=%s, nt=[%s..%s], aa=[%s..%s]',
                repA['contig'], repA['start'], repA['stop'], repA['strand'], gene, repA['product'], repA['nt'][:10], repA['nt'][-10:], repA['aa'][:10], repA['aa'][-10:]
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
        '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'qstart', 'qend', 'sstart', 'send', 'full_sseq',
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
            (aa_identifier, cluster_id, identity, alignment_length, query_start, query_end, subject_start, subject_end, subject_sequence) = line.rstrip('\n').split('\t')
            cds = cds_by_hexdigest[aa_identifier]
            query_cov = int(alignment_length) / len(cds['aa'])
            subject_cov = int(alignment_length) / len(subject_sequence)
            identity = float(identity) / 100
            if query_cov >= bc.MIN_PSEUDOGENE_QUERY_COVERAGE and identity >= bc.MIN_PSEUDOGENE_IDENTITY \
                    and bc.MIN_PSEUDOGENE_SUBJECT_COVERAGE <= subject_cov < bc.MIN_PSC_COVERAGE:
                cds['pseudo-inference'] = {
                    DB_PSC_COL_UNIREF90: cluster_id,
                    'query-cov': query_cov,
                    'subject-cov': subject_cov,
                    'identity': identity,
                    'gene_start': int(query_start),
                    'gene_end': int(query_end),
                    'reference_start': int(subject_start),
                    'reference_end': int(subject_end),
                    'reference_sequence': subject_sequence
                }
                pseudo_candidates.append(cds)
                log.debug(
                    'pseudogene-candidate: contig=%s, start=%i, stop=%i, strand=%s, aa-length=%i, query-cov=%0.3f, subject-cov=%0.3f, identity=%0.3f, UniRef90=%s',
                    cds['contig'], cds['start'], cds['stop'], cds['strand'], len(cds['aa']), query_cov, subject_cov, identity, cluster_id
                )
    log.info('found: pseudogene-candidates=%i', len(pseudo_candidates))
    return pseudo_candidates


def detect_pseudogenes(candidates: Sequence[dict], cdss: Sequence[dict], genome: dict) -> Sequence[dict]:
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
    contigs = {c['id']: c for c in genome['contigs']}
    candidates_extended_positions = {}
    with candidates_elongated_sequences_path.open(mode='w') as fh:
        for cds in candidates:
            contig = contigs[cds['contig']]
            cds_elongated = get_elongated_cds(cds, contig)
            seq = bu.extract_feature_sequence(cds_elongated, contig)
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

                    if alignment_length == len(cds['aa']):  # skip non-extended genes (full match)
                        log.debug(
                            'no pseudogene (full match): contig=%s, start=%i, stop=%i, strand=%s',
                            cds['contig'], cds['start'], cds['stop'], cds['strand']
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
                            'paralog': is_paralog(uniref90_by_hexdigest, aa_identifier, cluster_id)
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
                            truncation = bc.FEATURE_END_BOTH
                        elif bc.FEATURE_END_5_PRIME in directions:
                            truncation = bc.FEATURE_END_5_PRIME if cds['strand'] == bc.STRAND_FORWARD else bc.FEATURE_END_3_PRIME
                        elif bc.FEATURE_END_3_PRIME in directions:
                            truncation = bc.FEATURE_END_3_PRIME if cds['strand'] == bc.STRAND_FORWARD else bc.FEATURE_END_5_PRIME
                        cds['truncated'] = truncation

                        cds['pseudo'] = True
                        cds[bc.PSEUDOGENE] = pseudogene
                        cds.pop('hypothetical')
                        pseudogenes.append(cds)
                        log.info(
                            'pseudogene: contig=%s, start=%i, stop=%i, strand=%s, insertions=%s, deletions=%s, mutations=%s, effect=%s',
                            cds['contig'], cds['start'], cds['stop'], cds['strand'], observations.get(bc.PSEUDOGENE_CAUSE_INSERTION, []), observations.get(bc.PSEUDOGENE_CAUSE_DELETION, []), observations.get(bc.PSEUDOGENE_CAUSE_MUTATION, []), effects
                        )

                    elif observations[bc.PSEUDOGENE_EXCEPTION_SELENOCYSTEINE] or observations[bc.PSEUDOGENE_EXCEPTION_PYROLYSINE]:
                        # TODO handle translation exceptions, correct annotation
                        pass

    for cds in candidates:
        cds.pop('pseudo-inference')
    log.info('found: pseudogenes=%i', len(pseudogenes))
    return pseudogenes


def get_elongated_cds(cds: dict, contig: dict, offset: int = bc.PSEUDOGENE_OFFSET) -> Dict[str, Union[int, str, bool]]:
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

    contig_length = len(contig['sequence'])
    if contig['topology'] == 'circular' and elongated_cds['start'] - offset < 0:
        elongated_cds['start'] = contig_length + elongated_cds['start'] - offset
        elongated_cds['edge'] = True
    elif elongated_cds['start'] - offset < 0:
        elongated_cds['start'] = 1
        elongated_cds['elongation_upstream'] = cds['start']
    else:
        elongated_cds['start'] = elongated_cds['start'] - offset

    if contig['topology'] == 'circular' and elongated_cds['stop'] + offset > contig_length:
        elongated_cds['stop'] = elongated_cds['stop'] + offset - contig_length
        elongated_cds['edge'] = True
    elif elongated_cds['stop'] + offset > contig_length:
        elongated_cds['stop'] = contig_length
        elongated_cds['elongation_downstream'] = contig_length - cds['stop']
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
                    'pseudogene observation: contig=%s, start=%i, stop=%i, strand=%s, original start=%i',
                    cds['contig'], cds['start'], cds['stop'], cds['strand'], cds['start'] + genome_position
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
                'pseudogene observation: contig=%s, start=%i, stop=%i, strand=%s, cause=insertion, position=%i',
                cds['contig'], cds['start'], cds['stop'], cds['strand'], genome_position
            )
            alignment_position += 1
        elif char == '/':  # deletion
            deletions += 1
            genome_position = get_abs_position(cds, start, alignment_position, edge)
            observations[bc.PSEUDOGENE_CAUSE_DELETION].add(genome_position)
            observations['directions'].add(get_direction(alignment_position, edge))
            log.info(
                'pseudogene observation: contig=%s, start=%i, stop=%i, strand=%s, cause=deletion, position=%i',
                cds['contig'], cds['start'], cds['stop'], cds['strand'], genome_position
            )
        elif char == '*':  # stop codon, selenocysteine, pyrolysine
            if ref_char == 'U':  # selenocysteine
                genome_position = get_abs_position(cds, start, alignment_position, edge)
                observations[bc.PSEUDOGENE_EXCEPTION_SELENOCYSTEINE].add(genome_position)
                log.info(
                    'pseudogene observation: contig=%s, start=%i, stop=%i, strand=%s, exception=selenocysteine, position=%i',
                    cds['contig'], cds['start'], cds['stop'], cds['strand'], genome_position
                )
            elif ref_char == 'O':  # pyrolysine
                genome_position = get_abs_position(cds, start, alignment_position, edge)
                observations[bc.PSEUDOGENE_EXCEPTION_PYROLYSINE].add(genome_position)
                log.info(
                    'pseudogene observation: contig=%s, start=%i, stop=%i, strand=%s, exception=pyrolysin, position=%i',
                    cds['contig'], cds['start'], cds['stop'], cds['strand'], genome_position
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
                    'pseudogene observation: contig=%s, start=%i, stop=%i, strand=%s, effect=stop%s, position=%i',
                    cds['contig'], cds['start'], cds['stop'], cds['strand'], mutation, genome_position
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
