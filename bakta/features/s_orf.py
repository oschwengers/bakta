import concurrent.futures as cf
import logging
import math
import subprocess as sp

from collections import OrderedDict
from typing import Sequence

from Bio.Seq import Seq

import bakta.config as cfg
import bakta.constants as bc
import bakta.features.orf as orf
import bakta.utils as bu
import bakta.psc as psc
import bakta.pscc as pscc
import bakta.so as so


log = logging.getLogger('S_ORF')


def extract(data: dict):
    """Predict open reading frames in mem via BioPython."""
    orfs = []
    for seq in data['sequences']:
        nt_seq = Seq(seq['nt'])
        for strand, strand_nt_seq in [(bc.STRAND_FORWARD, nt_seq), (bc.STRAND_REVERSE, nt_seq.reverse_complement())]:  # strands +/-
            for frame in range(3):  # frames 1/2/3 -> 0, 1, 2
                seq_frame = strand_nt_seq[frame:]

                # remove non-triplet tail nucleotides
                residue = len(seq_frame) % 3
                if(residue != 0):
                    seq_frame = seq_frame[:-residue]

                aa_seq = str(seq_frame.translate(table=cfg.translation_table, stop_symbol='*', to_stop=False, cds=False))
                aa_start = aa_seq.find('M')
                aa_end = aa_seq.find('*', aa_start)
                while aa_start > -1 and aa_end > -1:
                    orf_length = aa_end - aa_start
                    if(orf_length >= bc.MIN_SORF_LENGTH and orf_length < bc.MAX_SORF_LENGTH):  # get all CDS starts (M)
                        aa = aa_seq[aa_start:aa_end]
                        (aa_digest, aa_hexdigest) = bu.calc_aa_hash(aa)
                        if(strand == bc.STRAND_FORWARD):
                            dna_start = aa_start * 3 + frame + 1  # +1: 0 based idx to 1 based
                            dna_stop = aa_end * 3 + 2 + frame + 1
                        else:
                            dna_start = len(strand_nt_seq) - frame - (aa_end + 1) * 3 + 1
                            dna_stop = len(strand_nt_seq) - frame - aa_start * 3
                        
                        sorf = OrderedDict()
                        sorf['type'] = bc.FEATURE_SORF
                        sorf['sequence'] = seq['id']
                        sorf['start'] = dna_start
                        sorf['stop'] = dna_stop
                        sorf['strand'] = strand
                        sorf['gene'] = None
                        sorf['product'] = None
                        sorf['frame'] = frame + 1
                        sorf['db_xrefs'] = [so.SO_SORF.id]
                        sorf['aa'] = aa
                        sorf['aa_digest'] = aa_digest
                        sorf['aa_hexdigest'] = aa_hexdigest
                        
                        nt = bu.extract_feature_sequence(sorf, seq)  # extract nt sequences
                        sorf['nt'] = nt

                        orfs.append(sorf)
                        log.debug(
                            'seq=%s, start=%i, stop=%i, strand=%s, frame=%i, aa-length=%i, aa=%s, nt=[%s..%s]',
                            seq['id'], sorf['start'], sorf['stop'], strand, frame, len(aa), aa, nt[:10], nt[-10:]
                        )
                    aa_start = aa_seq.find('M', aa_start + 1)
                    if(aa_start > aa_end):
                        aa_end = aa_seq.find('*', aa_start)

    log.info('predicted=%i', len(orfs))
    return orfs


def get_feature_start(feature: dict) -> int:
    return feature['start'] if feature['strand'] == bc.STRAND_FORWARD else feature['stop']


def get_feature_stop(feature: dict) -> int:
    return feature['stop'] if feature['strand'] == bc.STRAND_FORWARD else feature['start']


def overlap_filter(data: dict, orfs_raw: Sequence[dict]):
    """Filter in-mem ORFs by overlapping CDSs."""
    t_rnas_per_sequence = {seq['id']: [] for seq in data['sequences']}
    for t_rna in [feat for feat in data['features'] if feat['type'] == bc.FEATURE_T_RNA]:
        t_rnas = t_rnas_per_sequence[t_rna['sequence']]
        t_rnas.append(t_rna)
    for tm_rna in [feat for feat in data['features'] if feat['type'] == bc.FEATURE_TM_RNA]:
        t_rnas = t_rnas_per_sequence[tm_rna['sequence']]
        t_rnas.append(tm_rna)

    r_rna_per_sequence = {seq['id']: [] for seq in data['sequences']}
    for r_rna in [feat for feat in data['features'] if feat['type'] == bc.FEATURE_R_RNA]:
        r_rnas = r_rna_per_sequence[r_rna['sequence']]
        r_rnas.append(r_rna)

    # nc_rnas_per_sequence = {seq['id']: [] for seq in data['sequences']}
    # for nc_rna in [feat for feat in data['features'] if feat['type'] == bc.FEATURE_NC_RNA]:
    #     nc_rnas = nc_rnas_per_sequence[nc_rna['sequence']]
    #     nc_rnas.append(nc_rna)
    # for nc_rna in [feat for feat in data['features'] if feat['type'] == bc.FEATURE_NC_RNA_REGION]:
    #     nc_rnas = nc_rnas_per_sequence[nc_rna['sequence']]
    #     nc_rnas.append(nc_rna)

    crispr_arrays_per_sequence = {seq['id']: [] for seq in data['sequences']}
    for crispr_array in [feat for feat in data['features'] if feat['type'] == bc.FEATURE_CRISPR]:
        crispr_arrays = crispr_arrays_per_sequence[crispr_array['sequence']]
        crispr_arrays.append(crispr_array)

    cdss_per_sequence = {seq['id']: [] for seq in data['sequences']}
    for cds in [feat for feat in data['features'] if feat['type'] == bc.FEATURE_CDS]:
        cdss = cdss_per_sequence[cds['sequence']]
        cdss.append(cds)

    sorfs_per_sequence = {seq['id']: [] for seq in data['sequences']}
    for sorf in orfs_raw:
        orfs = sorfs_per_sequence[sorf['sequence']]
        orfs.append(sorf)

    discarded_sorf_keys = set()
    with cf.ProcessPoolExecutor(max_workers=cfg.threads) as tpe:
        futures = []
        for seq in data['sequences']:
            sequence_sorfs = sorfs_per_sequence[seq['id']]
            log.debug('filter: seq=%s, # sORFs=%i', seq['id'], len(sequence_sorfs))
            if(len(sequence_sorfs) < 100):  # execute sORF filter task
                sorf_keys = filter_sorf(sequence_sorfs, cdss_per_sequence[seq['id']], r_rna_per_sequence[seq['id']], t_rnas_per_sequence[seq['id']], crispr_arrays_per_sequence[seq['id']])
                for sorf_key in [sk for sk in sorf_keys if sk is not None]:
                    discarded_sorf_keys.add(sorf_key)
            elif(len(sequence_sorfs) < 1000):  # submit sORF filter task to thread pool
                futures.append(tpe.submit(filter_sorf, sequence_sorfs, cdss_per_sequence[seq['id']], r_rna_per_sequence[seq['id']], t_rnas_per_sequence[seq['id']], crispr_arrays_per_sequence[seq['id']]))
            else:  # submit sORF chunk filter tasks to thread pool
                chunk_size = math.ceil(len(sequence_sorfs) / cfg.threads) if (len(sequence_sorfs) >= cfg.threads * 1000) else 1000
                log.debug('filter: chunk-size=%i', chunk_size)
                for i in range(0, len(sequence_sorfs), chunk_size):
                    sorf_chunk = sequence_sorfs[i:i + chunk_size]
                    log.debug('filter chunk: i=%i, chunk-elements=%i', i, len(sorf_chunk))
                    futures.append(tpe.submit(filter_sorf, sorf_chunk, cdss_per_sequence[seq['id']], r_rna_per_sequence[seq['id']], t_rnas_per_sequence[seq['id']], crispr_arrays_per_sequence[seq['id']]))
        for f in futures:
            for sorf_key in [sk for sk in f.result() if sk is not None]:
                discarded_sorf_keys.add(sorf_key)

    valid_sorfs = []
    discarded_sorfs = []
    for sorfs in sorfs_per_sequence.values():
        for sorf in sorfs:
            key = orf.get_orf_key(sorf)
            if(key in discarded_sorf_keys):
                discarded_sorfs.append(sorf)
            else:
                valid_sorfs.append(sorf)

    log.info('sORF filter: # valid=%i, # discarded=%i', len(valid_sorfs), len(discarded_sorfs))
    return valid_sorfs, discarded_sorfs


def filter_sorf(sorf_chunk: Sequence[dict], sequence_cdss: Sequence[dict], sequence_r_rnas: Sequence[dict], sequence_t_rnas: Sequence[dict], sequence_crispr_arrays: Sequence[dict]):
    discarded_sorf_keys = []
    for sorf in sorf_chunk:
        break_flag = False
        # filter CDS overlapping ORFs
        for cds in sequence_cdss:
            # log.debug('filter short ORFs by CDS: %s%i[%i->%i]', cds['strand'], cds['frame'], cds['start'], cds['stop'])
            if(sorf['strand'] == cds['strand']):
                if(sorf['frame'] == cds['frame']):
                    # fast/simple overlap detection for in frame orientation
                    if(cds['stop'] < sorf['start'] or cds['start'] > sorf['stop']):
                        continue
                    else:
                        discarded_sorf_keys.append(orf.get_orf_key(sorf))
                        break_flag = True
                        break
                else:
                    # if(sorf['start'] < cds['start'] and sorf['stop'] > cds['start']):  # ToDo: add max overlap threshold
                    #     # out-of-frame sorf partially overlapping CDS upstream
                    #     continue
                    # elif(sorf['start'] >= cds['start'] and sorf['stop'] <= cds['stop']):
                    if(sorf['start'] >= cds['start'] and sorf['stop'] <= cds['stop']):
                        # out-of-frame sorf completely overlapped by CDS
                        discarded_sorf_keys.append(orf.get_orf_key(sorf))
                        break_flag = True
                        break
                    # elif(sorf['start'] < cds['stop'] and sorf['stop'] > cds['stop']):  # ToDo: add max overlap threshold
                    #     # out-of-frame sorf partially overlapping CDS downstream
                    #     continue
            else:
                if(sorf['frame'] == cds['frame']):
                    # fast/simple overlap detection for in frame orientation
                    if(cds['stop'] < sorf['start'] or cds['start'] > sorf['stop']):
                        continue
                    else:
                        discarded_sorf_keys.append(orf.get_orf_key(sorf))
                        break_flag = True
                        break
                else:
                    # ToDo: add out-of-frame filters
                    # if(sorf['start'] < cds['start'] and sorf['stop'] > cds['start']):  # ToDo: add max overlap threshold
                    #     # out-frame sorf partially overlapping CDS upstream
                    #     continue
                    # elif(sorf['start'] >= cds['start'] and sorf['stop'] <= cds['stop']):
                    if(sorf['start'] >= cds['start'] and sorf['stop'] <= cds['stop']):
                        # out-frame sorf completely overlapped by CDS
                        discarded_sorf_keys.append(orf.get_orf_key(sorf))
                        break_flag = True
                        break
                    # elif(sorf['start'] < cds['stop'] and sorf['stop'] > cds['stop']):  # ToDo: add max overlap threshold
                    #     # out-frame sorf partially overlapping CDS downstream
                    #     continue
        if(break_flag):
            continue

        # filter rRNA overlapping ORFs
        for r_rna in sequence_r_rnas:
            # log.debug('filter short ORFs by rRNA: %s[%i->%i]', r_rna['strand'], r_rna['start'], r_rna['stop'])
            # fast/simple overlap detection for rRNAs
            if(sorf['stop'] < r_rna['start'] or sorf['start'] > r_rna['stop']):
                continue
            else:
                discarded_sorf_keys.append(orf.get_orf_key(sorf))
                break_flag = True
                break
        if(break_flag):
            continue

        # filter tRNA overlapping ORFs
        # log.debug('filter short ORFs by tRNA: %s[%i->%i]', t_rna['strand'], t_rna['start'], t_rna['stop'])
        for t_rna in sequence_t_rnas:
            # fast/simple overlap detection for tRNAs
            if(sorf['stop'] < t_rna['start'] or sorf['start'] > t_rna['stop']):
                continue
            else:
                discarded_sorf_keys.append(orf.get_orf_key(sorf))
                break_flag = True
                break
        if(break_flag):
            continue

        # filter CRISPR array overlapping ORFs
        # log.debug('filter short ORFs by CRISPR: [%i->%i]', crispr['start'], crispr['stop'])
        for crispr in sequence_crispr_arrays:
            # fast/simple overlap detection for CRISPR
            if(sorf['stop'] < crispr['start'] or sorf['start'] > crispr['stop']):
                continue
            else:
                discarded_sorf_keys.append(orf.get_orf_key(sorf))
                break_flag = True
                break
        if(break_flag):
            continue

    return discarded_sorf_keys


def annotation_filter(sorfs: Sequence[dict]) -> Sequence[dict]:
    """Filter sORFs according to available annotations."""
    valid_sorfs = []
    for sorf in sorfs:
        gene = None
        product = None

        ips = sorf.get('ips', None)
        if(ips is not None):
            tmp = ips.get('gene', '')
            if(tmp != ''):
                gene = tmp
            tmp = ips.get('product', '')
            if(tmp != ''):
                product = tmp

        psc = sorf.get('psc', None)
        if(psc is not None):
            tmp = psc.get('gene', '')
            if(tmp != ''):
                gene = tmp
            tmp = psc.get('product', '')
            if(tmp != ''):
                product = tmp

        pscc = sorf.get('pscc', None)
        if(pscc is not None):
            tmp = pscc.get('product', '')
            if(tmp != ''):
                product = tmp

        if(gene is None and product is None):
            sorf['hypothetical'] = True
        else:
            valid_sorfs.append(sorf)

    return valid_sorfs


def search_pscs(sorfs: Sequence[dict]):
    """Conduct homology search of sORFs against the full db PSCs"""
    return search(sorfs, 'full')


def search_psccs(sorfs: Sequence[dict]):
    """Conduct homology search of sORFs against the light db PSCCs"""
    return search(sorfs, 'light')


def search(sorfs: Sequence[dict], cluster_type: str):
    """Conduct homology search of sORFs against sORF db."""
    sorf_aa_path = cfg.tmp_path.joinpath('sorf.faa')
    orf.write_internal_faa(sorfs, sorf_aa_path)
    diamond_output_path = cfg.tmp_path.joinpath('diamond.sorf.tsv')
    diamond_db_path = cfg.db_path.joinpath('sorf.dmnd')
    cmd = [
        'diamond',
        'blastp',
        '--db', str(diamond_db_path),
        '--query', str(sorf_aa_path),
        '--out', str(diamond_output_path),
        '--id', str(int(bc.MIN_SORF_IDENTITY * 100)),  # '90',
        '--query-cover', str(int(bc.MIN_SORF_COVERAGE * 100)),  # '90'
        '--subject-cover', str(int(bc.MIN_SORF_COVERAGE * 100)),  # '90'
        '--max-target-seqs', '1',  # single best output
        '--outfmt', '6', 'qseqid', 'sseqid', 'qlen', 'slen', 'length', 'pident', 'evalue', 'bitscore',
        '--threads', str(cfg.threads),
        '--tmpdir', str(cfg.tmp_path),  # use tmp folder
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
    if(proc.returncode != 0):
        log.debug('stdout=\'%s\', stderr=\'%s\'', proc.stdout, proc.stderr)
        log.warning('sORF failed! diamond-error-code=%d', proc.returncode)
        raise Exception(f'diamond error! error code: {proc.returncode}')

    sorf_by_aa_digest = orf.get_orf_dictionary(sorfs)
    with diamond_output_path.open() as fh:
        for line in fh:
            (sorf_hash, cluster_id, query_length, subject_length, alignment_length, identity, evalue, bitscore) = line.split('\t')
            sorf = sorf_by_aa_digest[sorf_hash]
            query_cov = int(alignment_length) / len(sorf['aa'])
            subject_cov = int(alignment_length) / int(subject_length)
            identity = float(identity) / 100
            bitscore = float(bitscore)
            evalue = float(evalue)
            if(query_cov >= bc.MIN_SORF_COVERAGE and subject_cov >= bc.MIN_SORF_COVERAGE and identity >= bc.MIN_SORF_IDENTITY):
                result = {
                    'query_cov': query_cov,
                    'subject_cov': subject_cov,
                    'identity': identity,
                    'score': bitscore,
                    'evalue': evalue
                }
                result[psc.DB_PSC_COL_UNIREF90 if cluster_type == 'full' else pscc.DB_PSCC_COL_UNIREF50] = cluster_id
                sorf['psc' if cluster_type == 'full' else 'pscc'] = result
                log.info(
                    'homology: seq=%s, start=%i, stop=%i, strand=%s, aa-length=%i, query-cov=%0.3f, subject-cov=%0.3f, identity=%0.3f, score=%0.1f, evalue=%1.1e, UniRef=%s',
                    sorf['sequence'], sorf['start'], sorf['stop'], sorf['strand'], len(sorf['aa']), query_cov, subject_cov, identity, bitscore, evalue, cluster_id
                )

    sorfs_found = []
    sorfs_not_found = []
    for sorf in sorfs:
        if('psc' in sorf  or  'pscc' in sorf):
            sorfs_found.append(sorf)
        else:
            sorfs_not_found.append(sorf)
    log.info('found=%i', len(sorfs_found))
    return sorfs_found, sorfs_not_found
