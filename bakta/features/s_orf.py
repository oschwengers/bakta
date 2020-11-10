
import logging
import subprocess as sp
from collections import OrderedDict

from Bio.Seq import Seq

import bakta.config as cfg
import bakta.constants as bc
import bakta.utils as bu
import bakta.psc as psc

log = logging.getLogger('S_ORF')


def extract(genome):
    """Predict open reading frames in mem via BioPython."""

    orfs = []
    for contig in genome['contigs']:
        dna_seq = Seq(contig['sequence'])
        for strand, seq in [(bc.STRAND_FORWARD, dna_seq), (bc.STRAND_REVERSE, dna_seq.reverse_complement())]:  # strands +/-
            for frame in range(3):  # frames 1/2/3 -> 0, 1, 2
                seq_frame = seq[frame:]

                # remove non-triplet tail nucleotides
                residue = len(seq_frame) % 3
                if(residue != 0):
                    seq_frame = seq_frame[:-residue]

                aa_seq = str(seq_frame.translate(table=11, stop_symbol='*', to_stop=False, cds=False))
                aa_start = aa_seq.find('M')
                aa_end = aa_seq.find('*', aa_start)
                while aa_start > -1 and aa_end > -1:
                    orf_length = aa_end - aa_start
                    if(orf_length >= bc.MIN_SORF_LENGTH and orf_length <= bc.MAX_SORF_LENGTH):  # get all CDS starts (M)
                        if(strand == '+'):
                            dna_start = aa_start * 3 + frame + 1  # +1: 0 based idx to 1 based
                            dna_stop = aa_end * 3 + 2 + frame + 1
                        else:
                            dna_start = len(seq) - frame - (aa_end + 1) * 3 + 1
                            dna_stop = len(seq) - frame - aa_start * 3
                        sequence = aa_seq[aa_start:aa_end]

                        test_dna_seq = Seq(contig['sequence'][dna_start - 1:dna_stop])
                        if(strand == '-'):
                            test_dna_seq = test_dna_seq.reverse_complement()
                        test_seq = test_dna_seq.translate(table=11, stop_symbol='*', to_stop=False, cds=True)
                        assert sequence == test_seq, f'seqs not equal! a={sequence}, b={test_seq}'
                        (aa_digest, aa_hexdigest) = bu.calc_aa_hash(sequence)

                        sorf = OrderedDict()
                        sorf['type'] = bc.FEATURE_SORF
                        sorf['contig'] = contig['id']
                        sorf['start'] = dna_start
                        sorf['stop'] = dna_stop
                        sorf['strand'] = strand
                        sorf['product'] = None
                        sorf['gene'] = None
                        sorf['frame'] = frame + 1
                        sorf['db_xrefs'] = []
                        sorf['sequence'] = sequence
                        sorf['aa_digest'] = aa_digest
                        sorf['aa_hexdigest'] = aa_hexdigest
                        
                        orfs.append(sorf)
                        log.debug(
                            'contig=%s, start=%i, stop=%i, strand=%s, frame=%i, length=%i, seq=%s',
                            contig['id'], sorf['start'], sorf['stop'], strand, frame, len(sequence), sequence
                        )
                    aa_start = aa_seq.find('M', aa_start + 1)
                    if(aa_start > aa_end):
                        aa_end = aa_seq.find('*', aa_start)

    log.info('predicted=%i', len(orfs))
    return orfs


def get_feature_start(feature):
    return feature['start'] if feature['strand'] == '+' else feature['stop']


def get_feature_stop(feature):
    return feature['stop'] if feature['strand'] == '+' else feature['start']


def overlap_filter(genome, orfs_raw):
    """Filter in-mem ORFs by overlapping CDSs."""

    contig_t_rnas = {k['id']: [] for k in genome['contigs']}
    for t_rna in genome['features'].get(bc.FEATURE_T_RNA, []):
        t_rnas = contig_t_rnas[t_rna['contig']]
        t_rnas.append(t_rna)
    for tm_rna in genome['features'].get(bc.FEATURE_TM_RNA, []):
        t_rnas = contig_t_rnas[tm_rna['contig']]
        t_rnas.append(tm_rna)

    contig_r_rnas = {k['id']: [] for k in genome['contigs']}
    for r_rna in genome['features'].get(bc.FEATURE_R_RNA, []):
        r_rnas = contig_r_rnas[r_rna['contig']]
        r_rnas.append(r_rna)

    # contig_nc_rnas = {k['id']: [] for k in genome['contigs']}
    # for nc_rna in genome['features'].get(bc.FEATURE_NC_RNA, []):
    #     nc_rnas = contig_nc_rnas[nc_rna['contig']]
    #     nc_rnas.append(nc_rna)
    # for nc_rna in genome['features'].get(bc.FEATURE_NC_RNA_REGION, []):
    #     nc_rnas = contig_nc_rnas[nc_rna['contig']]
    #     nc_rnas.append(nc_rna)

    contig_crispr_arrays = {k['id']: [] for k in genome['contigs']}
    for crispr_array in genome['features'].get(bc.FEATURE_CRISPR, []):
        crispr_arrays = contig_crispr_arrays[crispr_array['contig']]
        crispr_arrays.append(crispr_array)

    contig_cdss = {k['id']: [] for k in genome['contigs']}
    for cds in genome['features'].get(bc.FEATURE_CDS, []):
        cdss = contig_cdss[cds['contig']]
        cdss.append(cds)

    contig_sorfs = {k['id']: [] for k in genome['contigs']}
    for sorf in orfs_raw:
        orfs = contig_sorfs[sorf['contig']]
        orfs.append(sorf)

    discarded_sorfs = {}
    for contig in genome['contigs']:
        log.debug('filter short ORFs on contig: %s', contig['id'])
        sorfs = contig_sorfs[contig['id']]

        # filter CDS overlapping ORFs
        for sorf in sorfs:
            for cds in contig_cdss[contig['id']]:
                # log.debug('filter short ORFs by CDS: %s%i[%i->%i]', cds['strand'], cds['frame'], cds['start'], cds['stop'])
                if(sorf['strand'] == cds['strand']):
                    if(sorf['frame'] == cds['frame']):
                        # fast/simple overlap detection for in frame orientation
                        if(cds['stop'] < sorf['start'] or cds['start'] > sorf['stop']):
                            continue
                        else:
                            key = f"{sorf['contig']}-{sorf['start']}-{sorf['stop']}-{sorf['strand']}-{sorf['aa_hexdigest']}"
                            discarded_sorfs[key] = sorf
                    else:
                        if(sorf['start'] < cds['start'] and sorf['stop'] > cds['start']):
                            # out-of-frame sorf partially overlapping CDS upstream
                            # ToDo: add max overlap threshold
                            continue
                        elif(sorf['start'] >= cds['start'] and sorf['stop'] <= cds['stop']):
                            # out-of-frame sorf completely overlapped by CDS
                            key = f"{sorf['contig']}-{sorf['start']}-{sorf['stop']}-{sorf['strand']}-{sorf['aa_hexdigest']}"
                            discarded_sorfs[key] = sorf
                        elif(sorf['start'] < cds['stop'] and sorf['stop'] > cds['stop']):
                            # out-of-frame sorf partially overlapping CDS downstream
                            # ToDo: add max overlap threshold
                            continue
                else:
                    if(sorf['frame'] == cds['frame']):
                        # fast/simple overlap detection for in frame orientation
                        if(cds['stop'] < sorf['start'] or cds['start'] > sorf['stop']):
                            continue
                        else:
                            key = f"{sorf['contig']}-{sorf['start']}-{sorf['stop']}-{sorf['strand']}-{sorf['aa_hexdigest']}"
                            discarded_sorfs[key] = sorf
                    else:
                        # ToDo: add out-of-frame filters
                        if(sorf['start'] < cds['start'] and sorf['stop'] > cds['start']):
                            # out-frame sorf partially overlapping CDS upstream
                            # ToDo: add max overlap threshold
                            continue
                        elif(sorf['start'] >= cds['start'] and sorf['stop'] <= cds['stop']):
                            # out-frame sorf completely overlapped by CDS
                            key = f"{sorf['contig']}-{sorf['start']}-{sorf['stop']}-{sorf['strand']}-{sorf['aa_hexdigest']}"
                            discarded_sorfs[key] = sorf
                        elif(sorf['start'] < cds['stop'] and sorf['stop'] > cds['stop']):
                            # out-frame sorf partially overlapping CDS downstream
                            # ToDo: add max overlap threshold
                            continue

            # filter rRNA overlapping ORFs
            for r_rna in contig_r_rnas[contig['id']]:
                # log.debug('filter short ORFs by rRNA: %s[%i->%i]', r_rna['strand'], r_rna['start'], r_rna['stop'])
                # fast/simple overlap detection for rRNAs
                if(sorf['stop'] < r_rna['start'] or sorf['start'] > r_rna['stop']):
                    continue
                else:
                    key = f"{sorf['contig']}-{sorf['start']}-{sorf['stop']}-{sorf['strand']}-{sorf['aa_hexdigest']}"
                    discarded_sorfs[key] = sorf

            # filter tRNA overlapping ORFs
            # log.debug('filter short ORFs by tRNA: %s[%i->%i]', t_rna['strand'], t_rna['start'], t_rna['stop'])
            for t_rna in contig_t_rnas[contig['id']]:
                # fast/simple overlap detection for tRNAs
                if(sorf['stop'] < t_rna['start'] or sorf['start'] > t_rna['stop']):
                    continue
                else:
                    key = f"{sorf['contig']}-{sorf['start']}-{sorf['stop']}-{sorf['strand']}-{sorf['aa_hexdigest']}"
                    discarded_sorfs[key] = sorf
        
            # filter CRISPR array overlapping ORFs
            # log.debug('filter short ORFs by CRISPR: [%i->%i]', crispr['start'], crispr['stop'])
            for crispr in contig_crispr_arrays[contig['id']]:
                # fast/simple overlap detection for CRISPR
                if(sorf['stop'] < crispr['start'] or sorf['start'] > crispr['stop']):
                    continue
                else:
                    key = f"{sorf['contig']}-{sorf['start']}-{sorf['stop']}-{sorf['strand']}-{sorf['aa_hexdigest']}"
                    discarded_sorfs[key] = sorf
        

    valid_sorfs = []
    for sorfs in contig_sorfs.values():
        for sorf in sorfs:
            key = f"{sorf['contig']}-{sorf['start']}-{sorf['stop']}-{sorf['strand']}-{sorf['aa_hexdigest']}"
            if(key not in discarded_sorfs):
                valid_sorfs.append(sorf)

    log.info('sORF filter: # valid=%i, # discarded=%i', len(valid_sorfs), len(discarded_sorfs.values()))
    return valid_sorfs, discarded_sorfs


def annotation_filter(sorfs):
    """Filter sORFs according to available annotations."""
    valid_sorfs = []
    for sorf in sorfs:
        gene = None
        product = None
        added = False
        
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

        if(gene is None and product is None):
            sorf['hypothetical'] = True
        else:
            valid_sorfs.append(sorf)
    
    return valid_sorfs
    

def search_pscs(sorfs):
    """Conduct homology search of sORFs against sORF db."""
    sorf_fasta_path = cfg.tmp_path.joinpath('sorf.faa')
    with sorf_fasta_path.open(mode='w') as fh:
        for sorf in sorfs:
            fh.write(f">{sorf['aa_hexdigest']}\n{sorf['sequence']}\n")
    diamond_output_path = cfg.tmp_path.joinpath('diamond.sorf.tsv')
    diamond_db_path = cfg.db_path.joinpath('sorf.dmnd')
    cmd = [
        'diamond',
        'blastp',
        '--db', str(diamond_db_path),
        '--query', str(sorf_fasta_path),
        '--out', str(diamond_output_path),
        '--id', str(int(bc.MIN_PROTEIN_IDENTITY * 100)), # '90',
        '--query-cover', str(int(bc.MIN_SORF_COVERAGE * 100)), # '90'
        '--subject-cover', str(int(bc.MIN_SORF_COVERAGE * 100)), # '90'
        '--max-target-seqs', '1',  # single best output
        '--outfmt', '6',
        '--threads', str(cfg.threads),
        '--tmpdir', str(cfg.tmp_path),  # use tmp folder
        '--block-size', '3'  # slightly increase block size for faster executions
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

    sorf_by_aa_digest = {sorf['aa_hexdigest']: sorf for sorf in sorfs}
    with diamond_output_path.open() as fh:
        for line in fh:
            (sorf_hash, cluster_id, identity, alignment_length, align_mismatches,
                align_gaps, query_start, query_end, subject_start, subject_end,
                evalue, bitscore) = line.split('\t')
            sorf = sorf_by_aa_digest[sorf_hash]
            query_cov = int(alignment_length) / len(sorf['sequence'])
            identity = float(identity) / 100
            if(query_cov >= bc.MIN_PROTEIN_COVERAGE and identity >= bc.MIN_PROTEIN_IDENTITY):
                sorf['psc'] = {
                    psc.DB_PSC_COL_UNIREF90: cluster_id,
                    'query-cov': query_cov,
                    'identity': identity
                }
                log.info(
                    'homology: contig=%s, start=%i, stop=%i, strand=%s, aa-length=%i, query-cov=%0.3f, identity=%0.3f, UniRef90=%s',
                    sorf['contig'], sorf['start'], sorf['stop'], sorf['strand'], len(sorf['sequence']), query_cov, identity, cluster_id
                )

    pscs_found = []
    pscs_not_found = []
    for sorf in sorfs:
        if('psc' in sorf):
            pscs_found.append(sorf)
        else:
            pscs_not_found.append(sorf)
    log.info('found=%i', len(pscs_found))
    return pscs_found, pscs_not_found
    
