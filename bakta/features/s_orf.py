
import logging
import subprocess as sp

from Bio.Seq import Seq

import bakta.config as cfg
import bakta.constants as bc
import bakta.utils as bu

log = logging.getLogger('features:sorf')


def extract_sorfs(contigs):
    """Predict open reading frames in mem via BioPython."""

    orfs = []
    for contig in contigs:
        dna_seq = Seq(contig['sequence'])
        for strand, seq in [('+', dna_seq), ('-', dna_seq.reverse_complement())]:  # strands +/-
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
                    orf_length = aa_end - aa_start + 1
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
                        assert sequence == test_seq, "seqs not equal! a=%s, b=%s" % (sequence, test_seq)

                        orf = {
                            'type': bc.FEATURE_SORF,
                            'contig': contig['id'],
                            'start': dna_start,
                            'stop': dna_stop,
                            'strand': strand,
                            'frame': frame + 1,
                            'sequence': sequence,
                            'aa_hash': bu.calc_aa_hash(sequence)
                        }
                        orfs.append(orf)
                        log.debug(
                            'ORF: contig=%s, start=%i, stop=%i, strand=%s, frame=%i, length=%i, seq=%s',
                            contig['id'], orf['start'], orf['stop'], strand, frame, len(sequence), sequence
                        )
                    aa_start = aa_seq.find('M', aa_start + 1)
                    if(aa_start > aa_end):
                        aa_end = aa_seq.find('*', aa_start)

    log.info('# short ORFs: # %i', len(orfs))
    return orfs


def get_feature_start(feature):
    return feature['start'] if feature['strand'] == '+' else feature['stop']


def get_feature_stop(feature):
    return feature['stop'] if feature['strand'] == '+' else feature['start']


def overlap_filter_sorfs(data, orfs_raw):
    """Filter in-mem ORFs by overlapping CDSs."""

    contig_t_rnas = {k['id']: [] for k in data['contigs']}
    for t_rna in data[bc.FEATURE_T_RNA]:
        t_rnas = contig_t_rnas[t_rna['contig']]
        t_rnas.append(t_rna)
    for tm_rna in data[bc.FEATURE_TM_RNA]:
        t_rnas = contig_t_rnas[tm_rna['contig']]
        t_rnas.append(tm_rna)

    contig_r_rnas = {k['id']: [] for k in data['contigs']}
    for r_rna in data[bc.FEATURE_R_RNA]:
        r_rnas = contig_r_rnas[r_rna['contig']]
        r_rnas.append(r_rna)

    contig_nc_rnas = {k['id']: [] for k in data['contigs']}
    for nc_rna in data[bc.FEATURE_NC_RNA]:
        nc_rnas = contig_nc_rnas[nc_rna['contig']]
        nc_rnas.append(nc_rna)
    for nc_rna in data[bc.FEATURE_NC_RNA_REGION]:
        nc_rnas = contig_nc_rnas[nc_rna['contig']]
        nc_rnas.append(nc_rna)

    contig_cdss = {k['id']: [] for k in data['contigs']}
    for cds in data[bc.FEATURE_CDS]:
        cdss = contig_cdss[cds['contig']]
        cdss.append(cds)

    contig_orfs = {k['id']: [] for k in data['contigs']}
    for orf in orfs_raw:
        orfs = contig_orfs[orf['contig']]
        orfs.append(orf)

    discarded_orfs = []
    for contig in data['contigs']:
        log.debug('filter short ORFs on contig: %s', contig['id'])
        orfs = contig_orfs[contig['id']]

        # filter CDS overlapping ORFs
        for cds in contig_cdss[contig['id']]:
            log.debug('filter short ORFs by CDS: %s%i[%i->%i]', cds['strand'], cds['frame'], cds['start'], cds['stop'])
            for orf in orfs[:]:
                if(orf['strand'] == cds['strand']):
                    if(orf['frame'] == cds['frame']):
                        if(get_feature_stop(orf) == get_feature_stop(cds)):
                            # ORF and CDS share identical stop codon position
                            orfs.remove(orf)
                            discarded_orfs.append(orf)
                        elif(orf['start'] < cds['start'] and orf['stop'] > cds['start']):
                            # in-frame ORF partially overlapping CDS upstream
                            orfs.remove(orf)
                            discarded_orfs.append(orf)
                        elif(orf['start'] >= cds['start'] and orf['stop'] <= cds['stop']):
                            # in-frame ORF completely overlapped by CDS
                            orfs.remove(orf)
                            discarded_orfs.append(orf)
                        elif(orf['start'] < cds['stop'] and orf['stop'] > cds['stop']):
                            # in-frame ORF partially overlapping CDS downstream
                            orfs.remove(orf)
                            discarded_orfs.append(orf)
                    else:
                        if(orf['start'] < cds['start'] and orf['stop'] > cds['start']):
                            # out-of-frame ORF partially overlapping CDS upstream
                            # ToDo: add max overlap threshold
                            pass
                        elif(orf['start'] >= cds['start'] and orf['stop'] <= cds['stop']):
                            # out-of-frame ORF completely overlapped by CDS
                            orfs.remove(orf)
                            discarded_orfs.append(orf)
                        elif(orf['start'] < cds['stop'] and orf['stop'] > cds['stop']):
                            # out-of-frame ORF partially overlapping CDS downstream
                            # ToDo: add max overlap threshold
                            pass
                else:
                    if(orf['frame'] == cds['frame']):
                        if(orf['start'] < cds['start'] and orf['stop'] > cds['start']):
                            # in-frame ORF partially overlapping CDS upstream
                            # ToDo: add max overlap threshold
                            pass
                        elif(orf['start'] >= cds['start'] and orf['stop'] <= cds['stop']):
                            # in-frame ORF completely overlapped by CDS
                            orfs.remove(orf)
                            discarded_orfs.append(orf)
                        elif(orf['start'] < cds['stop'] and orf['stop'] > cds['stop']):
                            # in-frame ORF partially overlapping CDS downstream
                            # ToDo: add max overlap threshold
                            pass
                    else:
                        # ToDo: add out-of-frame filters
                        if(orf['start'] < cds['start'] and orf['stop'] > cds['start']):
                            # out-frame ORF partially overlapping CDS upstream
                            # ToDo: add max overlap threshold
                            pass
                        elif(orf['start'] >= cds['start'] and orf['stop'] <= cds['stop']):
                            # out-frame ORF completely overlapped by CDS
                            orfs.remove(orf)
                            discarded_orfs.append(orf)
                        elif(orf['start'] < cds['stop'] and orf['stop'] > cds['stop']):
                            # out-frame ORF partially overlapping CDS downstream
                            # ToDo: add max overlap threshold
                            pass

        # filter rRNA overlapping ORFs
        for r_rna in contig_r_rnas[contig['id']]:
            log.debug('filter short ORFs by rRNA: %s[%i->%i]', r_rna['strand'], r_rna['start'], r_rna['stop'])
            for orf in orfs[:]:
                if(orf['start'] >= r_rna['start'] and orf['stop'] <= r_rna['stop']):
                    # ORF completely overlapped by rRNA
                    orfs.remove(orf)
                    discarded_orfs.append(orf)
                    continue

        # filter tRNA overlapping ORFs
        for t_rna in contig_t_rnas[contig['id']]:
            log.debug('filter short ORFs by tRNA: %s[%i->%i]', t_rna['strand'], t_rna['start'], t_rna['stop'])
            for orf in orfs[:]:
                if(orf['start'] < t_rna['start'] and orf['stop'] > t_rna['start']):
                    # ORF partially overlapping tRNA upstream
                    orfs.remove(orf)
                    discarded_orfs.append(orf)
                elif(orf['start'] >= t_rna['start'] and orf['stop'] <= t_rna['stop']):
                    # ORF completely overlapped by tRNA
                    orfs.remove(orf)
                    discarded_orfs.append(orf)
                elif(orf['start'] < t_rna['stop'] and orf['stop'] > t_rna['start']):
                    # ORF partially overlapping tRNA downstream
                    orfs.remove(orf)
                    discarded_orfs.append(orf)

        # filter ncRNA overlapping ORFs
        for nc_rna in contig_nc_rnas[contig['id']]:
            log.debug('filter short ORFs by ncRNA: %s[%i->%i]', nc_rna['strand'], nc_rna['start'], nc_rna['stop'])
            for orf in orfs[:]:
                if(orf['start'] < nc_rna['start'] and orf['stop'] > nc_rna['start']):
                    # ORF partially overlapping ncRNA upstream
                    # ToDo: add max overlap threshold
                    pass
                elif(orf['start'] >= nc_rna['start'] and orf['stop'] <= nc_rna['stop']):
                    # ORF completely overlapped by ncRNA
                    # ToDo: allow sORF / leader overlap, test other overlaps
                    pass
                elif(orf['start'] < nc_rna['stop'] and orf['stop'] > nc_rna['start']):
                    # ORF partially overlapping ncRNA downstream
                    # ToDo: add max overlap threshold
                    pass

    valid_orfs = []
    for orfs in contig_orfs.values():
        for orf in orfs:
            valid_orfs.append(orf)

    log.info('short ORF filter: # valid=%i, # discarded=%i', len(valid_orfs), len(discarded_orfs))
    return valid_orfs, discarded_orfs


def mark_hypotheticals(sorfs):
    for sorf in sorfs:
        ups = sorf['ups']  # must exist, otherwise filtered before
        product = ups.get('product', '')
        if(product == ''):
            psc = sorf.get('psc', None)
            if(psc):
                product = psc.get('product', '')
                if(product == '' or 'uncharacterized' in product.lower()):
                    sorf['hypothetical'] = True
            else:
                sorf['hypothetical'] = True
