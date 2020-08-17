
import logging
import subprocess as sp

from Bio import SeqIO
from Bio.Seq import Seq

import bakta.config as cfg
import bakta.constants as bc
import bakta.utils as bu

log = logging.getLogger('predictions')


def predict_t_rnas(data, contigs_path):
    """Search for tRNA sequences."""

    txt_output_path = cfg.tmp_path.joinpath('trna.tsv')
    fasta_output_path = cfg.tmp_path.joinpath('trna.fasta')
    cmd = [
        'tRNAscan-SE',
        '-B',
        '--output', str(txt_output_path),
        '--fasta', str(fasta_output_path),
        '--thread', str(cfg.threads),
        str(contigs_path)
    ]
    proc = sp.run(
        cmd,
        cwd=str(cfg.tmp_path),
        env=cfg.env,
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        universal_newlines=True
    )
    if(proc.returncode != 0):
        log.debug(
            'tRNAs: cmd=%s, stdout=\'%s\', stderr=\'%s\'',
            cmd, proc.stdout, proc.stderr
        )
        log.warning('tRNAs failed! tRNAscan-SE-error-code=%d', proc.returncode)
        raise Exception("tRNAscan-SE error! error code: %i" % proc.returncode)

    trnas = {}
    with txt_output_path.open() as fh:
        for line in fh.readlines()[3:]:  # skip first 3 lines
            (contig, trna_id, start, stop, type, anti_codon, intron_begin, bounds_end, score, note) = line.split('\t')

            if(type == 'Undet'):
                type = ''

            start, stop, strand = int(start), int(stop), '+'
            if(start > stop):  # reverse
                start, stop = stop, start
                strand = '-'

            contig = contig.strip()  # bugfix for extra single whitespace in tRNAscan-SE output

            so_terms = {
                'ala': 'SO:0000254',  # alanyl_tRNA
                'gln': 'SO:0000261',  # glycyl_tRNA
                'glu': 'SO:0000259',  # glutaminyl_tRNA
                'gly': 'SO:0000260',  # glutamyl_tRNA
                'pro': 'SO:0000268',  # prolyl_tRNA
                'met': 'SO:0000266',  # methionyl_tRNA
                'asp': 'SO:0000256',  # asparaginyl_tRNA
                'thr': 'SO:0000270',  # threonyl_tRNA
                'val': 'SO:0000273',  # valyl_tRNA
                'tyr': 'SO:0000272',  # tyrosyl_tRNA
                'cys': 'SO:0000258',  # cysteinyl_tRNA
                'ile': 'SO:0000263',  # isoleucyl_tRNA
                'ser': 'SO:0000269',  # seryl_tRNA
                'leu': 'SO:0000264',  # leucyl_tRNA
                'trp': 'SO:0000271',  # tryptophanyl_tRNA
                'lys': 'SO:0000265',  # lysyl_tRNA
                'asn': 'SO:0000257',  # aspartyl_tRNA
                'arg': 'SO:0001036',  # arginyl_tRNA
                'his': 'SO:0000262',  # histidyl_tRNA
                'phe': 'SO:0000267',  # phenylalanyl_tRNA
                'selcys': 'SO:0005857'  # selenocysteinyl_tRNA
            }

            trna = {
                'type': bc.INSDC_FEATURE_T_RNA,
                'gene': "%s_trna" % type,
                'product': "tRNA-%s" % type,
                'contig': contig.strip(),
                'start': start,
                'stop': stop,
                'strand': strand,
                'inference': 'tRNAscan-SE',
                'score': float(score),
                'pseudo': 'pseudo' in note,
                'notes': ["tRNA-%s(%s)" % (type, anti_codon)],
                'db_xrefs': so_terms.get(type.lower(), [])
            }

            key = "%s.trna%s" % (contig, trna_id)
            log.debug("key=%s", key)
            trnas[key] = trna
            log.info(
                'tRNA: contig=%s, gene=%s, start=%i, stop=%i, strand=%s',
                trna['contig'], trna['gene'], trna['start'], trna['stop'], trna['strand']
            )

    with fasta_output_path.open() as fh:
        for record in SeqIO.parse(fh, 'fasta'):
            log.debug("key=%s", record.id)
            trna = trnas[record.id]
            trna['sequence'] = str(record.seq)
    trnas = list(trnas.values())
    log.info('tRNAs: # %i', len(trnas))
    return trnas


def predict_tm_rnas(data, contigs_path):
    """Search for tmRNA sequences."""

    contigs = {c['id']: c for c in data['contigs']}
    txt_output_path = cfg.tmp_path.joinpath('tmrna.tsv')
    cmd = [
        'aragorn',
        '-m',  # detect tmRNAs
        '-gcbact',
        '-w',  # batch mode
        '-o', str(txt_output_path),
        str(contigs_path)
    ]
    if(cfg.complete == True):
        cmd.append('-c')  # complete circular sequence(s)
    else:
        cmd.append('-l')  # linear sequence(s)

    proc = sp.run(
        cmd,
        cwd=str(cfg.tmp_path),
        env=cfg.env,
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        universal_newlines=True
    )
    if(proc.returncode != 0):
        log.debug(
            'tmRNAs: cmd=%s, stdout=\'%s\', stderr=\'%s\'',
            cmd, proc.stdout, proc.stderr
        )
        log.warning('tmRNAs failed! aragorn-error-code=%d', proc.returncode)
        raise Exception("aragorn error! error code: %i" % proc.returncode)

    tmrnas = []
    with txt_output_path.open() as fh:
        for line in fh:
            line = line.strip()
            cols = line.split()
            if(line[0] == '>'):
                contig = cols[0][1:]
            elif( len(cols) == 5 ):
                (nr, type, location, tag_location, tag_aa) = line.split()
                strand = '+'
                if(location[0] == 'c'):
                    strand = '-'
                    location = location[1:]
                (start, stop) = location[1:-1].split(',')
                start = int(start)
                stop = int(stop)

                # extract sequence
                seq = contigs[contig]['sequence'][start:stop]
                if(strand == '-'):
                    seq = Seq(seq).reverse_complement()
                tmrna = {
                    'type': bc.INSDC_FEATURE_TM_RNA,
                    'gene': 'ssrA',
                    'product': 'transfer-messenger RNA, SsrA',
                    'contig': contig,
                    'start': start,
                    'stop': stop,
                    'strand': strand,
                    'sequence': seq,
                    'inference': 'aragorn',
                    'db_xrefs': ['SO:0000584']
                }
                tmrnas.append(tmrna)
                log.info(
                    'tmRNA: contig=%s, gene=%s, start=%i, stop=%i, strand=%s',
                    tmrna['contig'], tmrna['gene'], tmrna['start'], tmrna['stop'], tmrna['strand']
                )
    log.info('tmRNAs: # %i', len(tmrnas))
    return tmrnas


def predict_r_rnas(data, contigs_path):
    """Search for ribosomal RNA sequences."""

    output_path = cfg.tmp_path.joinpath('rrna.tsv')

    cmd = [
        'cmsearch',
        '-Z', str(data['genome_size'] // 1000000),
        '--noali',
        '--cut_tc',
        '--notrunc',
        '--cpu', str(cfg.threads),
        '--tblout', str(output_path),
        str(cfg.db_path.joinpath('rRNA')),
        str(contigs_path)
    ]
    proc = sp.run(
        cmd,
        cwd=str(cfg.tmp_path),
        env=cfg.env,
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        universal_newlines=True
    )
    if(proc.returncode != 0):
        log.warning('rRNAs failed! cmscan-error-code=%d', proc.returncode)
        log.debug(
            'rRNAs: cmd=%s, stdout=\'%s\', stderr=\'%s\'',
            cmd, proc.stdout, proc.stderr
        )
        raise Exception("cmsearch error! error code: %i" % proc.returncode)

    rrnas = []
    with output_path.open() as fh:
        for line in fh:
            if(line[0] != '#'):
                (contig, accession, subject, subject_id, mdl, mdl_from, mdl_to,
                    start, stop, strand, trunc, passed, gc, bias, score, evalue,
                    inc, description) = line.strip().split()
                db_xrefs = ['GO:0005840', 'GO:0003735']
                if(subject_id == 'RF00001'):
                    rrna_tag = '5S'
                    db_xrefs += ['RFAM:RF00001', 'SO:0000652']
                elif(subject_id == 'RF00177'):
                    rrna_tag = '16S'
                    db_xrefs += ['RFAM:RF00177', 'SO:0001000']
                elif(subject_id == 'RF02541'):
                    rrna_tag = '23S'
                    db_xrefs += ['RFAM:RF02541', 'SO:0001001']

                rrna = {
                    'type': bc.INSDC_FEATURE_R_RNA,
                    'gene': "%s_rrna" % rrna_tag,
                    'product': "%s ribosomal RNA" % rrna_tag,
                    'contig': contig,
                    'start': int(start),
                    'stop': int(stop),
                    'strand': strand,
                    'inference': 'infernal',
                    'score': float(score),
                    'evalue': float(evalue),
                    'db_xrefs': db_xrefs
                }
                rrnas.append(rrna)
                log.debug(
                    'rRNA: contig=%s, gene=%s, start=%i, stop=%i, strand=%s',
                    rrna['contig'], rrna['gene'], rrna['start'], rrna['stop'], rrna['strand']
                )
    log.info('rRNAs: # %i', len(rrnas))
    return rrnas


def predict_nc_rnas(data, contigs_path):
    """Search for non-coding RNA sequences."""

    output_path = cfg.tmp_path.joinpath('ncrna.tsv')
    cmd = [
        'cmsearch',
        '-Z', str(data['genome_size'] // 1000000),
        '--noali',
        '--cut_tc',
        '--notrunc',
        '--rfam',
        '--cpu', str(cfg.threads),
        '--tblout', str(output_path),
        str(cfg.db_path.joinpath('ncRNA')),
        str(contigs_path)
    ]
    proc = sp.run(
        cmd,
        cwd=str(cfg.tmp_path),
        env=cfg.env,
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        universal_newlines=True
    )
    if(proc.returncode != 0):
        log.warning('ncRNAs failed! cmscan-error-code=%d', proc.returncode)
        log.debug(
            'ncRNAs: cmd=%s, stdout=\'%s\', stderr=\'%s\'',
            cmd, proc.stdout, proc.stderr
        )
        raise Exception("cmsearch error! error code: %i" % proc.returncode)

    rfam2go = {}
    rfam2go_path = cfg.db_path.joinpath('rfam-go.tsv')
    with rfam2go_path.open() as fh:
        for line in fh:
            (rfam, go) = line.split('\t')
            if(rfam in rfam2go):
                rfam2go[rfam].append(go)
            else:
                rfam2go[rfam] = [go]

    ncrnas = []
    with output_path.open() as fh:
        for line in fh:
            if(line[0] != '#'):
                (contig, accession, subject, subject_id, mdl, mdl_from, mdl_to,
                    start, stop, strand, trunc, passed, gc, bias, score, evalue,
                    inc, description) = line.strip().split()
                rfam_id = "RFAM:%s" % subject_id
                db_xrefs = [rfam_id, 'SO:0001263']
                if(rfam_id in rfam2go):
                    db_xrefs += rfam2go[rfam_id]
                ncrna = {
                    'type': bc.INSDC_FEATURE_NC_RNA,
                    'gene': subject,
                    'contig': contig,
                    'start': int(start),
                    'stop': int(stop),
                    'strand': strand,
                    'inference': 'infernal',
                    'score': float(score),
                    'evalue': float(evalue),
                    'db_xrefs': db_xrefs
                }
                ncrnas.append(ncrna)
                log.debug(
                    'ncRNA: contig=%s, gene=%s, start=%i, stop=%i, strand=%s',
                    ncrna['contig'], ncrna['gene'], ncrna['start'], ncrna['stop'], ncrna['strand']
                )
    log.info('ncRNAs: # %i', len(ncrnas))
    return ncrnas


def predict_crispr(data, contigs_path):
    pass  # SO:0001459 <- Sequence Ontology


def predict_cdss(contigs, filtered_contigs_path):
    """Predict open reading frames with Prodigal."""

    proteins_path = cfg.tmp_path.joinpath('proteins.faa')
    gff_path = cfg.tmp_path.joinpath('prodigal.gff')
    cmd = [
        'prodigal',
        '-i', str(filtered_contigs_path),
        '-a', str(proteins_path),
        '-f', 'gff',  # GFF output
        '-o', str(gff_path)  # prodigal output
    ]
    if(cfg.complete == False):
        cmd.append('-c')  # closed ends
    
    proc = sp.run(
        cmd,
        cwd=str(cfg.tmp_path),
        env=cfg.env,
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        universal_newlines=True
    )
    if(proc.returncode != 0):
        log.warning(
            'ORFs failed! prodigal-error-code=%d', proc.returncode
        )
        log.debug(
            'ORFs: cmd=%s stdout=\'%s\', stderr=\'%s\'',
            cmd, proc.stdout, proc.stderr
        )
        raise Exception("prodigal error! error code: %i" % proc.returncode)

    # parse orfs
    # TODO: replace code by BioPython GFF3 parser
    contigs = {k['id']: k for k in contigs}
    cdss = {}
    cds_id = 1
    with gff_path.open() as fh:
        for line in fh:
            if(line[0] != '#'):
                (contig, inference, _, start, stop, score, strand, _, annotations_raw) = line.strip().split('\t')
                gff_annotations = split_gff_annotation(annotations_raw)
                contig_orf_id = gff_annotations['ID'].split('_')[1]
                cds = {
                    'type': bc.INSDC_FEATURE_CDS,
                    'contig': contig,
                    'start': int(start),
                    'stop': int(stop),
                    'strand': strand,
                    'inference': inference,
                    'tmp_id': cds_id,
                    'start_type': gff_annotations['start_type'],
                    'rbs_motif': gff_annotations['rbs_motif']
                }
                cds_id += 1
                if(cds['strand'] == '+'):
                    cds['frame'] = (cds['start'] - 1) % 3 + 1
                else:
                    cds['frame'] = (contigs[cds['contig']]['length'] - cds['stop']) % 3 + 1
                cdss["%s_%s" % (cds['contig'], contig_orf_id)] = cds
                log.info(
                    'cds: contig=%s, start=%i, stop=%i, strand=%s',
                    cds['contig'], cds['start'], cds['stop'], cds['strand']
                )

    # extract translated orf sequences
    with proteins_path.open() as fh:
        for record in SeqIO.parse(fh, 'fasta'):
            cds = cdss[record.id]
            seq = str(record.seq)[:-1]  # discard trailing asterisk
            cds['sequence'] = seq
            cds['aa_hash'] = bu.calc_aa_hash(seq)

    gff_path.unlink()
    proteins_path.unlink()
    log.info('CDSs: # %i', len(cdss))
    return list(cdss.values())


def split_gff_annotation(annotation_string):
    annotations = {}
    for expr in annotation_string.split(';'):
        if(expr != ''):
            try:
                key, value = expr.split('=')
                annotations[key] = value
            except:
                log.error('expr=%s' % expr)
    return annotations


def extract_orfs(contigs):
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

                # print("contig=%s, strand=%s, frame=%s, length=%i, dna=%s" % (contig['id'], strand, frame, len(seq_frame), seq_frame[:200]))
                aa_seq = str(seq_frame.translate(table=11, stop_symbol='*', to_stop=False, cds=False))
                # print("\taa length=%i" % len(aa_seq))
                # print("\taa seq=%s" % aa_seq[:200])
                aa_start = aa_seq.find('M')
                aa_end = aa_seq.find('*', aa_start)
                while aa_start > -1 and aa_end > -1:
                    orf_length = aa_end - aa_start + 1
                    # print("\tM: %i, *: %i" % (aa_start, aa_end))
                    if(orf_length >= bc.MIN_ORF_LENGTH and orf_length <= bc.MAX_ORF_LENGTH):  # get all CDS starts (M)
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
                            'type': 'orf',
                            'contig': contig['id'],
                            'start': dna_start,
                            'stop': dna_stop,
                            'strand': strand,
                            'frame': frame + 1,
                            'inference': 'bakta',
                            'sequence': sequence,
                            'aa_hash': bu.calc_aa_hash(sequence)
                        }
                        orfs.append(orf)
                        # print(
                        #     "\tORF: contig=%s, strand=%s, frame=%i, start=%i, stop=%i, seq=%s..%s, length=%i" %
                        #     (contig['id'], strand, frame, orf['start'], orf['stop'], sequence[:min(30, len(sequence) - 5)], sequence[-3:], len(sequence))
                        # )
                    aa_start = aa_seq.find('M', aa_start + 1)
                    if(aa_start > aa_end):
                        aa_end = aa_seq.find('*', aa_start)

    log.info('ORFs: # %i', len(orfs))
    return orfs


def get_feature_start(feature):
    return feature['start'] if feature['strand'] == '+' else feature['stop']


def get_feature_stop(feature):
    return feature['stop'] if feature['strand'] == '+' else feature['start']


def overlap_filter_orfs(data, orfs_raw):
    """Filter in-mem ORFs by overlapping CDSs."""
    contig_cdss = {k['id']: [] for k in data['contigs']}
    for cds in data['cdss']:
        cdss = contig_cdss[cds['contig']]
        cdss.append(cds)

    contig_rrnas = {k['id']: [] for k in data['contigs']}
    for rrna in data['r_rnas']:
        rrnas = contig_rrnas[rrna['contig']]
        rrnas.append(rrna)

    contig_trnas = {k['id']: [] for k in data['contigs']}
    for trna in data['t_rnas']:
        trnas = contig_trnas[trna['contig']]
        trnas.append(trna)

    contig_orfs = {k['id']: [] for k in data['contigs']}
    for orf in orfs_raw:
        orfs = contig_orfs[orf['contig']]
        orfs.append(orf)

    discarded_orfs = []
    for contig in data['contigs']:
        log.debug('filter ORFs on contig: %s', contig['id'])
        orfs = contig_orfs[contig['id']]

        # filter CDS overlapping ORFs
        for cds in contig_cdss[contig['id']]:
            log.debug('filter ORFs by CDS: %s%i[%i->%i]', cds['strand'], cds['frame'], cds['start'], cds['stop'])
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

        # filter rRNA overlapping ORFs
        for rrna in contig_rrnas[contig['id']]:
            log.debug('filter ORFs by rRNA: %s[%i->%i]', rrna['strand'], rrna['start'], rrna['stop'])
            for orf in orfs[:]:
                if(orf['start'] >= rrna['start'] and orf['stop'] <= rrna['stop']):
                    # ORF completely overlapped by rRNA
                    orfs.remove(orf)
                    discarded_orfs.append(orf)
                    continue

        # filter tRNA overlapping ORFs
        for trna in contig_trnas[contig['id']]:
            log.debug('filter ORFs by tRNA: %s[%i->%i]', trna['strand'], trna['start'], trna['stop'])
            for orf in orfs[:]:
                if(orf['start'] < trna['start'] and orf['stop'] > trna['start']):
                    # in-frame ORF partially overlapping tRNA upstream
                    orfs.remove(orf)
                    discarded_orfs.append(orf)
                elif(orf['start'] >= trna['start'] and orf['stop'] <= trna['stop']):
                    # in-frame ORF completely overlapped by tRNA
                    orfs.remove(orf)
                    discarded_orfs.append(orf)
                elif(orf['start'] < trna['stop'] and orf['stop'] > trna['start']):
                    # in-frame ORF partially overlapping tRNA downstream
                    orfs.remove(orf)
                    discarded_orfs.append(orf)

    valid_orfs = []
    for orfs in contig_orfs.values():
        for orf in orfs:
            valid_orfs.append(orf)

    log.info('ORF filter: # valid=%i, # discarded=%i', len(valid_orfs), len(discarded_orfs))
    return valid_orfs, discarded_orfs
