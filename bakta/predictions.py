
import logging
import subprocess as sp

from Bio import SeqIO
from Bio.Seq import Seq

import bakta.constants as bc
import bakta.utils as bu

log = logging.getLogger('predictions')


def predict_t_rnas(config, data, contigs_path):
    """Search for tRNA sequences."""

    txt_output_path = config['tmp'].joinpath('trna.tsv')
    fasta_output_path = config['tmp'].joinpath('trna.fasta')
    cmd = [
        'tRNAscan-SE',
        '-B',
        '--output', str(txt_output_path),
        '--fasta', str(fasta_output_path),
        '--thread', str(config['threads']),
        str(contigs_path)
    ]
    proc = sp.run(
        cmd,
        cwd=str(config['tmp']),
        env=config['env'],
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
                'db_xrefs': ['SO:0001272']
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


def predict_r_rnas(config, data, contigs_path):
    """Search for ribosomal RNA sequences."""

    output_path = config['tmp'].joinpath('rrna.tsv')

    cmd = [
        'cmsearch',
        '-Z', str(data['genome_size'] // 1000000),
        '--noali',
        '--cut_tc',
        '--notrunc',
        '--cpu', str(config['threads']),
        '--tblout', str(output_path),
        str(config['db'].joinpath('rRNA')),
        str(contigs_path)
    ]
    proc = sp.run(
        cmd,
        cwd=str(config['tmp']),
        env=config['env'],
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
    log.info('rRNAs: # %i', len(contig['rrnas']))
    return rrnas


def predict_nc_rnas(config, data, contigs_path):
    """Search for non-coding RNA sequences."""

    output_path = config['tmp'].joinpath('ncrna.tsv')
    cmd = [
        'cmsearch',
        '-Z', str(data['genome_size'] // 1000000),
        '--noali',
        '--cut_tc',
        '--notrunc',
        '--rfam',
        '--cpu', str(config['threads']),
        '--tblout', str(output_path),
        str(config['db'].joinpath('ncRNA')),
        str(contigs_path)
    ]
    proc = sp.run(
        cmd,
        cwd=str(config['tmp']),
        env=config['env'],
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
    rfam2go_path = config['db'].joinpath('rfam-go.tsv')
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
    log.info('ncRNAs: # %i', len(contig['ncRNAs']))
    return ncrnas


def predict_crispr(config, data, contigs_path):
    pass  # SO:0001459 <- Sequence Ontology


def predict_cdss(config, contigs, filtered_contigs_path):
    """Predict open reading frames with Prodigal."""

    proteins_path = config['tmp'].joinpath('proteins.faa')
    gff_path = config['tmp'].joinpath('prodigal.gff')
    cmd = [
        'prodigal',
        '-i', str(filtered_contigs_path),
        '-a', str(proteins_path),
        '-c',  # closed ends
        '-f', 'gff',  # GFF output
        '-o', str(gff_path)  # prodigal output
    ]
    proc = sp.run(
        cmd,
        cwd=str(config['tmp']),
        env=config['env'],
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
                cols = line.split('\t')
                annotations = split_gff_annotation(cols[8])
                # contig_orf_id = cols[8].split(';')[0].split('=')[1].split('_')[1]
                contig_orf_id = annotations['ID'].split['_'][1]
                cds = {
                    'type': bc.INSDC_FEATURE_CDS,
                    'contig': cols[0],
                    'start': int(cols[3]),
                    'stop': int(cols[4]),
                    'strand': cols[6],
                    'inference': 'prodigal',
                    'tmp_id': cds_id,
                    'start_type': annotations['start_type'],
                    'rbs_motif': annotations['rbs_motif']
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
            seq = str(record.seq)
            cds['sequence'] = seq
            cds['aa_hash'] = bu.calc_aa_hash(seq)

    gff_path.unlink()
    proteins_path.unlink()

    # for cds in sorted(filter(lambda k: len(k['sequence']) < 50, cdss.values()), key=lambda k: len(k['sequence']), reverse=True):
    #     print("CDS: contig=%s, start=%i, stop=%i, strand=%s, length=%i" % (cds['contig'], cds['start'], cds['stop'], cds['strand'], len(cds['sequence'])))

    log.info('CDSs: # %i', len(cdss))
    return list(cdss.values())


def split_gff_annotation(annotation_string):
    annotations = {}
    for expr in annotation_string.split(';'):
        key, value = expr.split('=')
        annotations[key] = value
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
                aa_seq = str(seq_frame.translate(table=11, stop_symbol='*', to_stop=False, cds=False, gap=None))
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
                        sequence = aa_seq[aa_start:aa_end + 1]

                        test_dna_seq = Seq(contig['sequence'][dna_start - 1:dna_stop])
                        if(strand == '-'):
                            test_dna_seq = test_dna_seq.reverse_complement()
                        test_seq = test_dna_seq.translate(table=11, stop_symbol='*', to_stop=False, cds=True, gap=None) + '*'
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

    contig_rrna = {k['id']: [] for k in data['contigs']}
    for rrna in data['rrnas']:
        rrnas = contig_rrna[rrna['contig']]
        rrnas.append(rrna)

    contig_trna = {k['id']: [] for k in data['contigs']}
    for trna in data['trnas']:
        trnas = contig_trna[trna['contig']]
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
                        elif(orf['start'] < cds['stop'] and orf['stop'] > cds['start']):
                            # in-frame ORF partially overlapping CDS downstream
                            orfs.remove(orf)
                            discarded_orfs.append(orf)
                    else:
                        pass
                else:
                    if(orf['frame'] == cds['frame']):
                        pass
                    else:
                        pass

        # filter rRNA overlapping ORFs
        for rrna in contig_cdss[contig['id']]:
            log.debug('filter ORFs by rRNA: %s[%i->%i]', rrna['strand'], rrna['start'], rrna['stop'])
            for orf in orfs[:]:
                if(orf['start'] >= rrna['start'] and orf['stop'] <= rrna['stop']):
                    # ORF completely overlapped by rRNA
                    orfs.remove(orf)
                    discarded_orfs.append(orf)
                    continue

        # filter tRNA overlapping ORFs
        for trna in contig_cdss[contig['id']]:
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

    log.info('ORF filter: # valid=%i, # discarded=%i', (len(valid_orfs), len(discarded_orfs)))
    return valid_orfs, discarded_orfs
