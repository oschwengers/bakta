
import logging
import subprocess as sp

from Bio import SeqIO

import bakta.config as cfg
import bakta.constants as bc

log = logging.getLogger('features:t_rna')


SO_TERMS = {
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
        log.warning('tRNAs failed! tRNAscan-SE-error-code=%d', proc.returncode)
        raise Exception("tRNAscan-SE error! error code: %i" % proc.returncode)

    trnas = {}
    with txt_output_path.open() as fh:
        for line in fh.readlines()[3:]:  # skip first 3 lines
            (contig, trna_id, start, stop, trna_type, anti_codon, intron_begin, bounds_end, score, note) = line.split('\t')

            if(trna_type == 'Undet'):
                trna_type = ''

            start, stop, strand = int(start), int(stop), '+'
            if(start > stop):  # reverse
                start, stop = stop, start
                strand = '-'

            contig = contig.strip()  # bugfix for extra single whitespace in tRNAscan-SE output

            trna = {
                'type': bc.FEATURE_T_RNA,
                'contig': contig.strip(),
                'start': start,
                'stop': stop,
                'strand': strand,
                'score': float(score),
                'db_xrefs': []
            }

            if(trna_type == 'Undet'):
                trna['gene'] = ''
                trna['product'] = 'tRNA-Xxx'
            else:
                trna['gene'] = "%s_trna" % trna_type
                trna['product'] = "tRNA-%s" % trna_type
                trna['anti_codon'] = anti_codon
                trna['notes'] = ["tRNA-%s(%s)" % (trna_type, anti_codon)]
                so_term = SO_TERMS.get(trna_type.lower(), None)
                if(so_term):
                    trna['db_xrefs'].append(so_term)
            
            if('pseudo' in note):
                trna['gene'] = ''
                trna['product'] = '(pseudo) %s' % trna['product']
                trna['pseudo'] = True

            key = "%s.trna%s" % (contig, trna_id)
            trnas[key] = trna
            log.info(
                'contig=%s, gene=%s, start=%i, stop=%i, strand=%s',
                trna['contig'], trna.get('gene', ''), trna['start'], trna['stop'], trna['strand']
            )

    with fasta_output_path.open() as fh:
        for record in SeqIO.parse(fh, 'fasta'):
            trna = trnas[record.id]
            trna['sequence'] = str(record.seq)
    trnas = list(trnas.values())
    log.info('# %i', len(trnas))
    return trnas