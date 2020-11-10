import logging
import subprocess as sp
from collections import OrderedDict

from Bio.Seq import Seq

import bakta.config as cfg
import bakta.constants as bc
import bakta.utils as bu
import bakta.psc as psc

log = logging.getLogger('ORF')

def detect_spurious(orfs):
    """Detect spurious ORFs with AntiFam"""
    orf_fasta_path = cfg.tmp_path.joinpath('sorf.faa')
    with orf_fasta_path.open(mode='w') as fh:
        for orf in orfs:
            fh.write(f">{orf['aa_hexdigest']}\n{orf['sequence']}\n")
    
    output_path = cfg.tmp_path.joinpath('cds.spurious.hmm.tsv')
    cmd = [
        'hmmsearch',
        '--noali',
        '--cut_ga',  # use gathering cutoff
        '--tblout', str(output_path),
        '-Z', str(len(orfs)),
        '--cpu', str(cfg.threads),
        str(cfg.db_path.joinpath('antifam')),
        str(orf_fasta_path)
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
        log.warning('spurious ORF detection failed! hmmsearch-error-code=%d', proc.returncode)
        raise Exception(f'hmmsearch error! error code: {proc.returncode}')

    discarded_orfs = []
    orf_by_aa_digest = {orf['aa_hexdigest']: orf for orf in orfs}
    with output_path.open() as fh:
        for line in fh:
            if(line[0] != '#'):
                (aa_hexdigest, _, subject_name, subject_id, evalue, bitscore, _) = line.strip().split(maxsplit=6)
                orf = orf_by_aa_digest[aa_hexdigest]
                evalue = float(evalue)
                bitscore = float(bitscore)
                if(evalue > 1E-5):
                    log.debug(
                        'discard low spurious E value: contig=%s, start=%i, stop=%i, strand=%s, subject=%s, evalue=%f, bitscore=%f',
                        orf['contig'], orf['start'], orf['stop'], orf['strand'], subject_name, evalue, bitscore
                    )
                else:
                    discard = OrderedDict()
                    discard['type'] = bc.DISCARD_TYPE_SPURIOUS
                    discard['description'] = f'(partial) homology to spurious sequence HMM (AntiFam:{subject_id})'
                    discard['score'] = bitscore
                    discard['evalue'] = evalue
                    
                    orf['discarded'] = discard
                    discarded_orfs.append(orf)
                    log.info(
                        'discard spurious: contig=%s, start=%i, stop=%i, strand=%s, homology=%s, evalue=%f, bitscore=%f',
                        orf['contig'], orf['start'], orf['stop'], orf['strand'], subject_name, evalue, bitscore
                    )
    log.info('discarded=%i', len(discarded_orfs))
    return discarded_orfs