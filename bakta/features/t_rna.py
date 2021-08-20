import logging
import subprocess as sp

from collections import OrderedDict

from Bio import SeqIO

import bakta.config as cfg
import bakta.constants as bc
import bakta.so as so
import bakta.utils as bu


log = logging.getLogger('T_RNA')


SO_TERMS = {
    'ala': so.SO_TRNA_ALA,
    'gln': so.SO_TRNA_GLN,
    'glu': so.SO_TRNA_GLU,
    'gly': so.SO_TRNA_GLY,
    'pro': so.SO_TRNA_PRO,
    'met': so.SO_TRNA_MET,
    'fmet': so.SO_TRNA_MET,
    'asp': so.SO_TRNA_ASP,
    'thr': so.SO_TRNA_THR,
    'val': so.SO_TRNA_VAL,
    'tyr': so.SO_TRNA_TYR,
    'cys': so.SO_TRNA_CYS,
    'ile': so.SO_TRNA_ILE,
    'ile2': so.SO_TRNA_ILE,
    'ser': so.SO_TRNA_SER,
    'leu': so.SO_TRNA_LEU,
    'trp': so.SO_TRNA_TRP,
    'lys': so.SO_TRNA_LYS,
    'asn': so.SO_TRNA_ASN,
    'arg': so.SO_TRNA_ARG,
    'his': so.SO_TRNA_HIS,
    'phe': so.SO_TRNA_PHE,
    'sec': so.SO_TRNA_SELCYS
}


def predict_t_rnas(genome, contigs_path):
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
        raise Exception(f'tRNAscan-SE error! error code: {proc.returncode}')

    trnas = {}
    contigs = {c['id']: c for c in genome['contigs']}
    with txt_output_path.open() as fh:
        for line in fh.readlines()[3:]:  # skip first 3 lines
            (contig_id, trna_id, start, stop, trna_type, anti_codon, intron_begin, bounds_end, score, note) = line.split('\t')

            start, stop, strand = int(start), int(stop), bc.STRAND_FORWARD
            if(start > stop):  # reverse
                start, stop = stop, start
                strand = bc.STRAND_REVERSE
            contig_id = contig_id.strip()  # bugfix for extra single whitespace in tRNAscan-SE output

            trna = OrderedDict()
            trna['type'] = bc.FEATURE_T_RNA
            trna['contig'] = contig_id
            trna['start'] = start
            trna['stop'] = stop
            trna['strand'] = strand
            trna['gene'] = None
            trna['product'] = 'tRNA-Xxx'
            if(trna_type != 'Undet' and trna_type != 'Sup'):
                trna['gene'] = f'{trna_type}_trna'
                trna['product'] = f'tRNA-{trna_type}'
                trna['amino_acid'] = trna_type
                trna['anti_codon'] = anti_codon.lower()

            if('pseudo' in note):
                trna['product'] = f"(pseudo) {trna['product']}"
                trna['pseudo'] = True

            trna['score'] = float(score)

            nt = bu.extract_feature_sequence(trna, contigs[contig_id])  # extract nt sequences
            trna['nt'] = nt

            trna['db_xrefs'] = []
            so_term = SO_TERMS.get(trna_type.lower(), None)
            if(so_term):
                trna['db_xrefs'].append(so_term.id)

            key = f'{contig_id}.trna{trna_id}'
            trnas[key] = trna
            log.info(
                'contig=%s, start=%i, stop=%i, strand=%s, gene=%s, product=%s, score=%1.1f, nt=[%s..%s]',
                trna['contig'], trna['start'], trna['stop'], trna['strand'], trna.get('gene', ''), trna['product'], trna['score'], nt[:10], nt[-10:]
            )

    with fasta_output_path.open() as fh:
        for record in SeqIO.parse(fh, 'fasta'):
            trna = trnas[record.id]
            nt = str(record.seq).upper()
            assert trna['nt'] == nt
            if('anti_codon' in trna and trna['amino_acid'].lower() not in ['fmet', 'ile2', 'sec', 'sup']):  # exclude fMet, Ile2 and Sec (INSDC wrong anticodon issue)
                anticodon_pos = trna['nt'].lower().find(trna['anti_codon'])
                if(anticodon_pos > -1):
                    if(trna['strand'] == bc.STRAND_FORWARD):
                        start = trna['start'] + anticodon_pos
                        stop = start + 2
                    else:
                        stop = trna['stop'] - anticodon_pos
                        start = stop - 2
                    trna['anti_codon_pos'] = (start, stop)
    trnas = list(trnas.values())
    log.info('predicted=%i', len(trnas))
    return trnas
