import atexit
import logging
import multiprocessing as mp
import os
import sys

from typing import Sequence
from pathlib import Path

import bakta
import bakta.constants as bc
import bakta.config as cfg
import bakta.db as db
import bakta.utils as bu
import bakta.expert.amrfinder as exp_amr
import bakta.expert.protein_sequences as exp_aa_seq
import bakta.features.annotation as anno
import bakta.features.orf as orf
import bakta.features.cds as feat_cds
import bakta.io.fasta as fasta
import bakta.io.tsv as tsv
import bakta.ups as ups
import bakta.ips as ips
import bakta.psc as psc
import bakta.pscc as pscc


log = logging.getLogger('PROTEINS')


def main():
    # parse options and arguments
    parser = bu.init_parser(sub_command='_proteins')
    parser.add_argument('input', metavar='<input>', help='Protein sequences in (zipped) fasta format')
    
    arg_group_io = parser.add_argument_group('Input / Output')
    arg_group_io.add_argument('--db', '-d', action='store', default=None, help='Database path (default = <bakta_path>/db). Can also be provided as BAKTA_DB environment variable.')
    arg_group_io.add_argument('--output', '-o', action='store', default=os.getcwd(), help='Output directory (default = current working directory)')
    arg_group_io.add_argument('--prefix', '-p', action='store', default=None, help='Prefix for output files')
    
    arg_group_annotation = parser.add_argument_group('Annotation')
    arg_group_annotation.add_argument('--proteins', action='store', default=None, dest='proteins', help='Fasta file of trusted protein sequences for annotation')
    
    arg_group_general = parser.add_argument_group('Runtime & auxiliary options')
    arg_group_general.add_argument('--help', '-h', action='help', help='Show this help message and exit')
    arg_group_general.add_argument('--verbose', '-v', action='store_true', help='Print verbose information')
    arg_group_general.add_argument('--threads', '-t', action='store', type=int, default=mp.cpu_count(), help='Number of threads to use (default = number of available CPUs)')
    arg_group_general.add_argument('--tmp-dir', action='store', default=None, dest='tmp_dir', help='Location for temporary files (default = system dependent auto detection)')
    arg_group_general.add_argument('--version', '-V', action='version', version=f'%(prog)s {bakta.__version__}')
    args = parser.parse_args()

    ############################################################################
    # Setup logging
    ############################################################################
    prefix = args.prefix if args.prefix else Path(args.input).stem
    output_path = cfg.check_output_path(args)
    
    bu.setup_logger(output_path, prefix, args)
    log.info('prefix=%s', prefix)
    log.info('output=%s', output_path)

    ############################################################################
    # Checks and configurations
    # - check parameters and setup global configuration
    # - test database
    # - test binary dependencies
    ############################################################################
    try:
        if(args.input == ''):
            raise ValueError('File path argument must be non-empty')
        aa_path = Path(args.input).resolve()
        cfg.check_readability('proteins', aa_path)
        cfg.check_content_size('proteins', aa_path)
    except:
        log.error('provided input proteins file not valid! path=%s', args.input)
        sys.exit(f'ERROR: input proteins file ({args.input}) not valid!')
    log.info('input-path=%s', aa_path)
    
    cfg.check_db_path(args)
    cfg.db_info = db.check(cfg.db_path)
    cfg.check_tmp_path(args)
    cfg.check_user_proteins(args)
    cfg.check_threads(args)
    cfg.skip_cds = False  # circumvent main config setup

    bu.test_dependencies()
    atexit.register(bu.cleanup, log, cfg.tmp_path)  # register cleanup exit hook

    ############################################################################
    # Setup logging
    ############################################################################
    try:
        aas = fasta.import_contigs(aa_path, False, False)
        log.info('imported sequences=%i', len(aas))
        print(f'\timported: {len(aas)}')
    except:
        log.error('wrong file format or unallowed characters in amino acid sequences!', exc_info=True)
        sys.exit('ERROR: wrong file format or unallowed characters in amino acid sequences!')
    mock_start = 1
    for aa in aas:  # rename and mock feature attributes to reuse existing functions
        aa['type'] = bc.FEATURE_CDS
        aa['aa'] = aa['sequence']
        aa['locus'] = aa['id']
        aa['contig'] = 'mock'
        aa['start'] = mock_start
        aa['stop'] = -1
        aa['strand'] = bc.STRAND_UNKNOWN
        aa['frame'] = 1
        mock_start += 100
    print('annotate protein sequences...')
    annotate_aa(aas)
    
    print('write results...')
    annotations_path = output_path.joinpath(f'{prefix}.tsv')
    header_columns = ['ID', 'Length', 'Gene', 'Product', 'EC', 'GO', 'COG', 'RefSeq', 'UniParc', 'UniRef']
    print(f'\tfull annotations (TSV): {annotations_path}')
    tsv.write_features(aas, header_columns, map_aa_columns, annotations_path)
    hypotheticals_path = output_path.joinpath(f'{prefix}.hypotheticals.tsv')
    header_columns = ['ID', 'Length', 'Mol Weight [kDa]', 'Iso El. Point', 'Pfam hits']
    hypotheticals = hypotheticals = [aa for aa in aas if 'hypothetical' in aa]
    print(f'\tinformation on hypotheticals (TSV): {hypotheticals_path}')
    tsv.write_features(hypotheticals, header_columns, map_hypothetical_columns, hypotheticals_path)
    aa_output_path = output_path.joinpath(f'{prefix}.faa')
    print(f'\tannotated sequences (Fasta): {aa_output_path}')
    fasta.write_faa(aas, aa_output_path)


def map_aa_columns(feat: dict) -> Sequence[str]:
    gene = feat.get('gene', None)
    if(gene is None):
        gene = ''
    return [
        feat['id'],
        str(feat['length']),
        gene,
        feat['product'],
        ','.join([k.replace('EC:', '') for k in feat['db_xrefs'] if 'EC:' in k]),
        ','.join([k for k in feat['db_xrefs'] if 'GO:' in k]),
        ','.join([k.replace('COG:', '') for k in feat['db_xrefs'] if 'COG:' in k]),
        ','.join([k.replace('RefSeq:', '') for k in feat['db_xrefs'] if 'RefSeq:' in k]),
        ','.join([k.replace('UniParc:', '') for k in feat['db_xrefs'] if 'UniParc:' in k]),
        ','.join([k.replace('UniRef:', '') for k in feat['db_xrefs'] if 'UniRef' in k])
    ]


def map_hypothetical_columns(feat: dict) -> Sequence[str]:
    return [
        feat['id'],
        str(feat['length']),
        f"{(feat['seq_stats']['molecular_weight']/1000):.1f}" if feat['seq_stats']['molecular_weight'] else 'NA'
        f"{feat['seq_stats']['isoelectric_point']:.1f}" if feat['seq_stats']['isoelectric_point'] else 'NA'
        ','.join([k.replace('PFAM:', '') for k in feat['db_xrefs'] if 'PFAM:' in k])
    ]


def annotate_aa(aas: Sequence[dict]):
    for aa in aas:
        aa['aa_digest'], aa['aa_hexdigest'] = bu.calc_aa_hash(aa['sequence'])
    log.debug('lookup AA UPS/IPS')
    aas_ups, aas_not_found = ups.lookup(aas)
    aas_ips, tmp = ips.lookup(aas_ups)
    aas_not_found.extend(tmp)
    print(f'\tdetected IPSs: {len(aas_ips)}')
    if(len(aas_not_found) > 0):
        log.debug('search CDS PSC')
        aas_psc, aas_pscc, aas_not_found = psc.search(aas_not_found)
        print(f'\tfound PSCs: {len(aas_psc)}')
        print(f'\tfound PSCCs: {len(aas_pscc)}')
    print('\tlookup annotations...')
    log.debug('lookup CDS PSCs')
    psc.lookup(aas)  # lookup PSC info
    pscc.lookup(aas)  # lookup PSCC inf
    print('\tconduct expert systems...')  # conduct expert systems annotation
    aa_path = cfg.tmp_path.joinpath('aa.faa')
    orf.write_internal_faa(aas, aa_path)
    log.debug('conduct expert system: amrfinder')
    cfg.translation_table = 11
    expert_amr_found = exp_amr.search(aas, aa_path)
    print(f'\t\tamrfinder: {len(expert_amr_found)}')
    log.debug('conduct expert system: aa seqs')
    diamond_db_path = cfg.db_path.joinpath('expert-protein-sequences.dmnd')
    expert_aa_found = exp_aa_seq.search(aas, aa_path, 'expert_proteins', diamond_db_path)
    print(f'\t\tprotein sequences: {len(expert_aa_found)}')
    if(cfg.user_proteins):
        log.debug('conduct expert system: user aa seqs')
        user_aa_path = cfg.tmp_path.joinpath('user-proteins.faa')
        exp_aa_seq.write_user_protein_sequences(user_aa_path)
        user_aa_found = exp_aa_seq.search(aas, aa_path, 'user_proteins', user_aa_path)
        print(f'\t\tuser protein sequences: {len(user_aa_found)}')
    print('\tcombine annotations and mark hypotheticals...')
    log.debug('combine CDS annotations')
    for aa in aas:
        anno.combine_annotation(aa)  # combine IPS & PSC annotations and mark hypothetical
    log.debug('analyze hypotheticals')
    hypotheticals = [aa for aa in aas if 'hypothetical' in aa]
    if(len(hypotheticals) > 0):
        print(f'\tanalyze hypothetical proteins: {len(hypotheticals)}')
        pfam_hits = feat_cds.predict_pfam(hypotheticals)
        print(f"\tdetected Pfam hits: {len(pfam_hits)}")
        feat_cds.analyze_proteins(hypotheticals)
        print('\tcalculated proteins statistics')


if __name__ == '__main__':
    main()
