import atexit
import json
import logging
import os
import sys

from datetime import datetime
from pathlib import Path

from xopen import xopen

import bakta
import bakta.constants as bc
import bakta.config as cfg
import bakta.utils as bu
import bakta.io.fasta as fasta
import bakta.io.tsv as tsv
import bakta.io.gff as gff
import bakta.io.insdc as insdc
import bakta.plot as plot


log = logging.getLogger('IO')


def main():
    # parse options and arguments
    parser = bu.init_parser(sub_command='_proteins')
    parser.add_argument('input', metavar='<input>', help='Bakta annotations in JSON format')
    
    arg_group_io = parser.add_argument_group('Input / Output')
    arg_group_io.add_argument('--output', '-o', action='store', default=os.getcwd(), help='Output directory (default = current working directory)')
    arg_group_io.add_argument('--prefix', '-p', action='store', default=None, help='Prefix for output files')
    arg_group_io.add_argument('--force', '-f', action='store_true', help='Force overwriting existing output folder')
    
    arg_group_general = parser.add_argument_group('General')
    arg_group_general.add_argument('--help', '-h', action='help', help='Show this help message and exit')
    arg_group_general.add_argument('--verbose', '-v', action='store_true', help='Print verbose information')
    arg_group_general.add_argument('--debug', action='store_true', help='Run Bakta in debug mode. Temp data will not be removed.')
    arg_group_general.add_argument('--tmp-dir', action='store', default=None, dest='tmp_dir', help='Location for temporary files (default = system dependent auto detection)')
    arg_group_general.add_argument('--version', '-V', action='version', version=f'%(prog)s {cfg.version}')
    args = parser.parse_args()

    ############################################################################
    # Setup logging
    ############################################################################
    cfg.prefix = args.prefix if args.prefix else Path(args.input).stem
    output_path = cfg.check_output_path(args.output, args.force)
    cfg.force = args.force
    log.info('force=%s', args.force)
    
    bu.setup_logger(output_path, f'{cfg.prefix}.reconstruction', args)
    log.info('prefix=%s', cfg.prefix)
    log.info('output=%s', output_path)

    ############################################################################
    # Checks and configurations
    # - check parameters and setup global configuration
    # - test database
    # - test binary dependencies
    ############################################################################
    try:
        if args.input == '':
            raise ValueError('File path argument must be non-empty')
        annotation_path = Path(args.input).resolve()
        cfg.check_readability('annotation', annotation_path)
        cfg.check_content_size('annotation', annotation_path)
    except:
        log.error('provided annotation file not valid! path=%s', args.input)
        sys.exit(f'ERROR: annotation file ({args.input}) not valid!')
    log.info('input-path=%s', annotation_path)
    
    cfg.check_tmp_path(args)
    cfg.debug = args.debug
    log.info('debug=%s', cfg.debug)
    cfg.verbose = True if cfg.debug else args.verbose
    log.info('verbose=%s', cfg.verbose)
    
    if(cfg.verbose):
        print(f'Bakta v{cfg.version}')
        print('Options and arguments:')
        print(f'\tinput: {annotation_path}')
        print(f'\toutput: {cfg.output_path}')
        print(f'\tprefix: {cfg.prefix}')
        if(cfg.force): print(f'\tforce: {cfg.force}')
    
    if(cfg.debug):
        print(f"\nBakta runs in DEBUG mode! Temporary data will not be destroyed at: {cfg.tmp_path}")
    else:
        atexit.register(bu.cleanup, log, cfg.tmp_path)  # register cleanup exit hook
    
    ############################################################################
    # Import annotations from JSON
    ############################################################################
    print('Parse genome annotations...')
    with xopen(str(annotation_path), threads=0) as fh:
        data = json.load(fh)
    features = data['features']
    features_by_sequence = {seq['id']: [] for seq in data['sequences']}
    for feature in data['features']:
        seq_id = feature['sequence'] if 'sequence' in feature else feature['contig']  # <1.10.0 compatibility
        sequence_features = features_by_sequence.get(seq_id)
        sequence_features.append(feature)
    
    # set global config objects based on information from imported JSON document
    cfg.version = data['version']['bakta']
    cfg.db_info = {
        'type': data['version']['db']['type'],
        'major': data['version']['db']['version'].split('.')[0],
        'minor': data['version']['db']['version'].split('.')[1]
    }
    cfg.translation_table = data['genome']['translation_table']
    cfg.run_start = datetime.strptime(data['run']['start'], '%Y-%m-%d %H:%M:%S')
    cfg.run_end = datetime.strptime(data['run']['end'], '%Y-%m-%d %H:%M:%S')

    ############################################################################
    # Write output files
    # - write optional output files in GFF3/GenBank/EMBL formats
    # - measure runtime
    # - write comprehensive annotation results as JSON
    # - remove temp directory
    ############################################################################
    print(f'\nExport annotation results to: {cfg.output_path}')
    print('\thuman readable TSV...')
    tsv_path = cfg.output_path.joinpath(f'{cfg.prefix}.tsv')
    tsv.write_features(data['sequences'], features_by_sequence, tsv_path)

    print('\tGFF3...')
    gff3_path = cfg.output_path.joinpath(f'{cfg.prefix}.gff3')
    gff.write_features(data, features_by_sequence, gff3_path)

    print('\tINSDC GenBank & EMBL...')
    genbank_path = cfg.output_path.joinpath(f'{cfg.prefix}.gbff')
    embl_path = cfg.output_path.joinpath(f'{cfg.prefix}.embl')
    insdc.write_features(data, features, genbank_path, embl_path)

    print('\tgenome sequences...')
    fna_path = cfg.output_path.joinpath(f'{cfg.prefix}.fna')
    fasta.export_sequences(data['sequences'], fna_path, description=True, wrap=True)

    print('\tfeature nucleotide sequences...')
    ffn_path = cfg.output_path.joinpath(f'{cfg.prefix}.ffn')
    fasta.write_ffn(features, ffn_path)

    print('\ttranslated CDS sequences...')
    faa_path = cfg.output_path.joinpath(f'{cfg.prefix}.faa')
    fasta.write_faa(features, faa_path)

    print('\tfeature inferences...')
    tsv_path = cfg.output_path.joinpath(f'{cfg.prefix}.inference.tsv')
    tsv.write_feature_inferences(data['sequences'], features_by_sequence, tsv_path)

    print('\tcircular genome plot...')
    plot.write(data, features, cfg.output_path)

    hypotheticals = [feat for feat in features if feat['type'] == bc.FEATURE_CDS and 'hypothetical' in feat]
    print('\thypothetical TSV...')
    tsv_path = cfg.output_path.joinpath(f'{cfg.prefix}.hypotheticals.tsv')
    tsv.write_hypotheticals(hypotheticals, tsv_path)
    
    print('\ttranslated hypothetical CDS sequences...')
    faa_path = cfg.output_path.joinpath(f'{cfg.prefix}.hypotheticals.faa')
    fasta.write_faa(hypotheticals, faa_path)

    print('\tGenome and annotation summary...')
    summary_path = cfg.output_path.joinpath(f'{cfg.prefix}.txt')
    with summary_path.open('w') as fh_out:
        fh_out.write('Sequence(s):\n')
        fh_out.write(f"Length: {data['stats']['size']:}\n")
        fh_out.write(f"Count: {len(data['sequences'])}\n")
        fh_out.write(f"GC: {100 * data['stats']['gc']:.1f}\n")
        fh_out.write(f"N50: {data['stats']['n50']:}\n")
        if('n90' in data['stats']):
            fh_out.write(f"N90: {data['stats']['n90']:}\n")
        fh_out.write(f"N ratio: {100 * data['stats']['n_ratio']:.1f}\n")
        fh_out.write(f"coding density: {100 * data['stats']['coding_ratio']:.1f}\n")
        fh_out.write('\nAnnotation:\n')
        fh_out.write(f"tRNAs: {len([feat for feat in features if feat['type'] == bc.FEATURE_T_RNA])}\n")
        fh_out.write(f"tmRNAs: {len([feat for feat in features if feat['type'] == bc.FEATURE_TM_RNA])}\n")
        fh_out.write(f"rRNAs: {len([feat for feat in features if feat['type'] == bc.FEATURE_R_RNA])}\n")
        fh_out.write(f"ncRNAs: {len([feat for feat in features if feat['type'] == bc.FEATURE_NC_RNA])}\n")
        fh_out.write(f"ncRNA regions: {len([feat for feat in features if feat['type'] == bc.FEATURE_NC_RNA_REGION])}\n")
        fh_out.write(f"CRISPR arrays: {len([feat for feat in features if feat['type'] == bc.FEATURE_CRISPR])}\n")
        cdss = [feat for feat in features if feat['type'] == bc.FEATURE_CDS]
        fh_out.write(f"CDSs: {len(cdss)}\n")
        fh_out.write(f"pseudogenes: {len([cds for cds in cdss if 'pseudogene' in cds])}\n")
        fh_out.write(f"hypotheticals: {len([cds for cds in cdss if 'hypothetical' in cds])}\n")
        fh_out.write(f"sORFs: {len([feat for feat in features if feat['type'] == bc.FEATURE_SORF])}\n")
        fh_out.write(f"gaps: {len([feat for feat in features if feat['type'] == bc.FEATURE_GAP])}\n")
        fh_out.write(f"oriCs: {len([feat for feat in features if feat['type'] == bc.FEATURE_ORIC])}\n")
        fh_out.write(f"oriVs: {len([feat for feat in features if feat['type'] == bc.FEATURE_ORIV])}\n")
        fh_out.write(f"oriTs: {len([feat for feat in features if feat['type'] == bc.FEATURE_ORIT])}\n")
        fh_out.write('\nBakta:\n')
        fh_out.write(f'Software: v{cfg.version}\n')
        fh_out.write(f"Database: v{cfg.db_info['major']}.{cfg.db_info['minor']}, {cfg.db_info['type']}\n")
        fh_out.write('DOI: 10.1099/mgen.0.000685\n')
        fh_out.write('URL: github.com/oschwengers/bakta\n')


if __name__ == '__main__':
    main()
