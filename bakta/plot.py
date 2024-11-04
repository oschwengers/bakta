import atexit
import json
import logging
import os
import subprocess as sp
import sys

from pathlib import Path

import yaml
import Bio as bp
import numpy as np

from Bio import SeqUtils
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from pycirclize import Circos
from pycirclize.utils import load_prokaryote_example_file

import bakta
import bakta.utils as bu
import bakta.config as cfg
import bakta.constants as bc


COLORS = {
    'backbone': '#000000',
    'gc-positive': '#33a02c',
    'gc-negative': '#e31a1c',
    'gc-skew-positive': '#fdbf6f',
    'gc-skew-negative': '#1f78b4',
    'features': {  # feature type colors
        bc.FEATURE_CDS: 'cccccc',
        bc.FEATURE_SORF: 'cccccc',
        bc.FEATURE_T_RNA: '#b2df8a',
        bc.FEATURE_TM_RNA: '#b2df8a',
        bc.FEATURE_R_RNA: '#fb8072',
        bc.FEATURE_NC_RNA: '#fdb462',
        bc.FEATURE_NC_RNA_REGION : '#80b1d3',
        bc.FEATURE_CRISPR: '#bebada',
        bc.FEATURE_GAP: '#000000',
        'misc': '#666666'
    },
    'cog-classes': {
        'A': '#c5000b',  # RNA processing and modification
        'B': '#6a3d9a',  # Chromatin structure and dynamics
        'C': '#bebada',  # Energy production and conversion
        'D': '#ff7f00',  # Cell cycle control, cell division, chromosome partitioning
        'E': '#ffffb3',  # Amino acid transport and metabolism
        'F': '#e31a1c',  # Nucleotide transport and metabolism
        'G': '#fb8072',  # Carbohydrate transport and metabolism
        'H': '#b2df8a',  # Coenzyme transport and metabolism
        'I': '#ccebc5',  # Lipid transport and metabolism
        'J': '#a6cee3',  # Translation, ribosomal structure and biogenesis
        'K': '#b3de69',  # Transcription
        'L': '#80b1d3',  # Replication, recombination and repair
        'M': '#bc80bd',  # Cell wall/membrane/envelope biogenesis
        'N': '#cab2d6',  # Cell motility
        'O': '#33a02c',  # Posttranslational modification, protein turnover, chaperones
        'P': '#fccde5',  # Inorganic ion transport and metabolism
        'Q': '#1f78b4',  # Secondary metabolites biosynthesis, transport and catabolism
        'R': '#8dd3c7',  # General function prediction only
        'S': '#222222',  # Function unknown
        'T': '#fdb462',  # Signal transduction mechanisms
        'U': '#fdbf6f',  # Intracellular trafficking, secretion, and vesicular transport
        'V': '#fb9a99',  # Defense mechanisms
        'W': '#0084d1',  # Extracellular structures
        'X': '#d9d9d9',  # Mobilome: prophages, transposons
        'Y': '#ffff38',  # Nuclear structure
        'Z': '#ffff99'   # Cytoskeleton
    }
}


log = logging.getLogger('PLOT')


def main():
    # parse options and arguments
    parser = bu.init_parser(sub_command='_plot')
    parser.add_argument('input', metavar='<input>', help='Bakta annotations in JSON format')
    
    arg_group_io = parser.add_argument_group('Input / Output')
    arg_group_io.add_argument('--config', '-c', action='store', default=None, help='Plotting configuration in YAML format')
    arg_group_io.add_argument('--output', '-o', action='store', default=os.getcwd(), help='Output directory (default = current working directory)')
    arg_group_io.add_argument('--prefix', '-p', action='store', default=None, help='Prefix for output files')
    arg_group_io.add_argument('--force', '-f', action='store_true', help='Force overwriting existing output folder')

    arg_group_plot = parser.add_argument_group('Plotting')
    arg_group_plot.add_argument('--sequences', action='store', default='all', help='Sequences to plot: comma separated number or name (default = all, numbers one-based)')
    arg_group_plot.add_argument('--type', action='store', type=str, default=bc.PLOT_FEATURES, choices=[bc.PLOT_FEATURES, bc.PLOT_COG], help=f'Plot type (default = {bc.PLOT_FEATURES})')

    arg_group_general = parser.add_argument_group('General')
    arg_group_general.add_argument('--help', '-h', action='help', help='Show this help message and exit')
    arg_group_general.add_argument('--verbose', '-v', action='store_true', help='Print verbose information')
    arg_group_general.add_argument('--debug', action='store_true', help='Run Bakta in debug mode. Temp data will not be removed.')
    arg_group_general.add_argument('--tmp-dir', action='store', default=None, dest='tmp_dir', help='Location for temporary files (default = system dependent auto detection)')
    arg_group_general.add_argument('--version', action='version', version=f'%(prog)s {bakta.__version__}')
    args = parser.parse_args()

    ############################################################################
    # Setup logging
    ############################################################################
    cfg.prefix = args.prefix if args.prefix else Path(args.input).stem
    output_path = cfg.check_output_path(args.output, args.force)
    cfg.force = args.force
    log.info('force=%s', args.force)
    
    bu.setup_logger(output_path, cfg.prefix, args)
    log.info('prefix=%s', cfg.prefix)
    log.info('output=%s', output_path)

    
    ############################################################################
    # Checks and configurations
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
    plot_type = args.type
    log.info('plot-type=%s', plot_type)

    # check and open configuration file
    if args.config is not None:
        config_path = None
        try:
            if(args.config == ''):
                raise ValueError('File path argument must be non-empty')
            config_path = Path(args.config).resolve()
            cfg.check_readability('config', config_path)
            cfg.check_content_size('config', config_path)
            log.info('config=%s', config_path)
            with config_path.open() as fh:
                config = yaml.load(fh, Loader=yaml.FullLoader)
        except:
            log.error('provided config file not valid! path=%s', args.config)
            sys.exit(f'ERROR: config file ({args.config}) not valid!')
    else:
        config = {}

    bu.test_dependencies()
    if(cfg.verbose):
        print(f'Bakta v{bakta.__version__}')
        print('Options and arguments:')
        print(f'\tinput: {annotation_path}')
        if(args.config): print(f'\tconfig: {args.config}')
        print(f'\toutput: {cfg.output_path}')
        if(cfg.force): print(f'\tforce: {cfg.force}')
        print(f'\ttmp directory: {cfg.tmp_path}')
        print(f'\tprefix: {cfg.prefix}')

    if(cfg.debug):
        print(f"\nBakta runs in DEBUG mode! Temporary data will not be destroyed at: {cfg.tmp_path}")
    else:
        atexit.register(bu.cleanup, log, cfg.tmp_path)  # register cleanup exit hook

    ############################################################################
    # Write plots
    ############################################################################

    # load genome annotations
    print('Parse genome annotations...')
    with annotation_path.open('r') as fh:
        annotation = json.load(fh)
    features = annotation['features']
    sequences = annotation['sequences']

    # load colors if specified
    colors = COLORS
    conf_colors = config.get('colors', None)
    if conf_colors is not None:
        colors = {**colors, **conf_colors}

    print('Draw plots...')
    if args.sequences == 'all':  # write whole genome plot
        print(f'\tdraw circular genome plot (type={plot_type}) containing all sequences...')
        write(features, sequences, output_path, colors, plot_type=plot_type)
    else:  # write genome plot containing provided sequences only
        plot_sequences = []
        sequence_identifiers = []
        for selected_sequence in args.sequences.split(','):
            for i, seq in enumerate(sequences):
                sequence_no = str(i + 1)
                if selected_sequence == sequence_no:
                    plot_sequences.append(seq)
                    sequence_identifiers.append(sequence_no)
                elif selected_sequence.lower() == seq['id'].lower():
                    plot_sequences.append(seq)
                    sequence_identifiers.append(seq['id'])
        if len(plot_sequences) > 0:
            print(f'\tdraw circular genome plot (type={plot_type}) containing sequences: {sequence_identifiers}...')
            plot_name_suffix = '_'.join(sequence_identifiers)
            plot_sequence_ids = [seq['id'] for seq in plot_sequences]
            features = [feat for feat in features if feat['sequence'] in plot_sequence_ids]
            write(features, plot_sequences, output_path, colors, plot_name_suffix=plot_name_suffix, plot_type=plot_type)


def write(features, sequences, output_path, colors=COLORS, plot_name_suffix=None, plot_type=bc.PLOT_FEATURES):
    # config paths
    circos_path = cfg.tmp_path.joinpath(f'circos')
    circos_path.mkdir(parents=True, exist_ok=True)

    # select style
    if plot_type == bc.PLOT_COG:
        plot = write_features_type_cog(features, sequences, circos_path, colors)
    else:
        plot = write_features_type_feature(features, sequences, circos_path, colors)

    max_gc, max_gc_skew = write_gc_content_skew(sequences, circos_path, colors)
    

    


def write_features_type_feature(features, sequences, circos_path, colors):
    features_plus = []
    features_minus = []
    sequence_ids = set([seq['id'] for seq in sequences])
    for feat in features:
        if feat['sequence'] not in sequence_ids:
            continue
        seq, start, stop, type = feat['sequence'], feat['start'], feat['stop'], feat['type']
        color = colors['features'].get(type, colors['features']['misc'])
        if feat['strand'] == bc.STRAND_FORWARD:
            features_plus.append(f"{seq} {start} {stop} {bc.STRAND_FORWARD} color={hex_to_rgb(color)}")
        else:
            features_minus.append(f"{seq} {start} {stop} {bc.STRAND_REVERSE} color={hex_to_rgb(color)}")
    return [features_plus_path, features_minus_path]


def write_features_type_cog(features, sequences, circos_path, colors):
    features_plus = []
    features_minus = []
    features_extra = []
    sequence_ids = set([seq['id'] for seq in sequences])
    for feat in features:
        if feat['sequence'] not in sequence_ids:
            continue
        seq, start, stop = feat['sequence'], feat['start'], feat['stop']
        if feat['type'] == bc.FEATURE_CDS:
            color = colors['features'][bc.FEATURE_CDS]
            psc = feat.get('psc', None)
            if psc is not None:
                cog = psc.get('cog_category', None)
                if cog is not None:
                    if len(cog) != 1:
                        cog = cog[:1]
                    color = colors['cog-classes'].get(cog.upper(), colors['cog-classes']['S'])
            if feat['strand'] == bc.STRAND_FORWARD:
                features_plus.append(f"{seq} {start} {stop} {feat['strand']} color={hex_to_rgb(color)}")
            else:
                features_minus.append(f"{seq} {start} {stop} {feat['strand']} color={hex_to_rgb(color)}")
        else:
            features_extra.append(f"{seq} {start} {stop} {feat['strand']} color={hex_to_rgb(colors['features']['misc'])}")


def write_gc_content_skew(sequences, circos_path, colors):
    sequence_length = sum([seq['length'] for seq in sequences])
    step_size = int(sequence_length / 3600)  # 10 * 360°
    if step_size < 3:
        step_size = 3
    window_size = 2 * step_size
    if window_size < 50:
        window_size = 50
    gc_contents = []
    gc_skews = []
    max_gc = 0
    max_gc_skew = 0
    if float(bp.__version__) >= 1.80:
        gc_mean = SeqUtils.gc_fraction(''.join([seq['nt'] for seq in sequences]))
    else:
        gc_mean = SeqUtils.GC(''.join([seq['nt'] for seq in sequences])) / 100
    for seq in sequences:
        nt = seq['nt']
        for w in range(0, len(nt), step_size):
            start = w - window_size
            if start < 0:
                start += len(nt)
            stop = w + window_size
            if stop > len(nt):
                stop -= len(nt)
            nt_subseq = nt[start:stop] if start < stop else nt[start:] + nt[:stop]
            if float(bp.__version__) >= 1.80:
                gc_value = gc_mean - SeqUtils.gc_fraction(nt_subseq)
            else:
                gc_value = gc_mean - (SeqUtils.GC(nt_subseq) / 100)
            if max_gc < abs(gc_value):
                max_gc = abs(gc_value)
            g, c = nt_subseq.count('G'), nt_subseq.count('C')
            gc_skew = gc_skew = (g - c) / float(g + c) if (g + c) > 0 else 0.0
            if max_gc_skew < abs(gc_skew):
                max_gc_skew = abs(gc_skew)

    log.debug('write gc config: seq-length=%i, step-size=%i, window-size=%i, max-gc=%i, max-gc-skew=%i', sequence_length, step_size, window_size, max_gc, max_gc_skew)
    return max_gc, max_gc_skew



if __name__ == '__main__':
    main()
