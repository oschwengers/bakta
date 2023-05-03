import atexit
import json
import logging
import os
import subprocess as sp
import sys

from pathlib import Path

import yaml

from Bio import SeqUtils

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
    with annotation_path.open('r') as fh:
        annotation = json.load(fh)
    features = annotation['features']
    contigs = annotation['sequences']

    # load colors if specified
    colors = COLORS
    conf_colors = config.get('colors', None)
    if conf_colors is not None:
        colors = {**colors, **conf_colors}

    if args.sequences == 'all':  # write whole genome plot
        print(f'draw circular genome plot (type={plot_type}) containing all sequences...')
        write_plot(features, contigs, output_path, colors, plot_type=plot_type)
    else:  # write genome plot containing provided sequences only
        plot_contigs = []
        sequence_identifiers = []
        for selected_sequence in args.sequences.split(','):
            for i, contig in enumerate(contigs):
                sequence_no = str(i + 1)
                if selected_sequence == sequence_no:
                    plot_contigs.append(contig)
                    sequence_identifiers.append(sequence_no)
                elif selected_sequence.lower() == contig['id'].lower():
                    plot_contigs.append(contig)
                    sequence_identifiers.append(contig['id'])
        if len(plot_contigs) > 0:
            print(f'draw circular genome plot (type={plot_type}) containing sequences: {sequence_identifiers}...')
            plot_name_suffix = '_'.join(sequence_identifiers)
            plot_contig_ids = [c['id'] for c in plot_contigs]
            features = [feat for feat in features if feat['contig'] in plot_contig_ids]
            write_plot(features, plot_contigs, output_path, colors, plot_name_suffix=plot_name_suffix, plot_type=plot_type)


def write_plot(features, contigs, output_path, colors=COLORS, plot_name_suffix=None, plot_type=bc.PLOT_FEATURES):
    # config paths
    circos_path = cfg.tmp_path.joinpath(f'circos')
    circos_path.mkdir(parents=True, exist_ok=True)

    # fix edge features because Circos cannot handle them correctly
    non_edge_features = [feat for feat in features if not feat.get('edge', False)]
    contigs_by_id = {c['id']: c for c in contigs}
    for feat in [feat for feat in features if feat.get('edge', False)]:
        contig = contigs_by_id[feat['contig']]
        log.info('split edge feature: contig=%s, start=%i, stop=%i, strand=%s, edge=%s', contig['id'], feat['start'], feat['stop'], feat['strand'], feat['edge'])
        feat_1 = feat.copy()
        feat_1['stop'] = contig['length']
        feat_1['edge'] = False
        non_edge_features.append(feat_1)
        feat_2 = feat.copy()
        feat_2['start'] = 1
        feat_2['edge'] = False
        non_edge_features.append(feat_2)
    features = non_edge_features

    # write feature files
    if plot_type == bc.PLOT_COG:
        feature_paths = write_features_type_cog(features, contigs, circos_path, colors)
    else:
        feature_paths = write_features_type_feature(features, contigs, circos_path, colors)

    # write gc content and gc skew files
    tracks_path = circos_path.joinpath('tracks.conf')
    
    gc_content_path, max_gc, gc_skew_path, max_gc_skew = write_gc_content_skew(contigs, circos_path, colors)
    write_tracks(tracks_path, feature_paths, gc_content_path, max_gc, gc_skew_path, max_gc_skew)

    # write main config
    file_name = cfg.prefix if plot_name_suffix is None else f'{cfg.prefix}_{plot_name_suffix}'
    main_conf_path = write_main_config(circos_path, output_path, tracks_path, contigs, file_name, colors)
    
    # execute Circos
    log.info('write circular genome plot: file-name=%s, output-dir=%s', file_name, output_path)
    cmd = [
        'circos',
        '-conf',
        main_conf_path
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
        log.warning('Circos failed! circos-error-code=%d', proc.returncode)
        raise Exception(f'circos error! error code: {proc.returncode}')


def write_features_type_feature(features, contigs, circos_path, colors):
    features_plus = []
    features_minus = []
    contig_ids = set([c['id'] for c in contigs])
    for feat in features:
        if feat['contig'] not in contig_ids:
            continue
        contig, start, stop, type = feat['contig'], feat['start'], feat['stop'], feat['type']
        color = colors['features'].get(type, colors['features']['misc'])
        if feat['strand'] == bc.STRAND_FORWARD:
            features_plus.append(f"{contig} {start} {stop} {bc.STRAND_FORWARD} color={hex_to_rgb(color)}")
        else:
            features_minus.append(f"{contig} {start} {stop} {bc.STRAND_REVERSE} color={hex_to_rgb(color)}")
    features_plus_path = circos_path.joinpath('features-plus.txt')
    with features_plus_path.open('w') as fh:
        fh.write('\n'.join(features_plus))
        fh.write('\n')
    features_minus_path = circos_path.joinpath('features-minus.txt')
    with features_minus_path.open('w') as fh:
        fh.write('\n'.join(features_minus))
        fh.write('\n')
    return [features_plus_path, features_minus_path]


def write_features_type_cog(features, contigs, circos_path, colors):
    features_plus = []
    features_minus = []
    features_extra = []
    contig_ids = set([c['id'] for c in contigs])
    for feat in features:
        if feat['contig'] not in contig_ids:
            continue
        contig, start, stop = feat['contig'], feat['start'], feat['stop']
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
                features_plus.append(f"{contig} {start} {stop} {feat['strand']} color={hex_to_rgb(color)}")
            else:
                features_minus.append(f"{contig} {start} {stop} {feat['strand']} color={hex_to_rgb(color)}")
        else:
            features_extra.append(f"{contig} {start} {stop} {feat['strand']} color={hex_to_rgb(colors['features']['misc'])}")
    features_plus_path = circos_path.joinpath('features-plus.txt')
    with features_plus_path.open('w') as fh:
        fh.write('\n'.join(features_plus))
        fh.write('\n')
    features_minus_path = circos_path.joinpath('features-minus.txt')
    with features_minus_path.open('w') as fh:
        fh.write('\n'.join(features_minus))
        fh.write('\n')
    features_extra_path = circos_path.joinpath('features-extra.txt')
    with features_extra_path.open('w') as fh:
        fh.write('\n'.join(features_extra))
        fh.write('\n')
    return [features_plus_path, features_minus_path, features_extra_path]


def write_gc_content_skew(contigs, circos_path, colors):
    sequence_length = sum([c['length'] for c in contigs])
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
    gc_mean = SeqUtils.gc_fraction(''.join([c['sequence'] for c in contigs]))
    for contig in contigs:
        seq = contig['sequence']
        for w in range(0, len(seq), step_size):
            start = w - window_size
            if start < 0:
                start += len(seq)
            stop = w + window_size
            if stop > len(seq):
                stop -= len(seq)
            subseq = seq[start:stop] if start < stop else seq[start:] + seq[:stop]
            gc_value = gc_mean - SeqUtils.gc_fraction(subseq)
            if max_gc < abs(gc_value):
                max_gc = abs(gc_value)
            gc_color = colors['gc-positive'] if gc_value >= 0 else colors['gc-negative']
            gc_contents.append(f"{contig['id']} {w} {w} {gc_value} fill_color={hex_to_rgb(gc_color)}")
            g, c = subseq.count('G'), subseq.count('C')
            gc_skew = gc_skew = (g - c) / float(g + c) if (g + c) > 0 else 0.0
            if max_gc_skew < abs(gc_skew):
                max_gc_skew = abs(gc_skew)
            gc_skew_color = colors['gc-skew-positive'] if gc_skew >= 0 else colors['gc-skew-negative']
            gc_skews.append(f"{contig['id']} {w} {w} {gc_skew} fill_color={hex_to_rgb(gc_skew_color)}")

    log.debug('write gc config: seq-length=%i, step-size=%i, window-size=%i, max-gc=%i, max-gc-skew=%i', sequence_length, step_size, window_size, max_gc, max_gc_skew)
    gc_content_path = circos_path.joinpath('gc_content.txt')
    with gc_content_path.open('w') as fh:
        fh.write('\n'.join(gc_contents))
        fh.write('\n')
    gc_skew_path = circos_path.joinpath('gc_skew.txt')
    with gc_skew_path.open('w') as fh:
        fh.write('\n'.join(gc_skews))
        fh.write('\n')
    return gc_content_path, max_gc, gc_skew_path, max_gc_skew


def write_tracks(tracks_path, feature_paths, gc_content_path, max_gc, gc_skew_path, max_gc_skew):  # write track configuration
    track_radius = 1.0
    track_texts = []
    for feature in feature_paths:
        track_text = f'''
<plot>
type             = tile
file             = {feature}
r1               = {track_radius}r
r0               = {track_radius - 0.1}r
orientation      = out
layers           = 1
margin           = 0.01u
thickness        = 100
padding          = 1
stroke_color     = black
stroke_thickness = 0
layers_overflow  = collapse
</plot>
        '''
        track_texts.append(track_text)
        track_radius -= 0.1
    gc_radius = 0.2
    for gc_path, maximum in [(gc_content_path, max_gc), (gc_skew_path, max_gc_skew)]:
        minimum = -maximum
        track_text = f'''
<plot>
type        = histogram
file        = {gc_path}
r1          = {track_radius}r
r0          = {track_radius-gc_radius}r
min         = {minimum}
max         = {maximum}
thickness   = 0
orientation = out
</plot>    
        '''
        track_texts.append(track_text)
        track_radius -= gc_radius
    with tracks_path.open('w') as fh:
        fh.write('\n'.join(track_texts))
        fh.write('\n')
        

def hex_to_rgb(hex_string):
    hex_string = hex_string.replace('#', '')
    rgb = []
    for i in (0, 2, 4):
        rgb.append(str(int(hex_string[i:i+2], 16)))
    return ','.join(rgb)


def write_main_config(circos_path, output_path, tracks_path, contigs, file_name, colors):
    karyotype_path = circos_path.joinpath('karyotype.txt')
    sequence_length = sum([c['length'] for c in contigs])
    chromosomes_units = round(sequence_length/(10**(len(str(sequence_length)) - 1)))*(10**(len(str(sequence_length)) - 1))

    if sequence_length > 10_000:
        label_prefix = 'kbp'
        multiplier = 0.001
    else:
        label_prefix = 'bp'
        multiplier = 1
    
    log.debug('write main config: file-name=%s, output-dir=%s, multiplier=%s, label-prefix%s', file_name, output_path, multiplier, label_prefix)
    main_config_text = f'''
karyotype                   = {karyotype_path}
chromosomes_units           = {chromosomes_units}
chromosomes_display_default = yes
<plots>
<<include {tracks_path}>>
</plots>
<image>
<<include image.conf>>
file*                       = {file_name}
dir*                        = {output_path}
</image> 

<ideogram>
<spacing>
default     = 0.01r
</spacing>
radius      = 0.85r
thickness   = 5p
fill        = yes
</ideogram>

show_ticks       = yes
show_tick_labels = yes       
<ticks>
label_separation = 20p
tick_separation  = 2p 
radius           = 1r
color            = black
thickness        = 4p
label_offset     = 10p
multiplier       = {multiplier}
orientation      = out
format           = %d {label_prefix}
<tick>
spacing          = 0.1u
show_label       = yes
label_size       = 40   
size             = 25p
thickness        = 4p
</tick>
<tick>
spacing          = 0.025u
show_label       = yes
label_size       = 20
size             = 15p
thickness        = 3p
</tick>
<tick>
spacing          = 0.0025u
size             = 5p
</tick>
</ticks>

<<include etc/colors_fonts_patterns.conf>>
<<include etc/housekeeping.conf>>
max_ideograms*   = 10000
    '''
    main_conf_path = circos_path.joinpath('main.conf')
    with main_conf_path.open('w') as fh:
        fh.write(main_config_text)
        fh.write('\n')

    # write karyotype file
    karyotypes = []
    for i, c in enumerate(contigs):
        karyotypes.append(f"chr - {c['id']} {i + 1} 0 {c['length']} {hex_to_rgb(colors['backbone'])}")
    with karyotype_path.open('w') as fh:
        fh.write('\n'.join(karyotypes))
        fh.write('\n')

    return main_conf_path


if __name__ == '__main__':
    main()
