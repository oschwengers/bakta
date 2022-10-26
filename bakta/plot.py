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


NT_FEATURE_COLOR = '#888888'

COG_COLORS = {
    'J': '#ff0000',  # Translation, ribosomal structure and biogenesis
    'A': '#c2af58',  # RNA processing and modification
    'K': '#ff9900',  # Transcription
    'L': '#ffff00',  # Replication, recombination and repair
    'B': '#ffc600',  # Chromatin structure and dynamics
    'D': '#99ff00',  # Cell cycle control, cell division, chromosome partitioning
    'Y': '#493126',  # Nuclear structure
    'V': '#ff008a',  # Defense mechanisms
    'T': '#0000ff',  # Signal transduction mechanisms
    'M': '#9ec928',  # Cell wall/membrane/envelope biogenesis
    'N': '#006633',  # Cell motility
    'Z': '#660099',  # Cytoskeleton
    'W': '#336699',  # Extracellular structures
    'U': '#33cc99',  # Intracellular trafficking, secretion, and vesicular transport
    'O': '#00ffff',  # Posttranslational modification, protein turnover, chaperones
    'C': '#9900ff',  # Energy production and conversion
    'G': '#805642',  # Carbohydrate transport and metabolism
    'E': '#ff00ff',  # Amino acid transport and metabolism
    'F': '#99334d',  # Nucleotide transport and metabolism
    'H': '#727dcc',  # Coenzyme transport and metabolism
    'I': '#5c5a1b',  # Lipid transport and metabolism
    'P': '#0099ff',  # Inorganic ion transport and metabolism
    'Q': '#ffcc99',  # Secondary metabolites biosynthesis, transport and catabolism
    'R': '#ff9999',  # General function prediction only
    'S': '#d6aadf'	 # Function unknown
    }
COG_DEFAULT_COLOR = '#000000'
POSITIVE_GC_COLOR = '#CC6458'
NEGATIVE_GC_COLOR = '#43CC85'
POSITIVE_GC_SKEW_COLOR = '#CCBE6C'
NEGATIVE_GC_SKEW_COLOR = '#5A4ECC'
BACKBONE_COLOR = '#000000'


log = logging.getLogger('PLOT')


def main():
    # parse options and arguments
    parser = bu.init_parser(sub_command='_plot')

    parser.add_argument('input', metavar='<input>', help='Bakta annotations in JSON format')
    parser.add_argument('--config', '-c', action='store', default=None, help='Plotting configuration in YAML format')
    
    arg_group_io = parser.add_argument_group('Input / Output')
    arg_group_io.add_argument('--output', '-o', action='store', default=os.getcwd(), help='Output directory (default = current working directory)')
    arg_group_io.add_argument('--prefix', '-p', action='store', default='Plot', help='Prefix for output files')

    arg_group_plot = parser.add_argument_group('Plotting')
    arg_group_plot.add_argument('--plots', action='store', default='all', nargs='+', help='Sequences to plot: no or name (default = all))')

    arg_group_general = parser.add_argument_group('Runtime & auxiliary options')
    arg_group_general.add_argument('--help', '-h', action='help', help='Show this help message and exit')
    arg_group_general.add_argument('--verbose', '-v', action='store_true', help='Print verbose information')
    arg_group_general.add_argument('--debug', action='store_true', help='Run Bakta in debug mode. Temp data will not be removed.')
    arg_group_general.add_argument('--tmp-dir', action='store', default=None, dest='tmp_dir', help='Location for temporary files (default = system dependent auto detection)')
    arg_group_general.add_argument('--version', '-V', action='version', version=f'%(prog)s {bakta.__version__}')
    args = parser.parse_args()

    ############################################################################
    # Setup logging
    ############################################################################
    cfg.prefix = args.prefix if args.prefix else Path(args.genome).stem
    output_path = cfg.check_output_path(args)
    
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
    bu.test_dependencies()

    if(cfg.debug):
        print(f"\nBakta runs in DEBUG mode! Temporary data will not be destroyed at: {cfg.tmp_path}")
    else:
        atexit.register(bu.cleanup, log, cfg.tmp_path)  # register cleanup exit hook

    with annotation_path.open('r') as fh:
        annotation = json.load(fh)
    features = annotation['features']
    contigs = annotation['sequences']

    if args.config is not None:  # check and open configuration file
        try:
            with open(args.config) as conf:
                config = yaml.load(conf, Loader=yaml.FullLoader)
        except FileNotFoundError:
            log.error(f'configuration file not found!, {args.config}')
            sys.exit(f'ERROR: please check your selected configuration file!, {args.config}')
    else:
        config = {}

    # load colors if specified
    positive_gc_color = config.get('pgcc', POSITIVE_GC_COLOR)
    negative_gc_color = config.get('ngcc', NEGATIVE_GC_COLOR)
    positive_gc_skew_color = config.get('pgcs', POSITIVE_GC_SKEW_COLOR)
    negative_gc_skew_color = config.get('ngcs', NEGATIVE_GC_SKEW_COLOR)

    if args.plots == 'all':  # 
        plot_name = '_all_sequences'
        print('draw plot...')
        write_plot(features, contigs, output_path, plot_name, positive_gc_color, negative_gc_color, positive_gc_skew_color, negative_gc_skew_color)
    else:  # write plot according to plot configuration
        for p in args.plots:
            contig_selection = False
            if 'all' in p:
                plot_count = 1
                for contig in contigs:
                    plot_contig = [contig]
                    plot_nr = f'_{plot_count}'
                    print(f'draw plot{plot_nr}...')
                    write_plot(features, plot_contig, output_path, plot_nr, positive_gc_color, negative_gc_color, positive_gc_skew_color, negative_gc_skew_color)
                    plot_count += 1

            p = p.split('/')
            plot_contig = []
            for contig in contigs:
                contig_nr = contig['id'][7:]
                if contig_nr in p:
                    plot_contig.append(contig)
                    contig_selection = True
            if contig_selection is True:
                plot_nr = ''
                for x in p:
                    plot_nr += f'_{x}'
                print(f'draw plot{plot_nr}...')
                write_plot(features, plot_contig, output_path, plot_nr, positive_gc_color, negative_gc_color, positive_gc_skew_color, negative_gc_skew_color)


def write_plot(features, contigs, output_path, plot_nr='', positive_gc_color=POSITIVE_GC_COLOR, negative_gc_color=NEGATIVE_GC_COLOR, positive_gc_skew_color=POSITIVE_GC_SKEW_COLOR, negative_gc_skew_color=NEGATIVE_GC_SKEW_COLOR):
    sequence_length = sum([c['length'] for c in contigs])
    sequences = ''.join([c['sequence'] for c in contigs])
    window_size = int(sequence_length/100) if sequence_length < 10000 else int(sequence_length/1000)
    window_size = 3 if window_size < 3 else window_size
    step_size = int(window_size * 0.2) if window_size >= 5 else 1
    track_radius = 1.0
    gc_radius = 0.2
    if sequence_length > 10_000_000:
        label_prefix = 'mbp'
        multiplier = 0.000001
    elif sequence_length > 100_000:
        label_prefix = 'kbp'
        multiplier = 0.001
    else:
        label_prefix = 'bp'
        multiplier = 1
    
    # config paths
    circos_path = cfg.tmp_path.joinpath(f'circos')
    circos_path.mkdir(parents=True, exist_ok=True)

    # write feature files
    cds_plus_features = []
    cds_minus_features = []
    noncoding_features = []
    for feat in features:
        contig, start, stop, color = feat['contig'], feat['start'], feat['stop'], NT_FEATURE_COLOR
        if feat['type'] == 'cds':
            psc = feat.get('psc', None)
            if psc is None:
                continue
            cog = psc.get('cog_category', None)
            if cog is None:
                continue
            if len(cog) != 1:
                cog = cog[:1]
            color = COG_COLORS.get(cog, COG_DEFAULT_COLOR)
            
            if feat['strand'] == bc.STRAND_FORWARD:
                cds_plus_features.append(f"{contig} {start} {stop} {feat['strand']} color={hex_to_rgb(color)}")
            else:
                cds_minus_features.append(f"{contig} {start} {stop} {feat['strand']} color={hex_to_rgb(color)}")
        else:
            noncoding_features.append(f"{contig} {start} {stop} {feat['strand']} color={hex_to_rgb(color)}")
    cds_plus_path = circos_path.joinpath('cds-plus.txt')
    with cds_plus_path.open('w') as fh:
        fh.write('\n'.join(cds_plus_features))
    cds_minus_path = circos_path.joinpath('cds-minus.txt')
    with cds_minus_path.open('w') as fh:
        fh.write('\n'.join(cds_minus_features))
    nc_path = circos_path.joinpath('nc.txt')
    with nc_path.open('w') as fh:
        fh.write('\n'.join(noncoding_features))

    # write gc content unf gc skew files
    gc_contents = []
    gc_skews = []
    max_gc_content = 0
    max_gc_skew = 0
    gc_mean = SeqUtils.GC(sequences)
    for contig in contigs:
        seq = contig['sequence']
        for w in range(0, len(seq), step_size):
            start = int(w - (window_size / 2))
            if start < 0:
                start += len(seq)
            stop = int(w + (window_size / 2))
            if stop > len(seq):
                stop -= len(seq)
            if start < stop:
                subseq = seq[start:stop]
            else:
                subseq = seq[start:] + seq[:stop]
            gc_value = gc_mean - SeqUtils.GC(subseq)
            if max_gc_content < abs(gc_value):
                max_gc_content = abs(gc_value)
            content_color = positive_gc_color if gc_value > 0 else negative_gc_color
            gc_contents.append(f"{contig['id']} {w} {w} {gc_value} fill_color={hex_to_rgb(content_color)}")
            g = subseq.count('G')
            c = subseq.count('C')
            if (g + c) > 0:
                gc_skew = (g - c) / float(g + c)
            else:
                gc_skew = 0.0
            if max_gc_skew < abs(gc_skew):
                max_gc_skew = abs(gc_skew)
            skew_color = positive_gc_skew_color if gc_skew > 0 else negative_gc_skew_color
            gc_skews.append(f"{contig['id']} {w} {w} {gc_skew} fill_color={hex_to_rgb(skew_color)}")

    gc_content_path = circos_path.joinpath('gc_content.txt')
    with gc_content_path.open('w') as fh:
        fh.write('\n'.join(gc_contents))
    gc_skew_path = circos_path.joinpath('gc_skew.txt')
    with gc_skew_path.open('w') as fh:
        fh.write('\n'.join(gc_skews))

    # write main config
    karyotype_path = circos_path.joinpath('karyotype.txt')
    ideogram_path = circos_path.joinpath('ideogram.conf')
    ticks_path = circos_path.joinpath('ticks.conf')
    tracks_path = circos_path.joinpath('tracks.conf')
    chromosomes_units = round(sequence_length/(10**(len(str(sequence_length)) - 1)))*(10**(len(str(sequence_length)) - 1))
    main_config_text = f'''
    karyotype                   = {karyotype_path}
    chromosomes_units           = {chromosomes_units}
    chromosomes_display_default = yes
    <plots>
    <<include {tracks_path}>>
    </plots>
    <image>
    <<include image.conf>>
    file*                       = {cfg.prefix}{plot_nr}
    dir*                        = {output_path}
    </image> 
    <<include {ideogram_path}>>
    <<include {ticks_path}>>
    <<include etc/colors_fonts_patterns.conf>>
    <<include etc/housekeeping.conf>>
    '''
    main_conf = circos_path.joinpath('main.conf')
    with main_conf.open('w') as fh:
        fh.write(main_config_text)

    # write karyotype file
    karyotypes = []
    for i, c in enumerate(contigs):
        karyotypes.append(f"chr - {c['id']} {i + 1} 0 {c['length']} {hex_to_rgb(BACKBONE_COLOR)}")
    with karyotype_path.open('w') as fh:
        fh.write('\n'.join(karyotypes))

    # write ideogram config
    ideogram_text = '''
    <ideogram>
    <spacing>
    default = 0.01r
    </spacing>
    radius           = 0.85r
    thickness        = 5p
    fill             = yes
    </ideogram>
    '''
    with ideogram_path.open('w') as fh:
        fh.write(ideogram_text)

    # write ticks config:
    ticks_text = f'''
    show_ticks       = yes
    show_tick_labels = yes       
    <ticks>
    label_separation = 20p
    tick_separation      = 2p 
    radius      = 1r
    color       = black
    thickness    = 4p
    label_offset = 10p
    multiplier  = {multiplier}
    orientation = out
    format      = %d{label_prefix}
    <tick>
    spacing      = 0.1u
    show_label   = yes
    label_size   = 40   
    size         = 25p
    thickness    = 4p
    </tick>
    <tick>
    spacing      = 0.025u
    show_label   = yes
    label_size   = 20
    size         = 15p
    thickness    = 3p
    </tick>
    <tick>
    spacing      = 0.0025u
    size         = 5p
    </tick>
    </ticks>
    '''
    with ticks_path.open('w') as fh:
        fh.write(ticks_text)

    # write track configuration
    track_texts = []
    for feature in [cds_plus_path, cds_minus_path, nc_path]:
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
    for gc_path in [gc_content_path, gc_skew_path]:
        maximum = max_gc_content if gc_path == gc_content_path else max_gc_skew
        minimum = -maximum
        track_text = f'''
        <plot>
        type = histogram
        file = {gc_path}
        r1 = {track_radius}r
        r0 = {track_radius-gc_radius}r
        min = {minimum}
        max = {maximum}
        thickness = 0
        orientation = out
        </plot>    
        '''
        track_texts.append(track_text)
        track_radius -= gc_radius
    with tracks_path.open('w') as fh:
        fh.write('\n'.join(track_texts))
    
    # run Circos
    cmd = [
        'circos',
        '-conf',
        main_conf
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
        log.debug('stdout=\'%s\', stderr=\'%s\'', proc.stdout, proc.stderr)
        log.warning('Circos failed! circos-error-code=%d', proc.returncode)
        raise Exception(f'circos error! error code: {proc.returncode}')


    # if Path(f'{output_path}/{cfg.prefix}{plot_nr}.png').is_file():
    #     pass
    # else:
    #     print('ERROR: Circos failed to draw your image! Please check your installation!')
    # return


def hex_to_rgb(hex_string):
    hex_string = hex_string.replace('#', '')
    rgb = []
    for i in (0, 2, 4):
        rgb.append(str(int(hex_string[i:i+2], 16)))
    return ','.join(rgb)


def khp(hex_with_prefix):
    hex_without_prefix = hex_with_prefix[2:]
    return hex_without_prefix


if __name__ == '__main__':
    main()
