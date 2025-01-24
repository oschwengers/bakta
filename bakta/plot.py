import atexit
import copy
import json
import logging
import os
import sys

from datetime import datetime
from pathlib import Path

import numpy as np
import yaml

from Bio.SeqFeature import FeatureLocation, CompoundLocation, AfterPosition, BeforePosition
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from pycirclize import Circos
from xopen import xopen

import bakta
import bakta.utils as bu
import bakta.config as cfg
import bakta.constants as bc
import bakta.io.insdc as insdc


COLORS = {
    'backbone': '#000000',
    'gc-positive': '#33a02c',
    'gc-negative': '#e31a1c',
    'gc-skew-positive': '#fdbf6f',
    'gc-skew-negative': '#1f78b4',
    'features': {  # feature type colors
        bc.FEATURE_CDS: '#cccccc',
        bc.FEATURE_SORF: '#cccccc',
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
    arg_group_plot.add_argument('--label', action='store', type=str, default=None, help=f"Plot center label (for line breaks use '|')")
    arg_group_plot.add_argument('--size', action='store', type=int, default=8, choices=[4, 8, 16], help='Plot size in inches: 4/8/16 (default = 8)')
    arg_group_plot.add_argument('--dpi', action='store', type=int, default=300, choices=[150, 300, 600], help='Plot resolution as dots per inch: 150/300/600 (default = 300)')

    arg_group_general = parser.add_argument_group('General')
    arg_group_general.add_argument('--help', '-h', action='help', help='Show this help message and exit')
    arg_group_general.add_argument('--verbose', '-v', action='store_true', help='Print verbose information')
    arg_group_general.add_argument('--debug', action='store_true', help='Run Bakta in debug mode. Temp data will not be removed.')
    arg_group_general.add_argument('--tmp-dir', action='store', default=None, dest='tmp_dir', help='Location for temporary files (default = system dependent auto detection)')
    arg_group_general.add_argument('--version', action='version', version=f'%(prog)s {cfg.version}')
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
    plot_label = args.label
    log.info('plot-label=%s', plot_label)
    plot_size = args.size
    log.info('plot-size=%s', plot_size)
    plot_dpi = args.dpi
    log.info('plot-dpi=%s', plot_dpi)

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
        print(f'Bakta v{cfg.version}')
        print('Options and arguments:')
        print(f'\tinput: {annotation_path}')
        if(args.config): print(f'\tconfig: {args.config}')
        print(f'\toutput: {cfg.output_path}')
        if(cfg.force): print(f'\tforce: {cfg.force}')
        print(f'\ttmp directory: {cfg.tmp_path}')
        print(f'\tprefix: {cfg.prefix}')
        print(f'\tlabel: {plot_label}')
        print(f'\tsize: {plot_size}')
        print(f'\tDPI: {plot_dpi}')

    if(cfg.debug):
        print(f"\nBakta runs in DEBUG mode! Temporary data will not be destroyed at: {cfg.tmp_path}")
    else:
        atexit.register(bu.cleanup, log, cfg.tmp_path)  # register cleanup exit hook

    ############################################################################
    # Write plots
    ############################################################################

    # load genome annotations
    print('Parse genome annotations...')
    with xopen(str(annotation_path), threads=0) as fh:
        data = json.load(fh)
    
    # set global config objects based on information from imported JSON document
    cfg.db_info = {
        'type': data['version']['db']['type'],
        'major': data['version']['db']['version'].split('.')[0],
        'minor': data['version']['db']['version'].split('.')[1]
    }
    cfg.translation_table = data['genome']['translation_table']
    cfg.run_start = datetime.strptime(data['run']['start'], '%Y-%m-%d %H:%M:%S')
    cfg.run_end = datetime.strptime(data['run']['end'], '%Y-%m-%d %H:%M:%S')

    features = data['features']
    sequences = data['sequences']

    # load colors if specified
    colors = COLORS
    conf_colors = config.get('colors', None)
    if conf_colors is not None:
        colors = {**colors, **conf_colors}

    print('Draw plots...')
    if args.sequences == 'all':  # write whole genome plot
        print(f'\tdraw circular genome plot (type={plot_type}) containing all sequences...')
        write(data, features, output_path, colors, plot_type=plot_type, plot_label=plot_label, plot_size=plot_size, plot_dpi=plot_dpi)
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
            data['features'] = [feat for feat in features if feat['sequence'] in plot_sequence_ids]  # reduce feature list in data object
            data['sequences'] = [seq for seq in sequences if seq['id'] in plot_sequence_ids]  # reduce sequence list in data object
            write(data, features, output_path, colors, plot_name_suffix=plot_name_suffix, plot_type=plot_type, plot_label=plot_label, plot_size=plot_size, plot_dpi=plot_dpi)


def write(data, features, output_path, colors=COLORS, plot_name_suffix=None, plot_type=bc.PLOT_FEATURES, plot_label=None, plot_size=8, plot_dpi=300):
    sequence_list = insdc.build_biopython_sequence_list(data, features)
    clipped = False
    if(len(sequence_list) > 20):
        sequence_list = sorted(sequence_list, key=lambda a: len(a.seq), reverse=True)[:20]  # select longest 20 sequences in draft mode
        clipped = True
    for seq in sequence_list:  # fix edge features because PyCirclize cannot handle them correctly
        seq.features = [feat for feat in seq.features if feat.type != 'gene' and feat.type != 'source']
        for feat in seq.features:
            if isinstance(feat.location, CompoundLocation):
                feat_loc = feat.location
                log.debug('split edge feature: seq=%s, start=%i, stop=%i, strand=%s', seq.id, feat_loc.start, feat_loc.end, '+' if feat_loc.strand==1 else '-')
                if(feat_loc.strand == +1):
                    feat.location = FeatureLocation(feat_loc.parts[0].start, len(seq.seq), strand=feat_loc.strand)
                    feat_2 = copy.deepcopy(feat)
                    feat_2.location = FeatureLocation(0, feat_loc.parts[1].end, strand=feat_loc.strand)
                    seq.features.append(feat_2)
                elif(feat_loc.strand == -1):
                    feat.location = FeatureLocation(0, feat_loc.parts[0].end, strand=feat_loc.strand)
                    feat_2 = copy.deepcopy(feat)
                    feat_2.location = FeatureLocation(feat_loc.parts[1].start,  len(seq.seq), strand=feat_loc.strand)
                    seq.features.append(feat_2)
            else:
                if isinstance(feat.location.start, AfterPosition) or isinstance(feat.location.start, BeforePosition):
                    feat.location = FeatureLocation(int(str(feat.location.start)[1:]), feat.location.end, strand=feat.location.strand)
                if isinstance(feat.location.end, AfterPosition) or isinstance(feat.location.end, BeforePosition):
                    feat.location = FeatureLocation(feat.location.start, int(str(feat.location.end)[1:]), strand=feat.location.strand)

    # build lable
    plot_label = build_label(data) if plot_label is None else plot_label.replace('|', '\n')

    # select style
    if plot_type == bc.PLOT_COG:
        plot = build_features_type_cog(data, sequence_list, clipped, colors, plot_label, plot_size, plot_dpi)
    else:
        plot = build_features_type_feature(data, sequence_list, clipped, colors, plot_label, plot_size, plot_dpi)
    file_name = cfg.prefix if plot_name_suffix is None else f'{cfg.prefix}_{plot_name_suffix}'
    for file_type in ['png', 'svg']:
        file_path = output_path.joinpath(f'{file_name}.{file_type}')
        plot.savefig(file_path)


def build_label(data):
    genus = data['genome'].get('genus', '') if data['genome'].get('genus', None) else ''
    species = data['genome'].get('species', '') if data['genome'].get('species', None) else ''
    taxon = ' '.join([genus, species]).replace('  ', ' ')
    if taxon == '' or taxon == ' ':
        taxon = None
    strain = data['genome'].get('strain', '') if data['genome'].get('strain', None) else ''
    genome_size = sum([len(seq['nt']) for seq in data['sequences']]) if 'nt' in data['sequences'][0] else sum([len(seq['sequence']) for seq in data['sequences']])  # <1.10.0 compatibility
    genome_size_lable = f'{(genome_size/(10**6)):0.1f} Mbp' if genome_size >= 10**6 else f'{(genome_size/(10**3)):0.1f} kbp'
    label_list = []
    sequence_number = len(data['sequences'])
    if sequence_number == 1:
        seq_type = bc.REPLICON_CONTIG
        for seq in data['sequences']:
            if seq['id'] == data['sequences'][0]['id']:
                seq_type = seq['type']
        if seq_type == bc.REPLICON_CHROMOSOME:
            label_list.append(taxon)
            label_list.append(strain)
            label_list.append(f'chromosome, {genome_size_lable}')
        elif seq_type == bc.REPLICON_PLASMID:
            plasmid_name = data['genome'].get('plasmid', None)
            label_list.append(plasmid_name if plasmid_name else seq['id'])
            label_list.append(f'plasmid, {genome_size_lable}')
        else:  # a single contig
            label_list.append(taxon)
            label_list.append(strain)
            label_list.append(f"{seq['id']}, {genome_size_lable}")
    else:  # full (draft) genome
        label_list.append(taxon)
        label_list.append(strain)
        label_list.append(f'{sequence_number} sequences, {genome_size_lable}')
    return '\n'.join([lable for lable in label_list if lable is not None])


def build_features_type_feature(data, sequence_list, clipped, colors, plot_label, plot_size, plot_dpi):
    # Get contig genome seqid & size, features dict
    total_sequence_length = sum([len(seq['nt']) for seq in data['sequences']]) if 'nt' in data['sequences'][0] else sum([len(seq['sequence']) for seq in data['sequences']])  # <1.10.0 compatibility
    seqid2seq = {rec.id:rec.seq for rec in sequence_list}
    seqid2size = {rec.id:len(rec.seq) for rec in sequence_list}
    seqid2features = {rec.id:rec.features for rec in sequence_list}

    if plot_size == 4:
        text_size = 6
    elif plot_size == 8:
        text_size = 12
    else:
        text_size = 30
    
    if(clipped):
        circos = Circos(seqid2size, space=2, start=10, end=350)
        circos.text('...', r=99.5, color=colors['backbone'], ha='center', va='center', size=10)
    else:
        circos = Circos(seqid2size, space=2)
    circos.text(plot_label, r=7, size=text_size, linespacing=1.5)
    for sector in circos.sectors:
        # build tracks
        outer_track = sector.add_track((99.5, 100))
        feature_forward_track = sector.add_track((87, 97), r_pad_ratio=0.1)
        feature_reverse_track = sector.add_track((77, 87), r_pad_ratio=0.1)
        gc_content_track = sector.add_track((62, 72))
        gc_skew_track = sector.add_track((50, 60))

        # plot outer track
        build_sequence_backbone_track(sector, outer_track, total_sequence_length, colors, plot_size)

        # plot feature tracks
        for feature in seqid2features[sector.name]:
            track = feature_forward_track if feature.location.strand == 1 else feature_reverse_track
            if feature.type == bc.INSDC_FEATURE_CDS:
                track.genomic_features([feature], fc=colors['features'][bc.FEATURE_CDS])
            elif feature.type == bc.INSDC_FEATURE_T_RNA:
                track.genomic_features([feature], fc=colors['features'][bc.FEATURE_T_RNA])
            elif feature.type == bc.INSDC_FEATURE_TM_RNA:
                track.genomic_features([feature], fc=colors['features'][bc.FEATURE_TM_RNA])
            elif feature.type == bc.INSDC_FEATURE_R_RNA:
                track.genomic_features([feature], fc=colors['features'][bc.FEATURE_R_RNA])
            elif feature.type == bc.INSDC_FEATURE_NC_RNA:
                track.genomic_features([feature], fc=colors['features'][bc.FEATURE_NC_RNA])
            elif feature.type == bc.INSDC_FEATURE_REGULATORY:
                track.genomic_features([feature], fc=colors['features'][bc.FEATURE_NC_RNA_REGION])
            elif feature.type == bc.INSDC_FEATURE_REPEAT_REGION:
                track.genomic_features([feature], fc=colors['features'][bc.FEATURE_CRISPR])
            elif feature.type == bc.INSDC_FEATURE_GAP:
                track.genomic_features([feature], fc=colors['features'][bc.FEATURE_GAP])
            elif feature.type == bc.INSDC_FEATURE_ORIGIN_REPLICATION:
                gc_skew_track.xticks([(feature.location.start + feature.location.end)/2], outer=False, label_size=text_size/2, labels=['oriC'], label_orientation='vertical')  # oriC/V
            elif feature.type == bc.INSDC_FEATURE_ORIGIN_TRANSFER:
                gc_skew_track.xticks([(feature.location.start + feature.location.end)/2], outer=False, label_size=text_size/2, labels=['oriT'], label_orientation='vertical')  # oriT
            else:
                track.genomic_features([feature], fc=colors['features']['misc'])
    
        # plot GC content and GC skew
        seq = str(seqid2seq[sector.name])
        build_gc_content_skew(seq, colors, gc_content_track, gc_skew_track)

    fig = circos.plotfig(dpi=plot_dpi, figsize=(plot_size,plot_size))
    build_legend(circos, colors, plot_size)
    return fig


def build_features_type_cog(data, sequence_list, clipped, colors, plot_label, plot_size, plot_dpi):
    # Get contig genome seqid & size, features dict
    total_sequence_length = sum([len(seq['nt']) for seq in data['sequences']])
    seqid2seq = {rec.id:rec.seq for rec in sequence_list}
    seqid2size = {rec.id:len(rec.seq) for rec in sequence_list}
    seqid2features = {rec.id:rec.features for rec in sequence_list}

    if plot_size == 4:
        text_size = 6
    elif plot_size == 8:
        text_size = 12
    else:
        text_size = 30
    
    if(clipped):
        circos = Circos(seqid2size, space=2, start=10, end=350)
        circos.text('...', r=99.5, color=colors['backbone'], ha='center', va='center', size=6)
    else:
        circos = Circos(seqid2size, space=2)
    circos.text(plot_label, r=7, size=text_size, linespacing=1.5)
    for sector in circos.sectors:
        # build tracks
        outer_track = sector.add_track((99.5, 100))
        feature_forward_track = sector.add_track((87, 97), r_pad_ratio=0.1)
        feature_reverse_track = sector.add_track((77, 87), r_pad_ratio=0.1)
        non_cds_feature_track = sector.add_track((67, 77), r_pad_ratio=0.1)
        gc_content_track = sector.add_track((57, 67))
        gc_skew_track = sector.add_track((45, 55))

        # plot outer track
        build_sequence_backbone_track(sector, outer_track, total_sequence_length, colors, plot_size)

        # plot feature tracks
        for feature in seqid2features[sector.name]:
            if feature.type == bc.INSDC_FEATURE_CDS:
                color = colors['features'][bc.FEATURE_CDS]
                dbxrefs = feature.qualifiers.get('db_xref', None)
                if dbxrefs is not None:
                    for dbxref in dbxrefs:
                        if bc.DB_PREFIX_COG in dbxref:
                            cog = dbxref.split(':')[1]
                            if len(cog) != 7:  # skip COG clusters
                                if len(cog) != 1:  # select first of two COG categories
                                    cog = cog[:1]
                                color = colors['cog-classes'].get(cog.upper(), colors['cog-classes']['S'])
                track = feature_forward_track if feature.location.strand == 1 else feature_reverse_track
                track.genomic_features([feature], fc=color)
            elif feature.type == bc.INSDC_FEATURE_T_RNA:
                non_cds_feature_track.genomic_features([feature], fc=colors['features'][bc.FEATURE_T_RNA])
            elif feature.type == bc.INSDC_FEATURE_TM_RNA:
                non_cds_feature_track.genomic_features([feature], fc=colors['features'][bc.FEATURE_TM_RNA])
            elif feature.type == bc.INSDC_FEATURE_R_RNA:
                non_cds_feature_track.genomic_features([feature], fc=colors['features'][bc.FEATURE_R_RNA])
            elif feature.type == bc.INSDC_FEATURE_NC_RNA:
                non_cds_feature_track.genomic_features([feature], fc=colors['features'][bc.FEATURE_NC_RNA])
            elif feature.type == bc.INSDC_FEATURE_REGULATORY:
                non_cds_feature_track.genomic_features([feature], fc=colors['features'][bc.FEATURE_NC_RNA_REGION])
            elif feature.type == bc.INSDC_FEATURE_REPEAT_REGION:
                non_cds_feature_track.genomic_features([feature], fc=colors['features'][bc.FEATURE_CRISPR])
            elif feature.type == bc.INSDC_FEATURE_GAP:
                non_cds_feature_track.genomic_features([feature], fc=colors['features'][bc.FEATURE_GAP])
            elif feature.type == bc.INSDC_FEATURE_ORIGIN_REPLICATION:
                gc_skew_track.xticks([(feature.location.start + feature.location.end)/2], outer=False, label_size=text_size/2, labels=['oriC'], label_orientation='vertical')  # oriC/V
            elif feature.type == bc.INSDC_FEATURE_ORIGIN_TRANSFER:
                gc_skew_track.xticks([(feature.location.start + feature.location.end)/2], outer=False, label_size=text_size/2, labels=['oriT'], label_orientation='vertical')  # oriT
            else:
                non_cds_feature_track.genomic_features([feature], fc=colors['features']['misc'])
    
        # plot GC content and GC skew
        seq = str(seqid2seq[sector.name])
        build_gc_content_skew(seq, colors, gc_content_track, gc_skew_track)

    fig = circos.plotfig(dpi=plot_dpi, figsize=(plot_size,plot_size))
    build_legend(circos, colors, plot_size)
    return fig


def build_sequence_backbone_track(sector, outer_track, total_sequence_length, colors, plot_size):
    outer_track.axis(fc=colors['backbone'])
    label_formatter=lambda v: f'{v / 1000:.1f} kbp'
    if total_sequence_length >= 1_000_000:
        major_interval = 500_000
        minor_interval = int(major_interval / 5)
        label_formatter=lambda v: f'{v / 1_000_000:.1f} Mbp'
    elif total_sequence_length >= 100_000:
        major_interval = 100_000
        minor_interval = int(major_interval / 10)
    elif total_sequence_length >= 10_000:
        major_interval = 10_000
        minor_interval = int(major_interval / 10)
    elif total_sequence_length >= 1_000:
        major_interval = 1_000
        minor_interval = int(major_interval / 10)
    else:
        major_interval = 100
        minor_interval = int(major_interval / 10)
    
    if plot_size == 4:
        text_size = 4
    elif plot_size == 8:
        text_size = 8
    else:
        text_size = 16
    
    if sector.size > minor_interval:
        outer_track.xticks_by_interval(major_interval, label_formatter=label_formatter, label_size=text_size)
        outer_track.xticks_by_interval(minor_interval, tick_length=1, show_label=False, label_size=text_size)


def build_gc_content_skew(sequence, colors, gc_content_track, gc_skew_track):
    pos_list, gc_contents = calc_gc_content(sequence)
    gc_contents = gc_contents - gc_fraction(sequence) * 100
    positive_gc_contents = np.where(gc_contents > 0, gc_contents, 0)
    negative_gc_contents = np.where(gc_contents < 0, gc_contents, 0)
    abs_max_gc_content = np.max(np.abs(gc_contents))
    gc_content_track.fill_between(pos_list, positive_gc_contents, 0, vmin=-abs_max_gc_content, vmax=abs_max_gc_content, color=colors['gc-positive'])
    gc_content_track.fill_between(pos_list, negative_gc_contents, 0, vmin=-abs_max_gc_content, vmax=abs_max_gc_content, color=colors['gc-negative'])

    pos_list, gc_skews = calc_gc_skew(sequence)
    positive_gc_skews = np.where(gc_skews > 0, gc_skews, 0)
    negative_gc_skews = np.where(gc_skews < 0, gc_skews, 0)
    abs_max_gc_skew = np.max(np.abs(gc_skews))
    gc_skew_track.fill_between(pos_list, positive_gc_skews, 0, vmin=-abs_max_gc_skew, vmax=abs_max_gc_skew, color=colors['gc-skew-positive'])
    gc_skew_track.fill_between(pos_list, negative_gc_skews, 0, vmin=-abs_max_gc_skew, vmax=abs_max_gc_skew, color=colors['gc-skew-negative'])


def build_legend(circos, colors, plot_size):

    if plot_size == 4:
        text_size = 3
        marker_size = 1
    elif plot_size == 8:
        text_size = 6
        marker_size = 5
    else:
        text_size = 12
        marker_size = 9

    handles=[
        Patch(color=colors['features'][bc.FEATURE_CDS], label='CDS'),
        Patch(color=colors['features'][bc.FEATURE_T_RNA], label='tRNA'),
        Patch(color=colors['features'][bc.FEATURE_R_RNA], label='rRNA'),
        Patch(color=colors['features'][bc.FEATURE_NC_RNA], label='ncRNA'),
        Patch(color=colors['features'][bc.FEATURE_NC_RNA_REGION], label='ncRNA reg'),
        Patch(color=colors['features'][bc.FEATURE_CRISPR], label='CRISPR'),
        Line2D([], [], color=colors['gc-positive'], label="+ GC", marker="^", ms=marker_size, ls="None"),
        Line2D([], [], color=colors['gc-negative'], label="- GC", marker="v", ms=marker_size, ls="None"),
        Line2D([], [], color=colors['gc-skew-positive'], label="+ GC Skew", marker="^", ms=marker_size, ls="None"),
        Line2D([], [], color=colors['gc-skew-negative'], label="- GC Skew", marker="v", ms=marker_size, ls="None")
    ]
    _ = circos.ax.legend(
        handles=handles,
        bbox_to_anchor=(0.5, 0.4),
        loc='center',
        ncols=2,
        fontsize=text_size
    )


def calc_gc_content(seq: str):
    pos_list, gc_content_list = [], []
    window_size = int(len(seq) / 500)
    step_size = int(len(seq) / 1000)
    if window_size == 0 or step_size == 0:
        window_size, step_size = len(seq), int(len(seq) / 2)
    pos_list = list(range(0, len(seq), step_size)) + [len(seq)]
    for pos in pos_list:
        window_start_pos = pos - int(window_size / 2)
        window_end_pos = pos + int(window_size / 2)
        if window_start_pos < 0:
            window_start_pos = 0
        if window_end_pos > len(seq):
            window_end_pos = len(seq)
        subseq = seq[window_start_pos:window_end_pos]
        gc_content = gc_fraction(subseq) * 100
        gc_content_list.append(gc_content)
    return np.array(pos_list).astype(np.int64), np.array(gc_content_list).astype(np.float64)


def calc_gc_skew(seq: str):
    pos_list, gc_skew_list = [], []
    window_size = int(len(seq) / 500)
    step_size = int(len(seq) / 1000)
    if window_size == 0 or step_size == 0:
        window_size, step_size = len(seq), int(len(seq) / 2)
    pos_list = list(range(0, len(seq), step_size)) + [len(seq)]
    for pos in pos_list:
        window_start_pos = pos - int(window_size / 2)
        window_end_pos = pos + int(window_size / 2)
        window_start_pos = 0 if window_start_pos < 0 else window_start_pos
        window_end_pos = len(seq) if window_end_pos > len(seq) else window_end_pos
        subseq = seq[window_start_pos:window_end_pos]
        g, c = subseq.count('G'), subseq.count('C')
        gc_skew = (g - c) / float(g + c) if (g + c) > 0 else 0.0
        gc_skew_list.append(gc_skew)
    return np.array(pos_list).astype(np.int64), np.array(gc_skew_list).astype(np.float64)


def gc_fraction(seq: str):
    gc = seq.count('C') + seq.count('G') + seq.count('S')
    length = len(seq)
    return 0 if length == 0 else gc / length

if __name__ == '__main__':
    main()
