import os
import json
import yaml
import tempfile
import sys
import logging
import subprocess as sp

from pathlib import Path
from Bio import SeqUtils

import bakta.utils as bu
import bakta.config as cfg
import bakta.constants as bc



cog_colors = {
    'J':	'#ff0000',	#Translation, ribosomal structure and biogenesis
    'A':	'#c2af58',	#RNA processing and modification
    'K':	'#ff9900',	#Transcription
    'L':	'#ffff00',	#Replication, recombination and repair
    'B':	'#ffc600',	#Chromatin structure and dynamics
    'D':	'#99ff00',	#Cell cycle control, cell division, chromosome partitioning
    'Y':	'#493126',	#Nuclear structure
    'V':	'#ff008a',	#Defense mechanisms
    'T':	'#0000ff',	#Signal transduction mechanisms
    'M':	'#9ec928',	#Cell wall/membrane/envelope biogenesis
    'N':	'#006633',	#Cell motility
    'Z':	'#660099',	#Cytoskeleton
    'W':	'#336699',	#Extracellular structures
    'U':	'#33cc99',	#Intracellular trafficking, secretion, and vesicular transport
    'O':	'#00ffff',	#Posttranslational modification, protein turnover, chaperones
    'C':	'#9900ff',	#Energy production and conversion
    'G':	'#805642',	#Carbohydrate transport and metabolism
    'E':	'#ff00ff',	#Amino acid transport and metabolism
    'F':	'#99334d',	#Nucleotide transport and metabolism
    'H':	'#727dcc',	#Coenzyme transport and metabolism
    'I':	'#5c5a1b',	#Lipid transport and metabolism
    'P':	'#0099ff',	#Inorganic ion transport and metabolism
    'Q':	'#ffcc99',	#Secondary metabolites biosynthesis, transport and catabolism
    'R':	'#ff9999',	#General function prediction only
    'S':	'#d6aadf'	#Function unknown
    }


def main():
    # parse options and arguments
    parser = bu.init_parser(sub_command='_plot')
    parser.add_argument('input', metavar='<input>', help='Bakta-Annotation in JSON-Format')
    parser.add_argument('--conf', '-c', action='store', default=None, dest='plot',
                                      help='YAML configuration file for plot configuration')
    arg_group_io = parser.add_argument_group('Input / Output')
    arg_group_io.add_argument('--output', '-o', action='store', default=os.getcwd(),
                              help='Output directory (default = current working directory)')
    arg_group_io.add_argument('--tmp_dir', '-t', action='store', default=Path(tempfile.mkdtemp()).resolve(),
                              help='Config directory (default = temporary directory')
    arg_group_io.add_argument('--prefix', '-p', action='store', default='Plot', help='Prefix for output files')
    args = parser.parse_args()

    log = logging.getLogger('PLOT')

    #Check Input
    try:
        if (args.input == ''):
            raise ValueError('File path argument must be non-empty')
        annotation_path = Path(args.input).resolve()
        cfg.check_readability('annotation', annotation_path)
        cfg.check_content_size('annotation', annotation_path)
    except:
        log.error('provided input annotation file not valid! path=%s', args.input)
        sys.exit(f'ERROR: input annotation file ({args.input}) not valid!')
    log.info('input-path=%s', annotation_path)
    cfg.check_tmp_path(args)

    #Import Annotation Information
    with open(args.input, 'r') as input:
        annotation = json.load(input)
    features = annotation['features']
    contigs = annotation['sequences']

    #Import Configuration Information
    if args.plot != None:
        try:
            with open(args.plot) as conf:
                config = yaml.load(conf, Loader=yaml.FullLoader)
        except:
            log.error('wrong configuration file format!')
            sys.exit('ERROR: wrong configuration file format!')

    try:
        pgcc = config['pgcc']
    except:
        pgcc = '#CC6458'
    try:
        ngcc = config['ngcc']
    except:
        ngcc = '#43CC85'
    try:
        pgcs = config['pgcs']
    except:
        pgcs = '#CCBE6C'
    try:
        ngcs = config['ngcs']
    except:
        ngcs = '#5A4ECC'

    try:
        plots_string = config['plot']
        plots = plots_string.replace('-', '')
        plots = plots.replace("'", '')
        plots = plots.replace('"', '')
        plots = plots.split()
        plot_count = 1
        contig_selection = False
        for p in plots:
            p = p.split('/')
            plot_contig = []
            for contig in contigs:
                if 'all' in p:
                    plot_contig = [contig]
                    print(f'drawing plot {plot_count}')
                    write_plot(features,
                               plot_contig,
                               args.output,
                               args.tmp_dir,
                               plot_count,
                               args.prefix,
                               pgcc,
                               ngcc,
                               pgcs,
                               ngcs)
                    plot_count += 1
                    continue
                if contig['id'] in p:
                    plot_contig.append(contig)
                    contig_selection = True
            if contig_selection == True:
                print(f'drawing plot {plot_count}')
                write_plot(features,
                                  plot_contig,
                                  args.output,
                                  args.tmp_dir,
                                  plot_count,
                                  args.prefix,
                                  pgcc,
                                  ngcc,
                                  pgcs,
                                  ngcs)
                plot_count += 1
        print('all plots are finished!')

    except:
        plot_count = 0
        print('drawing plot')
        write_plot(features,
                   contigs,
                   args.output,
                   args.tmp_dir,
                   plot_count,
                   args.prefix,
                   pgcc,
                   ngcc,
                   pgcs,
                   ngcs)
        print('plot has been drawn!')

def write_plot(features,
               contigs,
               outdir,
               config_dir = None,
               plot_count = 0,
               prefix = None,
               p_gc_content_color = '#CC6458',
               n_gc_content_color = '#43CC85',
               p_gc_skew_color = '#CCBE6C',
               n_gc_skew_color = '#5A4ECC'):
    #######################
    # Define Variables
    #######################
    if prefix == None:
        prefix = cfg.prefix
    added_sequence_length = 0
    added_sequences = ""
    for c in contigs:
        added_sequence_length += int(c['length'])
        added_sequences += c['sequence']
    window_size = added_sequence_length/100 if added_sequence_length < 2000 else added_sequence_length/1000
    window_size = 3 if window_size < 3 else window_size
    step_size = window_size* 0.2 if window_size >= 5 else 1
    track_radius = 1.0
    gc_radius = 0.2

    label_prefix = 'b'
    multiplier = 1
    if added_sequence_length > 1000:
        label_prefix = 'kb'
        multiplier = 0.001
    if added_sequence_length > 1000000:
        label_prefix = 'mb'
        multiplier = 0.000001


    ###########################
    # Config Paths
    ###########################
    if config_dir == None:
        config_dir = cfg.tmp_path.joinpath(f'circos_config_files')
    Path(config_dir).mkdir(parents=True, exist_ok=True)
    main_conf = f'{config_dir}/main.conf'
    karyotype_txt = f'{config_dir}/karyotype.txt'
    ideogram_conf = f'{config_dir}/ideogram.conf'
    ticks_conf = f'{config_dir}/ticks.conf'
    tracks_conf = f'{config_dir}/tracks.conf'

    f_cds_file = f'{config_dir}/f_cds.txt'
    r_cds_file = f'{config_dir}/r_cds.txt'
    nc_file = f'{config_dir}/nc.txt'
    gc_content_file = f'{config_dir}/gc_content.txt'
    gc_skew_file = f'{config_dir}/gc_skew.txt'

    #####################
    # write feature files
    #####################
    f_cds_contents = ""
    r_cds_contents = ""
    nc_contents = ""
    for f in features:
        contig, start, stop, color = f['contig'], f['start'], f['stop'], '#888888'
        if f['type'] == 'cds':
            try:
                psc = f['psc']
                cog = psc['cog_category']
                if len(cog) != 1:
                    cog = cog[:1]
                color = cog_colors[cog]
            except:
                pass
            if f['strand'] == bc.STRAND_FORWARD:
                f_cds_contents += f"{contig} {start} {stop} {f['strand']} color={hex_to_rgb(color)}\n"
            if f['strand'] == bc.STRAND_REVERSE:
                r_cds_contents += f"{contig} {start} {stop} {f['strand']} color={hex_to_rgb(color)}\n"
            else:
                continue
        else:
            nc_contents += f"{contig} {start} {stop} {f['strand']} color={hex_to_rgb(color)}\n"

    with open(f_cds_file, 'w') as f:
        f.write(f_cds_contents)
    with open(r_cds_file, 'w') as f:
        f.write(r_cds_contents)
    with open(nc_file, 'w') as f:
        f.write(nc_contents)

    #####################################
    # write gc content unf gc skew files
    #####################################
    gc_content_text = ""
    gc_skew_text = ""
    max_gc_content = 0
    max_gc_skew = 0
    gc_mean = SeqUtils.GC(added_sequences)
    for contig in contigs:
        seq = contig['sequence']
        step_size = int(step_size)
        window_size = int(window_size)
        for w in range(0, len(seq), step_size):

            start = w - (window_size / 2)
            start = start + len(seq) if start < 0 else start
            stop = w + (window_size / 2)
            stop = stop - len(seq) if stop > len(seq) else stop
            start = int(start)
            stop = int(stop)

            if start < stop:
                subseq = seq[start:stop]
            else:
                subseq = "".join([seq[start:], seq[:stop]])
            gc_value = gc_mean - SeqUtils.GC(subseq)
            max_gc_content = abs(gc_value) if max_gc_content < abs(gc_value) else max_gc_content
            content_color = p_gc_content_color if gc_value > 0 else n_gc_content_color
            gc_content_text += f"{contig['id']} {w} {w} {gc_value} fill_color={hex_to_rgb(content_color)}\n"
            G = float(subseq.count('G'))
            C = float(subseq.count('C'))
            try:
                 gc_skew = (G - C) / (G + C)
            except ZeroDivisionError:
                gc_skew = 0.0
            max_gc_skew = abs(gc_skew) if max_gc_skew < abs(gc_skew) else max_gc_skew
            skew_color = p_gc_skew_color if gc_skew > 0 else n_gc_skew_color
            gc_skew_text += f"{contig['id']} {w} {w} {gc_skew} fill_color={hex_to_rgb(skew_color)}\n"

    with open(gc_content_file, 'w') as f:
        f.write(gc_content_text)
    with open(gc_skew_file, 'w') as f:
        f.write(gc_skew_text)


    ##############################
    # write configurationfiles
    ##############################
    #write main config
    if plot_count == 0:
        plot_count = ""
    main_config_text = f'''
    karyotype                   = {karyotype_txt}
    chromosomes_units           = {added_sequence_length}
    chromosomes_display_default = yes
    <plots>
    <<include {tracks_conf}>>
    </plots>
    <image>
    <<include image.conf>>
    file*                       = {prefix}{plot_count}
    dir*                        = {outdir}
    </image> 
    <<include {ideogram_conf}>>
    <<include {ticks_conf}>>
    <<include etc/colors_fonts_patterns.conf>>
    <<include etc/housekeeping.conf>>
    '''
    with open(main_conf, 'w') as f:
        f.write(main_config_text)

    #write_karyotype_file
    karyotype_text = ""
    for i, c in enumerate(contigs):
        color = f'#{khp(str(hex(0x69 + (i%6)*0x28)))+khp(str(hex(0x69 + (i%6)*0x28)))+khp(str(hex(0x69 + (i%6)*0x28)))}'
        karyotype_text += f"chr - {c['id']} {i + 1} 0 {c['length']} {hex_to_rgb(color)}\n"
    with open(karyotype_txt, 'w') as f:
        f.write(karyotype_text)

    #write_ideogram_config():
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
    with open(ideogram_conf, 'w') as f:
        f.write(ideogram_text)

    #write_ticks_config:
    ticks_text = f'''
    show_ticks       = yes
    show_tick_labels = yes
    <ticks>
    radius      = 1r
    color       = black
    thickness   = 2p
    multiplier  = {multiplier}
    orientation = out
    format      = %.01f{label_prefix}
    <tick>
    spacing      = 0.15u
    show_label   = yes
    label_size   = 35
    label_offset = 10p
    size         = 25p
    thickness    = 4p
    </tick>
    <tick>
    spacing      = 0.015u
    show_label   = yes
    label_size   = 15
    size         = 15p
    thickness    = 3p
    </tick>
    <tick>
    spacing      = 0.0015u
    size         = 5p
    </tick>
    </ticks>
    '''
    with open(ticks_conf, 'w') as f:
            f.write(ticks_text)


    # write trackconfiguration
    tracks_text = ""
    feature_files = [f_cds_file, r_cds_file, nc_file]
    for feature in feature_files:
        tracks_text += f'''
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
        track_radius -= 0.1

    gc_files = [gc_content_file,gc_skew_file]
    for gc_file in gc_files:
        maximum = max_gc_content if gc_file == gc_content_file else max_gc_skew
        minimum = -maximum
        tracks_text += f'''
        <plot>
        type = histogram
        file = {gc_file}
        r1 = {track_radius}r
        r0 = {track_radius-gc_radius}r
        min = {minimum}
        max = {maximum}
        thickness = 0
        orientation = out
        </plot>    
        '''

        track_radius -= gc_radius
    with open(tracks_conf, 'w') as f:
        f.write(tracks_text)
    #############
    # run Circos
    #############
    sp.run(['circos', '-conf', main_conf],stdout=sp.DEVNULL, stderr=sp.STDOUT)
    return

def hex_to_rgb(hex):
    hex = hex.replace('#','')
    rgb = ''
    for i in (0,2,4):
        rgb += f',{str(int(hex[i:i+2],16))}'
    rgb = rgb.replace(',','',1)
    return rgb

def khp(hex_with_prefix):
    hex_without_prefix = hex_with_prefix[2:]
    return hex_without_prefix


if __name__ == '__main__':
    main()
