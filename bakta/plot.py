from pathlib import Path
from subprocess import run
from Bio import SeqUtils

import bakta.config as cfg


def write_plot(features,
               contigs,
               outdir,
               all_plots = False,
               plot_count = 0,
               p_gc_content_color = 'red',
               n_gc_content_color = 'yellow',
               p_gc_skew_color = 'blue',
               n_gc_skew_color = 'green'
                ):
    ##########################
    # 1Plot/sequence on/off
    ##########################
    if all_plots == True:
        outdir = outdir.joinpath('Plots')
        Path(outdir).mkdir(parents=True, exist_ok=True)
        for c in contigs:
            if c['complete'] == True:
                contig_features = []
                contig_contigs = [c]
                for f in features:
                    if f['contig'] == c['id']:
                        contig_features.append(f)
                plot_count += 1

                write_plot(contig_features,
                            contig_contigs,
                            outdir,
                            False,
                            plot_count,
                            p_gc_content_color,
                            p_gc_skew_color,
                            n_gc_content_color,
                            n_gc_skew_color)

    #######################
    # Define Variables
    #######################
    genome_lengths = 0
    window_sizes = []
    step_sizes = []
    for c in contigs:
        genome_lengths += int(c['length'])
        window_size = c['length']/100 if c['length'] < 2000 else c['length'] / 1000
        window_size = 10 if window_size < 10 else window_size
        window_sizes.append(window_size)
        step_size = window_size* 0.2 if window_size >= 5 else 1
        step_sizes.append(step_size)

    track_radius = 1.0

    ###########################
    # Config Paths
    ###########################
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
    # write featurefiles
    #####################
    f_cds_contents = ""
    r_cds_contents = ""
    nc_contents = ""
    for f in features:
        color = 'black'
        contig, start, stop = f['contig'], f['start'], f['stop']

        if f['type'] == 'cds':
            if f['strand'] == '+':
                f_cds_contents += f"{contig} {start} {stop} {f['strand']} color={color}\n"
            if f['strand'] == '-':
                r_cds_contents += f"{contig} {start} {stop} {f['strand']} color={color}\n"
            else:
                continue
        else:
            nc_contents += f"{contig} {start} {stop} {f['strand']} color={color}\n"

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
    x = zip(contigs, step_sizes, window_sizes)
    for contig, step_size, window_size in x:
        seq = contig['sequence']
        gc_mean = SeqUtils.GC(seq)
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
            gc_content_text += f"{contig['id']} {w} {w + step_size} {gc_value} fill_color={content_color}\n"
            G = float(subseq.count('G'))
            C = float(subseq.count('C'))
            try:
                 gc_skew = (G - C) / (G + C)
            except ZeroDivisionError:
                gc_skew = 0.0
            max_gc_skew = abs(gc_skew) if max_gc_skew < abs(gc_skew) else max_gc_skew
            skew_color = p_gc_skew_color if gc_skew > 0 else n_gc_skew_color
            gc_skew_text += f"{contig['id']} {w} {w + step_size} {gc_skew} fill_color={skew_color}\n"

    with open(gc_content_file, 'w') as f:
        f.write(gc_content_text)
    with open(gc_skew_file, 'w') as f:
        f.write(gc_skew_text)


    ##############################
    # write configurationfiles
    ##############################
    #write main config
    main_config_text = f'''
    karyotype                   = {karyotype_txt}
    chromosomes_units           = {10 ** (len(str(genome_lengths)) - 1)}
    chromosomes_display_default = yes
    <plots>
    <<include {tracks_conf}>>
    </plots>
    <image>
    <<include image.conf>>
    file*                       = {cfg.prefix}{plot_count}
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
    count = 0
    colors = [ 'dgrey', 'grey', 'red', 'black']
    for i, c in enumerate(contigs):
        karyotype_text += f"chr - {c['id']} {i + 1} {count} {len(c['sequence'])} {colors[i]}\n"

    with open(karyotype_txt, 'w') as f:
        f.write(karyotype_text)

    #write_ideogram_config():
    ideogram_text = '''
    <ideogram>
    <spacing>
    default = 0.005r
    </spacing>
    radius           = 0.80r
    thickness        = 15p
    fill             = yes
    stroke_color     = dgrey
    show_bands       = yes
    fill_bands       = yes
    stroke_thickness = 2p
    show_label       = no
    label_font       = default
    label_radius     = 1r + 75p
    label_size       = 30
    label_parallel   = yes
    </ideogram>
    '''
    with open(ideogram_conf, 'w') as f:
        f.write(ideogram_text)

    #write_ticks_config:
    ticks_text = '''
    show_ticks       = yes
    show_tick_labels = yes
    <ticks>
    radius      = 1r
    color       = black
    thickness   = 2p
    multiplier  = 1e-6
    orientation = out
    format      = %.1f Mb
    <tick>
    spacing      = 0.5u
    show_label   = yes
    label_size   = 35
    label_offset = 10p
    size         = 25p
    </tick>
    <tick>
    spacing      = 0.5u
    show_label   = no
    size         = 15p
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
        thickness        = 70
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
        r0 = {track_radius-0.15}r
        min = {minimum}
        max = {maximum}
        thickness = 0
        orientation = out
        </plot>    
        '''

        track_radius -= 0.15
    with open(tracks_conf, 'w') as f:
        f.write(tracks_text)

    #############
    # run Circos
    #############
    run(['circos', '-conf', main_conf])
    return
