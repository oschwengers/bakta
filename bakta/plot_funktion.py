from pathlib import Path
from subprocess import run
from Bio import SeqUtils

import bakta.config as cfg

genome_lengths = []
step_sizes = []


track_radius = 1.0


p_gc_content_color = 'red'
n_gc_content_color = 'yellow'
p_gc_skew_color = 'blue'
n_gc_skew_color = 'green'

def write_plot(features, contigs, plot_path):
    plot_data = {}

    plot_data['contigs'] = contigs
    plot_data['features'] = features
    plot_data['outdir'] = plot_path

    get_genome_lengths(contigs)
    # plot_data['genome_lengths'] = get_genome_lengths(contigs)
    plot_data['track_radius'] = 1.0

    plot_data['config_dir'] = cfg.tmp_path.joinpath(f'circos_config_files')
    Path(plot_data['config_dir']).mkdir(parents=True, exist_ok=True)
    config_files = {}
    config_files['main.conf'] = f'{plot_data["config_dir"]}/main.conf'
    config_files['karyotype.txt'] = f'{plot_data["config_dir"]}/karyotype.txt'
    config_files['ideogram.conf'] = f'{plot_data["config_dir"]}/ideogram.conf'
    config_files['ticks.conf'] = f'{plot_data["config_dir"]}/ticks.conf'
    config_files['tracks.conf'] = f'{plot_data["config_dir"]}/tracks.conf'
    plot_data['config_files'] = config_files
    track_files = {}
    track_files['f_cds_file'] = f'{plot_data["config_dir"]}/f_cds.txt'
    track_files['r_cds_file'] = f'{plot_data["config_dir"]}/r_cds.txt'
    track_files['nc_file'] = f'{plot_data["config_dir"]}/nc.txt'
    track_files['gc_content_file'] = f'{plot_data["config_dir"]}/gc_content.txt'
    track_files['gc_skew_file'] = f'{plot_data["config_dir"]}/gc_skew.txt'
    plot_data['track_files'] = track_files

    plot_data['window_sizes'] = get_window_sizes(contigs)
    plot_data['step_sizes'] = get_step_sizes(plot_data['window_sizes'])

    # write configurationfiles
    write_circos_config(plot_data)

    # write trackfiles
    write_feature_track(plot_data['config_files'], track_files['f_cds_file'])
    write_feature_track(plot_data['config_files'], track_files['r_cds_file'])
    write_feature_track(plot_data['config_files'], track_files['nc_file'])
    write_gc_track(plot_data['config_files'], track_files['gc_content_file'])
    write_gc_track(plot_data['config_files'], track_files['gc_skew_file'])

    # write featurefiles
    write_feature_files(plot_data)
    write_gc_files(plot_data)

    # run Circos
    run(['circos', '-conf', config_files['main.conf']])
    return




def write_circos_config(pd):
    config_files = pd['config_files']
    #write main config
    main_config_text = f'''
    karyotype                   = {config_files["karyotype.txt"]}
    chromosomes_units           = 1000000
    chromosomes_display_default = yes
    <image>
    <<include image.conf>>
    file*                       = {cfg.prefix}
    dir*                        = {pd["outdir"]}
    </image> 
    <plots>
    <<include {config_files["tracks.conf"]}>>
    </plots>
    <<include {config_files["ideogram.conf"]}>>
    <<include {config_files["ticks.conf"]}>>
    <<include etc/colors_fonts_patterns.conf>>
    <<include etc/housekeeping.conf>>
    '''
    with open(config_files['main.conf'], 'w') as f:
        f.write(main_config_text)

    #write_karyotype_file
    karyotype_text = ""
    count = 0
    colors = ['black', 'dgrey', 'grey', 'red']
    for i, c in enumerate(pd['contigs']):
        karyotype_text += f"chr - {c['id']} {i + 1} {count} {count + len(c['sequence'])} {colors[i]}\n"
        count += len('sequence')
    with open(config_files['karyotype.txt'], 'w') as f:
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
    with open(config_files['ideogram.conf'], 'w') as f:
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
    with open(config_files['ticks.conf'], 'w') as f:
            f.write(ticks_text)


def write_feature_track(config_files, track_file, track_radius):
    feature_track = f'''
    <plot>
    type             = tile
    file             = {track_file}
    r1               = {track_radius}r
    r0               = {track_radius-0.1}r
    orientation      = out
    layers           = 1
    margin           = 0.01u
    thickness        = 5
    padding          = 1
    stroke_color     = black
    stroke_thickness = 0
    layers_overflow  = collapse
    </plot>
    '''
    track_radius -= 0.1
    with open(config_files['tracks.conf'], 'a') as f:
        f.write(feature_track)
    return track_radius

def write_feature_files(pd):
    track_files = pd['track_files']
    f_cds_contents = ""
    r_cds_contents = ""
    nc_contents = ""
    for f in pd['features']:
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

    with open(track_files['f_cds_file'], 'w') as f:
        f.write(f_cds_contents)
    with open(track_files['r_cds_file'], 'w') as f:
        f.write(r_cds_contents)
    with open(track_files['nc_file'], 'w') as f:
        f.write(nc_contents)


def write_gc_track(config_files, track_file):
    track_radius = track_radius
    gc_track = f'''
    <plot>
    type = histogram
    file = {track_file}
    r1 = {track_radius}r
    r0 = {track_radius-0.1}r
    thickness = 0
    orientation = out
    </plot>    
    '''

    track_radius -= 0.1
    with open(config_files['tracks.conf'], 'a') as f:
        f.write(gc_track)
    return track_radius

def write_gc_files(pd):
    track_files = pd['track_files']
    gc_content_text = ""
    gc_skew_text = ""

    x = zip(pd['contigs'],pd['step_sizes'], pd['window_sizes'])
    for contig, step_size, window_size in x:
        seq = contig['sequence']
        gc_mean = SeqUtils.GC(seq)
        step_size = int(step_size)
        window_size = int(window_size)
        for w in range(0, len(seq), step_size):
            start = w - (window_size / 2)
            start = start + len(seq) if start < 0 else start
            end = w + (window_size / 2)
            end = end - len(seq) if end > len(seq) else end
            start = int(start)
            end = int(end)
            if start < end:
                subseq = seq[start:end]
            else:
                subseq = "".join([seq[start:], seq[:end]])
            gc_value = gc_mean - SeqUtils.GC(subseq)
            content_color = p_gc_content_color if gc_value > 0 else n_gc_content_color
            gc_content_text += f"{contig['id']} {w} {w + step_size} {gc_value} fill_color={content_color}\n"

            G = float(subseq.count('G'))
            C = float(subseq.count('C'))
            try:
                gc_skew = (G - C)/(G + C)
            except ZeroDivisionError:
                gc_skew = 0.0
            skew_color = p_gc_skew_color if gc_skew > 0 else n_gc_skew_color
            gc_skew_text += f"{contig['id']} {w} {w + step_size} {gc_skew} fill_color={skew_color}\n"
    with open(track_files['gc_content_file'], 'w') as f:
        f.write(gc_content_text)
    with open(track_files['gc_skew_file'], 'w') as f:
        f.write(gc_skew_text)


def get_window_sizes(contigs):
    window_sizes = []
    for i in contigs:
        window_sizes.append(2 if i['length'] < 2000 else i['length']/1000)
    return window_sizes

def get_step_sizes(window_sizes):
    for i in window_sizes:
        step_size = i*0.2 if i >= 5 else 1
        step_sizes.append(step_size)
        print(step_size)
    return step_sizes

def get_genome_lengths(contigs):
    for i in contigs:
        genome_lengths.append(i['length'])
    return genome_lengths





