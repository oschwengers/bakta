from subprocess import run
from pathlib import Path
from typing import List, Optional
from Bio import SeqUtils
import bakta.config as cfg


class Plot_Data:
    def __init__(self, features, contigs, plot_path):
        self.features = features
        self.contigs = contigs
        self.plot_path = plot_path
        self.outdir = plot_path  # cfg.tmp_path

        #Genom-Daten
        self.genom_length = 0
        self.get_genom_length()
        self.genom_seq = ""
        self.get_genom_seq()
        self.contig_seqs = []
        self.get_contig_seqs()
        self.average_gc = 0
        self.get_average_gc()

        #radius
        self._r = 1.0
        self.gc_skew_r = 0.15
        self.gc_content_r = 0.15
        #colors
        self.f_cds_color = 'black'
        self.r_cds_color = 'grey'
        self.nc_color = 'dgreen'
        self.gc_content_p_color = 'green'
        self.gc_content_n_color = 'blue'
        self.gc_skew_p_color = 'yellow'
        self.gc_skew_n_color = 'red'
        #trackconfig
        self._track_config = ""


        #TEMPFILES
        # Circos config directory
        self.outdir.mkdir(exist_ok=True)
        self.config_dir = self.outdir / "circos_config"
        self.config_dir.mkdir(exist_ok=True)
        self.ref_features_dir = self.config_dir / "reference_features"
        self.ref_features_dir.mkdir(exist_ok=True)
        # Circos config files
        self.config_file = self.config_dir / "circos.conf"
        self.ideogram_file = self.config_dir / "ideogram.conf"
        self.ticks_file = self.config_dir / "ticks.conf"
        self.karyotype_file = self.config_dir / "karyotype.txt"
        # Features config files
        self.f_cds_file = self.ref_features_dir / "forward_cds.txt"
        self.r_cds_file = self.ref_features_dir / "reverse_cds.txt"
        self.nc_file = self.ref_features_dir / "nc_file.txt"
        #GC content & GC skew config file
        self.gc_content_file = self.config_dir / "gc_content.txt"
        self.gc_skew_file = self.config_dir / "gc_skew.txt"





    ####################################################################
    ###Create Config-Files
    ####################################################################
    def write_config_file(self) -> None:
        #write config files
        self.write_ideogram_conf()
        self.write_ticks_conf()
        self.write_karyotype_file()

        #write and add feature files
        self.add_feature_track( self.f_cds_file, "cds", '+', self.f_cds_color)
        self.add_feature_track(self.r_cds_file, "cds", '-', self.r_cds_color)
        self.add_feature_track(self.nc_file, "else", '*', self.nc_color)
        self.add_gc_content_track()
        self.add_gc_skew_track()



        # write main config file
        config_contents = self._concat_lines(
            [
                "karyotype         = {0}".format(self.karyotype_file),
                "chromosomes_units = {0}".format(self._chromosome_units),

                "<<include {0}>>".format(self.ideogram_file),
                "<<include {0}>>".format(self.ticks_file),

                "<plots>",
                self._track_config.rstrip("\n"),
                "</plots>",

                "<image>",
                "<<include image.conf>>",
                "file* = {0}.png".format(cfg.prefix),
                "dir*  = {0}".format(self.plot_path),
                "</image>",

                "<<include colors_fonts_patterns.conf>> ",
                "<<include housekeeping.conf>> ",
            ]
        )
        with open(self.config_file, "w") as f:
            f.write(config_contents)

    #write karyotype file
    def write_karyotype_file(self) -> None:
        """Write karyotype txt"""
        contents = "" #f"chr - main 1 0 {self.genom_length} grey\n"
        colors = ["lgrey", "dgrey"]
        base_len = 0
        for idx, contig_seq in enumerate(self.contig_seqs):
            start, end, color = base_len, base_len + len(contig_seq), colors[idx % 2]
            contents += f"chr - contig_{idx + 1} {idx + 1} {start} {end} {color}\n"
            base_len += len(contig_seq)
        with open(self.karyotype_file, "w") as f:
            f.write(contents)

    #write ideogram file
    def write_ideogram_conf(self) -> None:
        """Write Circos Ideogram config"""
        contents = self._concat_lines(
            [
                "<ideogram>",
                "<spacing>",
                "default = 0.005r",
                "</spacing>",
                "radius           = 0.80r",
                "thickness        = 15p",
                "fill             = yes",
                "stroke_color     = dgrey",
                "show_bands       = yes",
                "fill_bands       = yes",
                "stroke_thickness = 2p",
                "show_label       = no",
                "label_font       = default",
                "label_radius     = 1r + 75p",
                "label_size       = 30",
                "label_parallel   = yes",
                "</ideogram>",
            ]
        )
        with open(self.ideogram_file, "w") as f:
            f.write(contents)

    #write ticks file
    def write_ticks_conf(self) -> None:
        """Write Circos ticks config"""
        contents = self._concat_lines(
            [
                "show_ticks       = yes",
                "show_tick_labels = yes",
                "<ticks>",
                "radius       = 1r",
                "color        = black",
                "thickness    = 2p",
                "multiplier   = {0}".format(self._ticks_multiplier),
                "orientation  = out",
                "format       = {0} {1}".format(self._ticks_format, self._ticks_unit),
                "<tick>",
                "spacing      = {0:.2f}u".format(self._largeticks_spacing),
                "show_label   = yes",
                "label_size   = 35p",
                "label_offset = 10p",
                "size         = 25p",
                "</tick>",
                "<tick>",
                "spacing      = {0:.2f}u".format(self._smallticks_spacing),
                "show_label   = no",
                "size         = 15p",
                "</tick>",
                "</ticks>",
            ]
        )
        with open(self.ticks_file, "w") as f:
            f.write(contents)



    ####################################################################
    ###Create Feature-Files
    ####################################################################
    #add feature track to config
    def add_feature_track(
        self,
        feature_file: Path,
        feature_type: str,
        target_strand: Optional[str] = None,
        color: str = "grey",
        feature_r: float = 0.07,
    ) -> None:
        """Add Feature track

        Args:
            feature_file (Path): Feature file to write
            feature_type (str): Feature type (e.g. 'CDS', 'rRNA', 'tRNA')
            target_strand (Optional[int]): Strand ('1', '-1', 'None')
            color (str): Feature color to be drawn
            feature_r (float): Feature radius size
        """
        self._write_feature_file(feature_file, feature_type, target_strand, color)
        self._track_config += self._concat_lines(
            [
                f"##### {feature_type} Feature Track #####",
                "<plot>",
                "type             = tile",
                "file             = {0}".format(feature_file),
                "r1               = {0:.3f}r".format(self._r),
                "r0               = {0:.3f}r".format(self._r - feature_r),
                "orientation      = out",
                "layers           = 1",
                "margin           = 0.01u",
                "thickness        = {0}".format(feature_r * 1000),
                "padding          = 1",
                "stroke_color     = black",
                "stroke_thickness = 0",
                "layers_overflow  = collapse",
                "</plot>",
            ]
        )
        self._r -= feature_r

    #write feature file
    def _write_feature_file(self,
        feature_file: Path,
        feature_type: str,
        target_strand: Optional[int],
        color: str
    ) -> None:

        contents = ""
        if feature_type != 'else':
            for f in self.features:
                if f['type'] != feature_type:
                    continue
                if f['strand'] != target_strand:
                    continue
                contig, start, end, strand = f['contig'], f['start'], f['stop'], f['strand']
                contents += f"{contig} {start} {end} {strand} color={color}\n"
        for f in self.features:
            if f['type'] == 'cds':
                continue
            contig, start, end, strand = f['contig'], f['start'], f['stop'], f['strand']
            contents += f"{contig} {start} {end} {strand} color={color}\n"
        with open(feature_file, "w") as f:
            f.write(contents)

    #add gc content track to config
    def add_gc_content_track(self) -> None:
        """Add GC Content track"""
        self._r = self._boundary if self._r > self._boundary else self._r
        abs_max_value = self._write_gc_content_file()
        self._track_config += self._concat_lines(
            [
                "##### GC Content Track #####",
                "<plot>",
                "type        = histogram",
                "file        = {0}".format(self.gc_content_file),
                "r1          = {0:.3f}r".format(self._r),
                "r0          = {0:.3f}r".format(self._r - self.gc_content_r),
                "min         = {0:.3f}".format(-abs_max_value),
                "max         = {0:.3f}".format(abs_max_value),
                "thickness   = 0",
                "orientation = out",
                "</plot>",
            ]
        )
        self._r -= self.gc_content_r

    #write gc content file
    def _write_gc_content_file(self) -> float:
        """Write GC Content file"""
        contents = ""
        max_gc_content = 0
        for c in self.contigs:
            gc_content_values = self.gc_content(c['sequence'], self._window_size, self._step_size)
            gc_content_values = [v - self.average_gc for v in gc_content_values]

            for i, value in enumerate(gc_content_values):
                start = i * self._step_size
                end = start + self._step_size
                end = self.genom_length if end > self.genom_length else end
                color = self.gc_content_p_color if value > 0 else self.gc_content_n_color
                contig = c['id']
                contents += f"{contig} {start} {end} {value} fill_color={color}\n"
                if abs(value) > max_gc_content:
                    max_gc_content = value
        with open(self.gc_content_file, "w") as f:
            f.write(contents)
        return max_gc_content

    #add gc skew track
    def add_gc_skew_track(self) -> None:
        """Add GC Skew track"""
        self._r = self._boundary if self._r > self._boundary else self._r
        abs_max_value = self._write_gc_skew_file()
        self._track_config += self._concat_lines(
            [
                "##### GC Skew Track #####",
                "<plot>",
                "type        = histogram",
                "file        = {0}".format(self.gc_skew_file),
                "r1          = {0:.3f}r".format(self._r),
                "r0          = {0:.3f}r".format(self._r - self.gc_skew_r),
                "min         = {0:.3f}".format(-abs_max_value),
                "max         = {0:.3f}".format(abs_max_value),
                "thickness   = 0",
                "orientation = out",
                "</plot>",
            ]
        )
        self._r -= self.gc_skew_r

    #write gc skew file
    def _write_gc_skew_file(self) -> float:
        """Write GC Skew file"""
        contents = ""
        max_gc_skew = 0
        for c in self.contigs:
            gc_skew_values = self.gc_skew(c['sequence'], self._window_size, self._step_size)
            for i, value in enumerate(gc_skew_values):
                start = i * self._step_size
                end = start + self._step_size
                end = self.genom_length if end > self.genom_length else end
                color = self.gc_skew_p_color if value > 0 else self.gc_skew_n_color
                contig = c['id']
                contents += f"{contig} {start} {end} {value} fill_color={color}\n"
                if abs(value) > max_gc_skew:
                    max_gc_skew = value
        with open(self.gc_skew_file, "w") as f:
            f.write(contents)
        return max_gc_skew

    ####################################################################
    ###Calculate Genom features
    ####################################################################
    def gc_content(self, contig, window_size: int = 5000, step_size: int = 2000) -> List[float]:
        """Calculate GC content in sliding window
        Args:
            window_size (int, optional): Window size
            step_size (int, optional): Step size
        Returns:
            List[float]: GC content values in sliding window"""
        gc_content_values = []
        seq = contig
        for i in range(0, len(seq), step_size):
            start_pos = i - int(window_size / 2)
            start_pos = 0 if start_pos < 0 else start_pos
            end_pos = i + int(window_size / 2)
            end_pos = len(seq) if end_pos > len(seq) else end_pos

            if start_pos > end_pos:
                subseq = ''.join([seq[start_pos:], seq[:end_pos]])
            else:
                subseq = seq[start_pos:end_pos]

            gc_content_values.append(SeqUtils.GC(subseq))

        return gc_content_values


    def gc_skew(self, contig, window_size: int = 5000, step_size: int = 2000) -> List[float]:
        """Calculate GC skew in sliding window
        Args:
            window_size (int, optional): Window size
            step_size (int, optional): Step size
        Returns:
            List[float]: GC skew values in sliding window"""
        gc_skew_values = []
        seq = contig
        for i in range(0, len(seq), step_size):
            start_pos = i - int(window_size / 2)
            start_pos = 0 if start_pos < 0 else start_pos
            end_pos = i + int(window_size / 2)
            end_pos = len(seq) if end_pos > len(seq) else end_pos

            subseq = seq[start_pos:end_pos]
            g = subseq.count("G") + subseq.count("g")
            c = subseq.count("C") + subseq.count("c")
            try:
                skew = (g - c) / float(g + c)
            except ZeroDivisionError:
                skew = 0.0
            gc_skew_values.append(skew)
        return gc_skew_values

    def get_genom_length(self):
        for i in self.contigs:
            self.genom_length += i['length']

    def get_genom_seq(self):
        for i in self.contigs:
            self.genom_seq += i['sequence']

    def get_contig_seqs(self):
        for i in self.contigs:
            self.contig_seqs.append(i['sequence'])

    def get_average_gc(self):
        self.average_gc = SeqUtils.GC(self.genom_seq)

    ####################################################################
    ###Config Essentials
    ####################################################################
    @property
    def _window_size(self) -> int:
        """Window size for GC content & GC skew calculation"""
        wisz = int(self.genom_length / 1000)
        if self.genom_length < 1000:
            wisz = self.genom_length / 100
        return wisz

    @property
    def _step_size(self) -> int:
        """Step size for GC content & GC skew calculation"""
        stsz = int(self._window_size * 0.4)
        if stsz == 0:
            stsz = 1
        return stsz

    @property
    def _chromosome_units(self) -> int:
        """Chromosome units"""
        return 10 ** (len(str(self.genom_length)) - 1)

    @property
    def _ticks_format(self) -> str:
        """Ticks format"""
        if self._chromosome_units >= 10**6:
            return "%.1f"
        else:
            return "%d"

    @property
    def _ticks_multiplier(self) -> float:
        """Ticks multiplier"""
        if self._chromosome_units >= 10**6:
            return 1e-6
        else:
            return 1e-3

    @property
    def _ticks_unit(self) -> str:
        """Ticks unit"""
        if self._chromosome_units >= 10**6:
            return "Mb"
        else:
            return "Kb"

    @property
    def _largeticks_spacing(self) -> float:
        """Largeticks spacing"""
        if self.genom_length / self._chromosome_units >= 2:
            return 0.5
        else:
            return 0.1

    @property
    def _smallticks_spacing(self) -> float:
        """Smallticks spacing"""
        if self.genom_length / self._chromosome_units >= 2:
            return 0.1
        else:
            return 0.01

    @property
    def _boundary(self) -> float:
        """Boundary radius for GC content & GC skew track"""
        return 0.7


    def _concat_lines(self, lines: List[str]) -> str:
        """Concatenate lines Args: lines (List[str]): Target lines
        Returns: str: Concatenated lines string"""
        return "\n".join(lines) + "\n"

    def _cog_colors(self):
        cog_colors = {
        'J':    '#ff0000',	#Translation, ribosomal structure and biogenesis
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
        'Q':    '#ffcc99',	#Secondary metabolites biosynthesis, transport and catabolism
        'R':    '#ff9999',	#General function prediction only
        'S':    '#d6aadf',	#Function unknown
        }



def write_plot(features, contigs, plot_path):
    plot_data = Plot_Data(features, contigs, plot_path)
    plot_data.write_config_file()

    for i in features:
        print(i)
    for i in contigs:
        print(i)

    #run Circos
    #run(['perl', '-v'])
    #run(['circos', f'-conf {plot_data.outdir}/circos_config/circos.conf'])
    return