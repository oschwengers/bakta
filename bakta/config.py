
import logging
import os
import sys
import tempfile
from pathlib import Path

log = logging.getLogger('config')

# runtime configurations
env = os.environ.copy()
bundled_binaries = False
threads = 1

# input / output path configurations
db_path = None
tmp_path = None
genome_path = None
output_path = None

# genome configurations
prefix = None
locus = None
locus_tag = None
genus = None
species = None
strain = None
keep_contig_headers = False
complete = False
plasmid = False

# input / output configurations
gff3 = False
genbank = False
embl = False


def setup(args):
    """Test environment and build a runtime configuration."""
    # runtime configurations
    global env, bundled_binaries, threads, verbose
    base_dir = Path(__file__).parent.parent
    share_dir = base_dir.joinpath('share')
    log.debug('base-dir=%s', base_dir)
    log.debug('share-dir=%s', share_dir)
    if(share_dir.is_dir() and os.access(str(share_dir), os.R_OK & os.X_OK)):
        env["PATH"] = str(share_dir) + ':' + env["PATH"]
        bundled_binaries = True
        log.debug('found bundled binaries')
    log.info('bundled-binaries=%s', bundled_binaries)
    threads = args.threads
    log.info('threads=%i', threads)
    verbose = args.verbose
    log.info('verbose=%s', verbose)

    # input / output path configurations
    global db_path, tmp_path, genome_path, output_path
    if(args.db):
        try:
            db_path = Path(args.db).resolve()
            if(not db_path.is_dir()):
                raise IOError('DB path is not a directory!')
        except:
            log.error('database directory (%s) not readable/accessible!', args.db)
            sys.exit('ERROR: database directory (%s) not readable/accessible!' % args.db)
    else:
        tmp = base_dir.joinpath('db')
        if(os.access(str(tmp), os.R_OK & os.X_OK)):
            log.debug('found bundled db=%s', tmp)
            db_path = tmp
    log.info('db-path=%s', db_path)

    if(args.tmp_dir):
        tmp_path = Path(args.tmp_dir)
        if(not tmp_path.exists()):
            log.debug('dedicated temp dir does not exist! tmp-dir=%s', tmp_path)
            sys.exit('ERROR: dedicated temp dir (%s) does not exist!' % tmp_path)
        else:
            log.info('use dedicated temp dir: path=%s', tmp_path)
            tmp_path = Path(tempfile.mkdtemp(dir=str(tmp_path)))
    else:
        tmp_path = Path(tempfile.mkdtemp())
    log.info('tmp-path=%s', tmp_path)

    try:
        genome_path = Path(args.genome).resolve()
        if(not os.access(str(genome_path), os.R_OK)):
            log.error('genome file not readable! path=%s', genome_path)
            sys.exit('ERROR: genome file (%s) not readable!' % genome_path)
        if(genome_path.stat().st_size == 0):
            log.error('empty genome file! path=%s', genome_path)
            sys.exit('ERROR: genome file (%s) is empty!' % genome_path)
    except:
        log.error('provided genome file not valid! path=%s', args.genome)
        sys.exit('ERROR: genome file (%s) not valid!' % args.genome)
    log.info('genome-path=%s', genome_path)

    try:
        output_path = Path(args.output) if args.output else Path.cwd()
        if(not output_path.exists()):
            output_path.mkdir(parents=True, exist_ok=True)
            log.info('created output dir: path=%s', output_path)
        else:
            log.debug('use existing output directory')
        output_path = output_path.resolve()
    except:
        log.error('could not set/create output directory! path=%s', args.output)
        sys.exit('ERROR: could not set/create output directory (%s)!' % args.output)
    log.info('output-path=%s', output_path)

    # input / output configurations
    global min_contig_length, output, prefix, tsv, gff3, genbank, embl, fna, faa
    min_contig_length = args.min_contig_length
    log.info('min_contig_length=%s', min_contig_length)
    prefix = args.prefix if args.prefix != '' else None
    log.info('prefix=%s', prefix)
    output = args.output
    log.info('output=%s', output)
    tsv = args.tsv
    log.info('tsv=%s', tsv)
    gff3 = args.gff3
    log.info('gff3=%s', gff3)
    genbank = args.genbank
    log.info('genbank=%s', genbank)
    embl = args.embl
    log.info('embl=%s', embl)
    faa = args.faa
    log.info('faa=%s', faa)
    fna = args.fna
    log.info('fna=%s', fna)

    # organism configurations
    global genus, species, strain, plasmid
    genus = args.genus if args.genus != '' else None
    log.info('genus=%s', genus)
    species = args.species if args.species != '' else None
    log.info('species=%s', species)
    strain = args.strain if args.strain != '' else None
    log.info('strain=%s', strain)
    plasmid = args.plasmid
    log.info('plasmid=%s', plasmid)

    # annotation configurations
    global prodigal_tf, translation_table, keep_contig_headers, locus, locus_tag, gram, complete, replicons
    prodigal_tf = args.prodigal_tf if args.prodigal_tf != '' else None
    if prodigal_tf is not None:
        try:
            prodigal_tf_path = Path(args.prodigal_tf).resolve()
            if(not os.access(str(prodigal_tf_path), os.R_OK)):
                log.error('prodigal training file not readable! path=%s', prodigal_tf_path)
                sys.exit('ERROR: Prodigal training file (%s) not readable!' % prodigal_tf_path)
            if(genome_path.stat().st_size == 0):
                log.error('empty prodigal training file! path=%s', prodigal_tf_path)
                sys.exit('ERROR: Prodigal training file (%s) is empty!' % prodigal_tf_path)
            prodigal_tf = prodigal_tf_path
        except:
            log.error('provided prodigal training file not valid! path=%s', args.prodigal_tf)
            sys.exit('ERROR: Prodigal training file (%s) not valid!' % args.prodigal_tf)
    log.info('prodigal_tf=%s', prodigal_tf)
    translation_table = args.translation_table
    log.info('translation_table=%s', translation_table)
    keep_contig_headers = args.keep_contig_headers
    log.info('keep_contig_headers=%s', keep_contig_headers)
    locus = args.locus if args.locus != '' else None
    log.info('locus=%s', locus)
    locus_tag = args.locus_tag if args.locus_tag != '' else None
    log.info('locus-tag=%s', locus_tag)
    gram = args.gram
    log.info('gram=%s', gram)
    complete = args.complete
    log.info('complete=%s', complete)
    replicons = args.replicons if args.replicons != '' else None
    if replicons is not None:
        try:
            replicon_table_path = Path(args.replicons).resolve()
            if(not os.access(str(replicon_table_path), os.R_OK)):
                log.error('replicon table not readable! path=%s', replicon_table_path)
                sys.exit('ERROR: replicon table file (%s) not readable!' % replicon_table_path)
            if(genome_path.stat().st_size == 0):
                log.error('empty replicon table file! path=%s', replicon_table_path)
                sys.exit('ERROR: replicon table file (%s) is empty!' % replicon_table_path)
            replicons = replicon_table_path
        except:
            log.error('provided replicon file not valid! path=%s', args.replicons)
            sys.exit('ERROR: replicon table file (%s) not valid!' % args.replicons)
    log.info('replicon-table=%s', replicons)
    
    # workflow configurations
    global skip_trna, skip_tmrna, skip_rrna, skip_ncrna, skip_ncrna_region, skip_crispr, skip_cds, skip_sorf, skip_gap, skip_ori
    skip_trna = args.skip_trna
    log.info('skip-tRNA=%s', skip_trna)
    skip_tmrna = args.skip_tmrna
    log.info('skip-tmRNA=%s', skip_tmrna)
    skip_rrna = args.skip_rrna
    log.info('skip-rRNA=%s', skip_rrna)
    skip_ncrna = args.skip_ncrna
    log.info('skip-ncRNA=%s', skip_ncrna)
    skip_ncrna_region = args.skip_ncrna_region
    log.info('skip-ncRNA-region=%s', skip_ncrna_region)
    skip_crispr = args.skip_crispr
    log.info('skip-CRISPR=%s', skip_crispr)
    skip_cds = args.skip_cds
    log.info('skip-CDS=%s', skip_cds)
    skip_sorf = args.skip_sorf
    log.info('skip-sORF=%s', skip_sorf)
    skip_gap = args.skip_gap
    log.info('skip-gap=%s', skip_gap)
    skip_ori = args.skip_ori
    log.info('skip-ori=%s', skip_ori)
