import logging
import multiprocessing as mp
import os
import sys
import tempfile

from pathlib import Path


log = logging.getLogger('CONFIG')


# runtime configurations
env = os.environ.copy()
threads = None
verbose = None

# input / output configuration
db_path = None
db_info = None
tmp_path = None
genome_path = None
min_contig_length = None
prefix = None
output_path = None

# organism configuration
genus = None
species = None
strain = None
plasmid = None
taxon = None

# annotation configuration
complete = None
prodigal_tf = None
translation_table = None
keep_contig_headers = None
locus = None
locus_tag = None
gram = None
replicons = None
compliant = None

# workflow configuration
skip_trna = None
skip_tmrna = None
skip_rrna = None
skip_ncrna = None
skip_ncrna_region = None
skip_crispr = None
skip_cds = None
skip_sorf = None
skip_gap = None
skip_ori = None


def setup(args):
    """Test environment and build a runtime configuration."""
    # runtime configurations
    global env, threads, verbose
    threads = args.threads
    if(threads <= 0):
        log.error("wrong argument for 'threads' parameter! threads=%i", threads)
        sys.exit(f"ERROR: wrong argument ({threads}) for 'threads' parameter! Value must be larger than 0.")
    elif(threads > mp.cpu_count()):
        log.error("wrong argument for 'threads' parameter! More threads requested than available: requested=%i, available=%i", threads, mp.cpu_count())
        sys.exit(f"ERROR: wrong argument ({threads}) for 'threads' parameter! More threads requested ({threads}) than available ({mp.cpu_count()}).")
    log.info('threads=%i', threads)
    verbose = args.verbose
    log.info('verbose=%s', verbose)

    # input / output path configurations
    global db_path, db_info, tmp_path, genome_path, min_contig_length, prefix, output_path
    if(args.db):
        db_dir = args.db
        log.debug('test parameter db: db_tmp=%s', db_dir)
        try:
            db_tmp_path = Path(db_dir).resolve()
            if(db_tmp_path.is_dir()):
                db_path = db_tmp_path
                log.info('database: type=parameter, path=%s', db_path)
            else:
                log.error('unvalid database path: type=parameter, path=%s', db_tmp_path)
                raise IOError()
        except:
            sys.exit(f'ERROR: wrong database path! --db={db_dir}')
    elif('BAKTA_DB' in env):
        db_dir = env['BAKTA_DB']
        log.debug('test env db: db_tmp=%s', db_dir)
        try:
            db_tmp_path = Path(db_dir).resolve()
            if(db_tmp_path.is_dir()):
                db_path = db_tmp_path
                log.info('database: type=environment, path=%s', db_path)
            else:
                log.error('unvalid database path: type=environment, path=%s', db_tmp_path)
                raise IOError()
        except:
            sys.exit(f'ERROR: wrong database path! BAKTA_DB={db_dir}')
    else:
        base_dir = Path(__file__).parent.parent
        db_tmp_path = base_dir.joinpath('db')
        log.debug('test base_dir db: db_tmp=%s', db_tmp_path)
        if(db_tmp_path.is_dir()):
            db_path = db_tmp_path
            log.info('database: type=base-dir, path=%s', db_path)
        else:
            log.error('unvalid database path: type=base-dir, path=%s', db_tmp_path)
            sys.exit('ERROR: database neither provided nor auto-detected!\nPlease, download the mandatory db and provide it via either the --db parameter, a BAKTA_DB environment variable or copy it into the Bakta base directory.\nFor further information please read the readme.md')

    if(args.tmp_dir):
        tmp_path = Path(args.tmp_dir)
        if(not tmp_path.exists()):
            log.debug('dedicated temp dir does not exist! tmp-dir=%s', tmp_path)
            sys.exit(f'ERROR: dedicated temporary directory ({tmp_path}) does not exist!')
        else:
            log.info('use dedicated temp dir: path=%s', tmp_path)
            tmp_path = Path(tempfile.mkdtemp(dir=str(tmp_path)))
    else:
        tmp_path = Path(tempfile.mkdtemp())
    log.info('tmp-path=%s', tmp_path)

    try:
        if(args.genome == ''):
            raise ValueError('File path argument must be non-empty')
        genome_path = Path(args.genome).resolve()
        check_readability('genome', genome_path)
        check_content_size('genome', genome_path)
    except:
        log.error('provided genome file not valid! path=%s', args.genome)
        sys.exit(f'ERROR: genome file ({args.genome}) not valid!')
    log.info('genome-path=%s', genome_path)

    # input / output configurations
    min_contig_length = args.min_contig_length
    if(min_contig_length <= 0):
        log.error("wrong argument for 'min-contig-length' parameter! min_contig_length=%s", min_contig_length)
        sys.exit(f"ERROR: wrong argument ({min_contig_length}) for 'min- contig-length' parameter! Value must be larger than 0")
    log.info('min_contig_length=%s', min_contig_length)
    log.info('prefix=%s', prefix)  # set in bakta.py before global logger config
    log.info('output-path=%s', output_path)

    # organism configurations
    global genus, species, strain, plasmid, taxon
    genus = args.genus
    if(genus):
        genus = genus.capitalize()
    log.info('genus=%s', genus)
    species = args.species
    log.info('species=%s', species)
    strain = args.strain
    log.info('strain=%s', strain)
    plasmid = args.plasmid
    log.info('plasmid=%s', plasmid)
    taxon = f'{genus} {species} {strain}'
    taxon = ' '.join(taxon.replace('None', '').split())
    if(taxon == ''):
        taxon = None

    # annotation configurations
    global complete, prodigal_tf, translation_table, keep_contig_headers, locus, locus_tag, gram, replicons, compliant
    complete = args.complete
    log.info('complete=%s', complete)
    prodigal_tf = args.prodigal_tf
    if(prodigal_tf is not None):
        try:
            if(prodigal_tf == ''):
                raise ValueError('File path argument must be non-empty')
            prodigal_tf_path = Path(args.prodigal_tf).resolve()
            check_readability('prodigal training', prodigal_tf_path)
            check_content_size('prodigal training', prodigal_tf_path)
            prodigal_tf = prodigal_tf_path
        except:
            log.error('provided prodigal training file not valid! path=%s', prodigal_tf)
            sys.exit(f'ERROR: Prodigal training file ({prodigal_tf}) not valid!')
    log.info('prodigal_tf=%s', prodigal_tf)
    translation_table = args.translation_table
    log.info('translation_table=%s', translation_table)
    gram = args.gram
    log.info('gram=%s', gram)
    locus = args.locus
    log.info('locus=%s', locus)
    locus_tag = args.locus_tag
    log.info('locus-tag=%s', locus_tag)
    keep_contig_headers = args.keep_contig_headers
    log.info('keep_contig_headers=%s', keep_contig_headers)
    replicons = args.replicons
    if(replicons is not None):
        try:
            if(replicons == ''):
                raise ValueError('File path argument must be non-empty')
            replicon_table_path = Path(args.replicons).resolve()
            check_readability('replicon table', replicon_table_path)
            check_content_size('replicon table', replicon_table_path)
            replicons = replicon_table_path
        except:
            log.error('provided replicon file not valid! path=%s', replicons)
            sys.exit(f'ERROR: replicon table file ({replicons}) not valid!')
    log.info('replicon-table=%s', replicons)
    compliant = args.compliant
    log.info('compliant=%s', compliant)
    if(compliant):
        min_contig_length = 200
        log.info('compliant mode! min_contig_length=%s', min_contig_length)

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


def check_readability(file_name, file_Path):
    if(not os.access(str(file_Path), os.R_OK)):
        log.error('%s file not readable! path=%s', file_name, file_Path)
        sys.exit(f'ERROR: {file_name} file ({file_Path}) not readable!')


def check_content_size(file_name, file_path):
    if(file_path.stat().st_size == 0):
        log.error('empty %s file! path=%s', file_name, file_path)
        sys.exit(f'ERROR: {file_name} file ({file_path}) is empty!')
