
import logging
import os
import sys
import tempfile
from pathlib import Path

log = logging.getLogger('CONFIG')

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
    threads = args.threads
    log.info('threads=%i', threads)
    verbose = args.verbose
    log.info('verbose=%s', verbose)

    # input / output path configurations
    global db_path, tmp_path, genome_path, output_path
    db_path = None
    if(args.db):
        db_dir = args.db
        log.debug('test parameter db: db_tmp=%s', db_dir)
        try:
            db_tmp_path = Path(db_dir).resolve()
            if(db_tmp_path.is_dir()):
                db_path = db_tmp_path
                log.info('database detected: type=parameter, path=%s', db_path)
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
                log.info('database detected: type=environment, path=%s', db_path)
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
            log.info('database detected: type=base-dir, path=%s', db_path)
        else:
            log.error('unvalid database path: type=base-dir, path=%s', db_tmp_path)
            sys.exit('ERROR: database neither auto-detected nor provided!\nPlease, download the mandatory db and provide it either via the --db parameter, via a BAKTA_DB environment variable or copy it into the Bakta base directory.\nFor further information please read the readme.md')

    if(args.tmp_dir):
        tmp_path = Path(args.tmp_dir)
        if(not tmp_path.exists()):
            log.debug('dedicated temp dir does not exist! tmp-dir=%s', tmp_path)
            sys.exit(f'ERROR: dedicated temp dir ({tmp_path}) does not exist!')
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
            sys.exit(f'ERROR: genome file ({genome_path}) not readable!')
        if(genome_path.stat().st_size == 0):
            log.error('empty genome file! path=%s', genome_path)
            sys.exit(f'ERROR: genome file ({genome_path}) is empty!')
    except:
        log.error('provided genome file not valid! path=%s', args.genome)
        sys.exit(f'ERROR: genome file ({args.genome}) not valid!')
    log.info('genome-path=%s', genome_path)

    log.info('output-path=%s', output_path)

    # input / output configurations
    global min_contig_length, prefix, taxon
    min_contig_length = args.min_contig_length
    log.info('min_contig_length=%s', min_contig_length)
    log.info('prefix=%s', prefix)

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
    
    taxon = f"{genus} {species} {strain}"
    taxon = ' '.join(taxon.replace('None', '').split())
    if(taxon == ''):
        taxon = None

    # annotation configurations
    global prodigal_tf, translation_table, keep_contig_headers, locus, locus_tag, gram, complete, replicons
    prodigal_tf = args.prodigal_tf if args.prodigal_tf != '' else None
    if prodigal_tf is not None:
        try:
            prodigal_tf_path = Path(args.prodigal_tf).resolve()
            if(not os.access(str(prodigal_tf_path), os.R_OK)):
                log.error('prodigal training file not readable! path=%s', prodigal_tf_path)
                sys.exit(f'ERROR: Prodigal training file ({prodigal_tf_path}) not readable!')
            if(genome_path.stat().st_size == 0):
                log.error('empty prodigal training file! path=%s', prodigal_tf_path)
                sys.exit(f'ERROR: Prodigal training file ({prodigal_tf_path}) is empty!')
            prodigal_tf = prodigal_tf_path
        except:
            log.error('provided prodigal training file not valid! path=%s', args.prodigal_tf)
            sys.exit(f'ERROR: Prodigal training file ({args.prodigal_tf}) not valid!')
    log.info('prodigal_tf=%s', prodigal_tf)
    translation_table = args.translation_table if args.translation_table != '' else None
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
                sys.exit(f'ERROR: replicon table file ({replicon_table_path}) not readable!')
            if(genome_path.stat().st_size == 0):
                log.error('empty replicon table file! path=%s', replicon_table_path)
                sys.exit(f'ERROR: replicon table file ({replicon_table_path}) is empty!')
            replicons = replicon_table_path
        except:
            log.error('provided replicon file not valid! path=%s', args.replicons)
            sys.exit(f'ERROR: replicon table file ({args.replicons}) not valid!')
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
