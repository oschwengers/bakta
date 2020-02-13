
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
keep_contig_names = False
complete = False
plasmid = False

# input / output configurations
pretty_json = False
gff3 = False
genbank = False
embl = False


def setup(args):
    """Test environment and build a runtime configuration."""
    # runtime configurations
    global env, bundled_binaries, threads
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

    # genome configurations
    global prefix, locus, locus_tag, genus, species, strain, keep_contig_names, complete, plasmid
    prefix = args.prefix if args.prefix != '' else None
    log.info('prefix=%s', prefix)
    locus = args.locus if args.locus != '' else None
    log.info('locus=%s', locus)
    locus_tag = args.locus_tag if args.locus_tag != '' else None
    log.info('locus-tag=%s', locus_tag)
    genus = args.genus if args.genus != '' else None
    log.info('genus=%s', genus)
    species = args.species if args.species != '' else None
    log.info('species=%s', species)
    strain = args.strain if args.strain != '' else None
    log.info('strain=%s', strain)
    keep_contig_names = args.keep_contig_names
    log.info('keep-contig-names=%s', keep_contig_names)
    complete = args.complete
    log.info('complete=%s', complete)
    plasmid = args.plasmid
    log.info('plasmid=%s', plasmid)

    # input / output configurations
    global pretty_json, gff3, genbank, embl
    pretty_json = args.pretty_json
    log.info('pretty-json=%s', pretty_json)
    gff3 = args.gff3
    log.info('gff3=%s', gff3)
    genbank = args.genbank
    log.info('genbank=%s', genbank)
    embl = args.embl
    log.info('embl=%s', embl)
