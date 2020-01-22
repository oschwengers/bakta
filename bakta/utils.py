
import hashlib
import logging
import os
import sys
import tempfile
import subprocess as sp
from pathlib import Path

log = logging.getLogger('utils')


def log_cmd(logger, args, config):
    logger.info('command line: %s', ' '.join(sys.argv))
    logger.info('configuration: db-path=%s', config['db'])
    logger.info('configuration: bundled binaries=%s', config['bundled-binaries'])
    logger.info('configuration: tmp-path=%s', config['tmp'])
    logger.info('parameters: genome=%s', config['genome_path'])
    logger.info('parameters: output=%s', config['output_path'])
    logger.info('parameters: prefix=%s', args.prefix)
    logger.info('parameters: keep-contig-names=%s', args.keep_contig_names)
    logger.info('parameters: locus=%s', args.locus)
    logger.info('parameters: locus-tag=%s', args.locus_tag)
    logger.info('parameters: genus=%s', args.genus)
    logger.info('parameters: species=%s', args.species)
    logger.info('parameters: strain=%s', args.strain)
    logger.info('parameters: complete=%s', args.complete)
    logger.info('options: gff3=%s', args.gff3)
    logger.info('options: genbank=%s', args.genbank)
    logger.info('options: embl=%s', args.embl)
    logger.info('options: threads=%d', args.threads)
    logger.info('options: pretty-json=%d', args.pretty_json)


def setup_configuration():
    """Test environment and build a runtime configuration."""

    config = {
        'env': os.environ.copy(),
        'tmp': Path(tempfile.mkdtemp()),
        'bundled-binaries': False
    }
    base_dir = Path(__file__).parent.parent
    share_dir = base_dir.joinpath('share')
    log.debug('config: base-dir=%s', base_dir)
    log.debug('config: share-dir=%s', share_dir)
    if(os.access(str(share_dir), os.R_OK & os.X_OK)):
        config['env']["PATH"] = str(share_dir) + ':' + config['env']["PATH"]
        config['bundled-binaries'] = True
    log.debug('config: bundled binaries=%s', config['bundled-binaries'])

    db_path = base_dir.joinpath('db')
    if(os.access(str(db_path), os.R_OK & os.X_OK)):
        config['db'] = db_path
        log.debug('config: use bundled db, db-path=%s', db_path)
    return config


def test_database(config):
    """Test if database directory exists, is accessible and contains necessary files."""

    if('db' not in config):
        log.error('database directory not detected!')
        sys.exit('ERROR: database directory not detected! Please provide a valid path to the database directory.')

    if(not os.access(str(config['db']), os.R_OK & os.X_OK)):
        log.error('database directory (%s) not readable/accessible!', config['db'])
        sys.exit('ERROR: database directory (%s) not readable/accessible!' % config['db'])

    file_names = [
        'rRNA.i1f',
        'rRNA.i1i',
        'rRNA.i1m',
        'rRNA.i1p',
        'ncRNA.i1f',
        'ncRNA.i1i',
        'ncRNA.i1m',
        'ncRNA.i1p',
        'rfam-go.tsv'
        # 'ups.db',
        # 'psc.db',
        # 'psc.inf',
        # 'psc_0.inf',
        # 'psc_0.nam',
        # 'psc_0.off',
        # 'psc_0.seq',
        # 'psc_0.src'
    ]

    for file_name in file_names:
        path = config['db'].joinpath(file_name)
        if(not os.access(str(path), os.R_OK) or not path.is_file()):
            log.error('database file not readable! file=%s', file_name)
            sys.exit('ERROR: database file (%s) not readable!' % file_name)
    return


def test_dependencies():
    """Test the proper installation of necessary 3rd party executables."""

    # test tRNAscan-SE
    try:
        sp.check_call(
            ['tRNAscan-SE', '-h'],
            stdout=sp.DEVNULL,
            stderr=sp.DEVNULL
        )
    except FileNotFoundError:
        log.exception('tRNAscan-SE not found!')
        sys.exit('ERROR: \'tRNAscan-SE\' not executable!')
    except:
        pass

    # test cmsearch
    try:
        sp.check_call(
            ['cmsearch', '-h'],
            stdout=sp.DEVNULL,
            stderr=sp.DEVNULL
        )
    except FileNotFoundError:
        log.exception('cmsearch not found!')
        sys.exit('ERROR: \'cmsearch\' not executable!')
    except:
        pass

    # test prodigal
    try:
        sp.check_call(
            ['prodigal', '-v'],
            stdout=sp.DEVNULL,
            stderr=sp.DEVNULL
        )
    except FileNotFoundError:
        log.exception('prodigal not found!')
        sys.exit('ERROR: \'prodigal\' not executable!')
    except:
        pass

    # test ghostz
    try:
        sp.check_call(
            ['ghostz'],
            stdout=sp.DEVNULL,
            stderr=sp.DEVNULL
        )
    except FileNotFoundError:
        log.exception('ghostz not found!')
        sys.exit('ERROR: \'ghostz\' not executable!')
    except:
        pass

    # test blastn
    try:
        sp.check_call(
            ['blastn', '-version'],
            stdout=sp.DEVNULL,
            stderr=sp.DEVNULL
        )
    except FileNotFoundError:
        log.exception('blastn not found!')
        sys.exit('ERROR: \'blastn\' not executable!')
    except:
        pass


def create_locus_prefix(args, contigs):
    """Create either genus/species or sequence MD5 hex based locus prefix.
    Max locus name length is 37 for GenBank -> 32 + _ + 4 digits"""
    if(args.genus != '' and args.species != ''):
        locus_prefix = args.genus[:1] + args.species[:2]
        return locus_prefix.upper()
    else:
        hash = hashlib.md5()
        for contig in contigs:
            hash.update(str.encode(contig['sequence']))
        return hash.hexdigest()[0:5]


def create_locus_tag_prefix(args, contigs):
    """Create either genus/species or sequence MD5 hex based locus tag prefix."""
    if(args.locus != ''):
        return args.locus
    elif(args.genus != '' and args.species != ''):
        locus_prefix = args.genus[:1] + args.species[:2]
        if(args.strain != ''):
            locus_prefix += args.strain[:2]
        return locus_prefix.upper()
    else:
        hash = hashlib.md5()
        for contig in contigs:
            hash.update(str.encode(contig['sequence']))
        return hash.hexdigest()[0:5]


def calc_aa_hash(seq):
    return hashlib.md5(seq.encode('utf-8')).hexdigest()


def has_annotation(feature, attribute):
    value = feature.get(attribute, None)
    if(value is not None and value != ''):
        return True
    else:
        return False
