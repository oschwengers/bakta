
import argparse
import hashlib
import logging
import multiprocessing as mp
import os
import sys
import subprocess as sp

import bakta
import bakta.config as cfg

log = logging.getLogger('utils')


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog='bakta',
        description='Comprehensive and rapid annotation of bacterial genomes.'
    )
    parser.add_argument('genome', metavar='<genome>', help='(Draft) genome in fasta format')
    parser.add_argument('--db', '-d', action='store', help='Database path (default = <bakta_path>/db)')
    parser.add_argument('--verbose', '-v', action='store_true', help='Print verbose information')
    parser.add_argument('--threads', '-t', action='store', type=int, default=mp.cpu_count(), help='Number of threads to use (default = number of available CPUs)')

    parser.add_argument('--min-contig-length', '-m', action='store', type=int, default=1, dest='min_contig_length', help='Minimum contig size (default = 1)')
    parser.add_argument('--prefix', '-p', action='store', default='', help='Prefix for output files')
    parser.add_argument('--output', '-o', action='store', default=os.getcwd(), help='Output directory (default = current working directory)')
    parser.add_argument('--gff3', action='store_true', help='Write GFF3 annotation file')
    parser.add_argument('--genbank', action='store_true', help='Write GenBank annotation file')
    parser.add_argument('--embl', action='store_true', help='Write EMBL annotation file')

    parser.add_argument('--keep-contig-names', action='store_true', dest='keep_contig_names', help='Keep original contig names')
    parser.add_argument('--locus', action='store', default='', help='Locus prefix')
    parser.add_argument('--locus-tag', action='store', default='', dest='locus_tag', help='Locus tag prefix')
    parser.add_argument('--genus', action='store', default='', help='Genus name')
    parser.add_argument('--species', action='store', default='', help='Species name')
    parser.add_argument('--strain', action='store', default='', help='Strain name')
    parser.add_argument('--plasmid', action='store', default='', help='Plasmid name')
    parser.add_argument('--gram', action='store', default='', choices=['+', '-'], help="Gram type: ''/+/- (default = '')")
    parser.add_argument('--complete', action='store_true', help="Replicons (chromosome/plasmid[s]) are complete")

    parser.add_argument('--version', action='version', version='%(prog)s ' + bakta.__version__)
    parser.add_argument('--citation', action='store_true', help='Print citation')
    return parser.parse_args()


def test_database():
    """Test if database directory exists, is accessible and contains necessary files."""

    if(cfg.db_path is None):
        log.error('database directory not provided nor detected!')
        sys.exit('ERROR: database directory not provided nor detected! Please provide a valid path to the database directory.')

    if(not os.access(str(cfg.db_path), os.R_OK & os.X_OK)):
        log.error('database directory (%s) not readable/accessible!', cfg.db_path)
        sys.exit('ERROR: database directory (%s) not readable/accessible!' % cfg.db_path)

    file_names = [
        'rRNA.i1f',
        'rRNA.i1i',
        'rRNA.i1m',
        'rRNA.i1p',
        'ncRNA.i1f',
        'ncRNA.i1i',
        'ncRNA.i1m',
        'ncRNA.i1p',
        'rfam-go.tsv',
        'bakta.db',
        'psc.dmnd'
    ]

    for file_name in file_names:
        path = cfg.db_path.joinpath(file_name)
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

    # test diamond
    try:
        sp.check_call(
            ['diamond', '--version'],
            stdout=sp.DEVNULL,
            stderr=sp.DEVNULL
        )
    except FileNotFoundError:
        log.exception('diamond not found!')
        sys.exit('ERROR: \'diamond\' not executable!')
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
