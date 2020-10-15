
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
        description='Comprehensive and rapid annotation of bacterial genomes.',
        add_help=False
    )
    parser.add_argument('genome', metavar='<genome>', help='(Draft) genome in fasta format')

    arg_group_io = parser.add_argument_group('Input / Output')
    arg_group_io.add_argument('--db', '-d', action='store', help='Database path (default = <bakta_path>/db)')
    arg_group_io.add_argument('--min-contig-length', '-m', action='store', type=int, default=1, dest='min_contig_length', help='Minimum contig size (default = 1)')
    arg_group_io.add_argument('--prefix', '-p', action='store', default='', help='Prefix for output files')
    arg_group_io.add_argument('--output', '-o', action='store', default=os.getcwd(), help='Output directory (default = current working directory)')
    arg_group_io.add_argument('--tsv', action='store_true', help='Write TSV annotation file')
    arg_group_io.add_argument('--gff3', action='store_true', help='Write GFF3 annotation file')
    arg_group_io.add_argument('--genbank', action='store_true', help='Write GenBank annotation file')
    arg_group_io.add_argument('--embl', action='store_true', help='Write EMBL annotation file')
    arg_group_io.add_argument('--faa', action='store_true', help='Write translated CDS sequences as fasta file')

    arg_group_organism = parser.add_argument_group('Organism')
    arg_group_organism.add_argument('--genus', action='store', default='', help='Genus name')
    arg_group_organism.add_argument('--species', action='store', default='', help='Species name')
    arg_group_organism.add_argument('--strain', action='store', default='', help='Strain name')
    arg_group_organism.add_argument('--plasmid', action='store', default='', help='Plasmid name')
    
    arg_group_annotation = parser.add_argument_group('Annotation')
    arg_group_annotation.add_argument('--prodigal-tf', action='store', dest='prodigal_tf', help='Path to existing Prodigal training file to use for CDS prediction')
    arg_group_annotation.add_argument('--translation-table', action='store', type=int, default=11, dest='translation_table', help='Translation table to use (default = 11)')
    arg_group_annotation.add_argument('--keep-contig-names', action='store_true', dest='keep_contig_names', help='Keep original contig names')
    arg_group_annotation.add_argument('--locus', action='store', default='', help='Locus prefix')
    arg_group_annotation.add_argument('--locus-tag', action='store', default='', dest='locus_tag', help='Locus tag prefix')
    arg_group_annotation.add_argument('--gram', action='store', default='?', choices=['+', '-', '?'], help="Gram type: +/-/? (default = '?')")
    arg_group_annotation.add_argument('--complete', action='store_true', help="Replicons (chromosome/plasmid[s]) are complete")

    arg_group_workflow = parser.add_argument_group('Workflow')
    arg_group_workflow.add_argument('--skip-trna', action='store_true', dest='skip_trna', help="Skip tRNA detection & annotation")
    arg_group_workflow.add_argument('--skip-tmrna', action='store_true', dest='skip_tmrna', help="Skip tmRNA detection & annotation")
    arg_group_workflow.add_argument('--skip-rrna', action='store_true', dest='skip_rrna', help="Skip rRNA detection & annotation")
    arg_group_workflow.add_argument('--skip-ncrna', action='store_true', dest='skip_ncrna', help="Skip ncRNA detection & annotation")
    arg_group_workflow.add_argument('--skip-ncrna-region', action='store_true', dest='skip_ncrna_region', help="Skip ncRNA region detection & annotation")
    arg_group_workflow.add_argument('--skip-crispr', action='store_true', dest='skip_crispr', help="Skip CRISPR array detection & annotation")
    arg_group_workflow.add_argument('--skip-cds', action='store_true', dest='skip_cds', help="Skip CDS detection & annotation")
    arg_group_workflow.add_argument('--skip-sorf', action='store_true', dest='skip_sorf', help="Skip sORF detection & annotation")
    arg_group_workflow.add_argument('--skip-gap', action='store_true', dest='skip_gap', help="Skip gap detection & annotation")

    arg_group_general = parser.add_argument_group('General')
    arg_group_general.add_argument('--help', '-h', action='help', help='Show this help message and exit')
    arg_group_general.add_argument('--verbose', '-v', action='store_true', help='Print verbose information')
    arg_group_general.add_argument('--threads', '-t', action='store', type=int, default=mp.cpu_count(), help='Number of threads to use (default = number of available CPUs)')
    arg_group_general.add_argument('--tmp-dir', action='store', default=None, dest='tmp_dir', help='Location for temporary files (default = system dependent auto detection)')
    arg_group_general.add_argument('--version', action='version', version='%(prog)s ' + bakta.__version__)
    arg_group_general.add_argument('--citation', action='store_true', help='Print citation')
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
        'ncRNA-genes.i1f',
        'ncRNA-genes.i1i',
        'ncRNA-genes.i1m',
        'ncRNA-genes.i1p',
        'ncRNA-regions.i1f',
        'ncRNA-regions.i1i',
        'ncRNA-regions.i1m',
        'ncRNA-regions.i1p',
        'rfam-go.tsv',
        'antifam.h3f',
        'antifam.h3i',
        'antifam.h3m',
        'antifam.h3p',
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

    # test pilercr
    try:
        sp.check_call(
            ['pilercr', '-options'],
            stdout=sp.DEVNULL,
            stderr=sp.DEVNULL
        )
    except FileNotFoundError:
        log.exception('pilercr not found!')
        sys.exit('ERROR: \'pilercr\' not executable!')
    except:
        pass


def create_locus_prefix(contigs):
    """Create either genus/species or sequence MD5 hex based locus prefix.
    Max locus name length is 37 for GenBank -> 32 + _ + 4 digits"""
    if(cfg.genus != '' and cfg.species != ''):
        locus_prefix = cfg.genus[:1] + cfg.species[:2]
        return locus_prefix.upper()
    else:
        hash = hashlib.md5()
        for contig in contigs:
            hash.update(str.encode(contig['sequence']))
        return hash.hexdigest()[0:5]


def create_locus_tag_prefix(contigs):
    """Create either genus/species or sequence MD5 hex based locus tag prefix."""
    hash = hashlib.md5()
    for contig in contigs:
        hash.update(str.encode(contig['sequence']))
    hexdigest = hash.hexdigest().upper()
    locus_prefix = []
    i = 0
    while i < 6:
        c = hexdigest[i]
        if(c >= '0' and c <= '9'):
            c = chr(ord('F') + int(c) + 1) 
        locus_prefix.append(c)
        i += 1
    return ''.join(locus_prefix)


def calc_aa_hash(seq):
    aa_hash = hashlib.md5(seq.encode('utf-8'))
    return (aa_hash.digest(), aa_hash.hexdigest())


def has_annotation(feature, attribute):
    value = feature.get(attribute, None)
    if(value is not None and value != ''):
        return True
    else:
        return False
