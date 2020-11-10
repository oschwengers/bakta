
import argparse
import hashlib
import logging
import multiprocessing as mp
import os
import re
import sys
import subprocess as sp

import bakta
import bakta.constants as bc
import bakta.config as cfg

log = logging.getLogger('UTILS')


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

    arg_group_organism = parser.add_argument_group('Organism')
    arg_group_organism.add_argument('--genus', action='store', default='', help='Genus name')
    arg_group_organism.add_argument('--species', action='store', default='', help='Species name')
    arg_group_organism.add_argument('--strain', action='store', default='', help='Strain name')
    arg_group_organism.add_argument('--plasmid', action='store', default='', help='Plasmid name')
    
    arg_group_annotation = parser.add_argument_group('Annotation')
    arg_group_annotation.add_argument('--prodigal-tf', action='store', dest='prodigal_tf', help='Path to existing Prodigal training file to use for CDS prediction')
    arg_group_annotation.add_argument('--translation-table', action='store', type=int, default=11, choices=[11, 4], dest='translation_table', help='Translation table to use: 11/4 (default = 11)')
    arg_group_annotation.add_argument('--complete', action='store_true', help="Replicons (chromosome/plasmid[s]) are complete")
    arg_group_annotation.add_argument('--gram', action='store', default='?', choices=['+', '-', '?'], help="Gram type: +/-/? (default = '?')")
    arg_group_annotation.add_argument('--locus', action='store', default='', help="Locus prefix (instead of 'contig')")
    arg_group_annotation.add_argument('--locus-tag', action='store', default='', dest='locus_tag', help='Locus tag prefix')
    arg_group_annotation.add_argument('--keep-contig-headers', action='store_true', dest='keep_contig_headers', help='Keep original contig headers')
    arg_group_annotation.add_argument('--replicons', '-r', action='store', dest='replicons', help="Replicon information table (TSV)")

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
    arg_group_workflow.add_argument('--skip-ori', action='store_true', dest='skip_ori', help="Skip oriC/oriT detection & annotation")

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
        sys.exit(f'ERROR: database directory ({cfg.db_path}) not readable/accessible!')

    file_names = [
        'antifam.h3f',
        'antifam.h3i',
        'antifam.h3m',
        'antifam.h3p',
        'bakta.db',
        'ncRNA-genes.i1f',
        'ncRNA-genes.i1i',
        'ncRNA-genes.i1m',
        'ncRNA-genes.i1p',
        'ncRNA-regions.i1f',
        'ncRNA-regions.i1i',
        'ncRNA-regions.i1m',
        'ncRNA-regions.i1p',
        'oric.fna',
        'orit.fna',
        'psc.dmnd',
        'rfam-go.tsv',
        'rRNA.i1f',
        'rRNA.i1i',
        'rRNA.i1m',
        'rRNA.i1p',
        'sorf.dmnd'
    ]

    for file_name in file_names:
        path = cfg.db_path.joinpath(file_name)
        if(not os.access(str(path), os.R_OK) or not path.is_file()):
            log.error('database file not readable! file=%s', file_name)
            sys.exit(f'ERROR: database file ({file_name}) not readable!')
    return


def test_dependencies():
    """Test the proper installation of necessary 3rd party executables."""

    # test tRNAscan-SE
    if(cfg.skip_trna is False):
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

    # test Aragorn
    if(cfg.skip_tmrna is False):
        try:
            sp.check_call(
                ['aragorn', '-h'],
                stdout=sp.DEVNULL,
                stderr=sp.DEVNULL
            )
        except FileNotFoundError:
            log.exception('aragorn not found!')
            sys.exit('ERROR: \'aragorn\' not executable!')
        except:
            pass

    # test cmscan
    if(cfg.skip_rrna is False or cfg.skip_ncrna is False or cfg.skip_ncrna_region is False):
        try:
            sp.check_call(
                ['cmscan', '-h'],
                stdout=sp.DEVNULL,
                stderr=sp.DEVNULL
            )
        except FileNotFoundError:
            log.exception('cmscan not found!')
            sys.exit('ERROR: \'cmscan\' not executable!')
        except:
            pass

    # test prodigal
    if(cfg.skip_cds is False):
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

    # test hmmsearch
    if(cfg.skip_cds is False or cfg.skip_sorf is False):
        try:
            sp.check_call(
                ['hmmsearch', '-h'],
                stdout=sp.DEVNULL,
                stderr=sp.DEVNULL
            )
        except FileNotFoundError:
            log.exception('hmmsearch not found!')
            sys.exit('ERROR: \'hmmsearch\' not executable!')
        except:
            pass

    # test diamond
    if(cfg.skip_cds is False or cfg.skip_sorf is False):
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
    if(cfg.skip_crispr is False):
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


def calc_genome_stats(genome, features):
    
    genome_size = genome['size']
    log.info('genome-size=%i', genome_size)

    # N50
    gc_sum = 0
    n_sum = 0
    n50 = 0
    contig_length_sum = 0
    for contig in genome['contigs']:
        seq = contig['sequence']
        gc_sum += seq.count('G') + seq.count('C')
        n_sum += seq.count('N')
        contig_length = len(seq)
        contig_length_sum += contig_length
        if(contig_length_sum >= genome_size / 2):
            n50 = contig_length
            break
    genome['n50'] = n50
    log.info('N50=%i', n50)

    gc_ratio = gc_sum / (genome_size - n_sum)
    genome['gc'] = gc_ratio
    log.info('GC=%0.3f', gc_ratio)
    
    n_ratio = n_sum / genome_size
    genome['n_ratio'] = n_ratio
    log.info('N=%0.3f', n_ratio)

    coding_nts = 0
    for feat in features:
        coding_nts += feat['stop'] - feat['start'] + 1  # feature coding nucleotides
    coding_ratio = coding_nts / (genome_size - n_sum)
    genome['coding_ratio'] = coding_ratio
    log.info('coding-ratio=%0.3f', coding_ratio)
    
    return {
        'gc': gc_ratio,
        'n_ratio': n_ratio,
        'n50': n50,
        'coding_ratio': coding_ratio
    }


def parse_replicon_table(replicon_table_path):
    replicons = {}
    try:
        with replicon_table_path.open() as fh:
            for line in fh:
                (original_locus_id, new_locus_id, replicon_type, topology, name) = line.strip().split('\t')
                #ToDO add locus id checks
                if(new_locus_id == '' or new_locus_id == ''):
                    new_locus_id = None
                replicon_type = replicon_type.lower()
                if(replicon_type == 'c' or 'chrom' in replicon_type):
                    replicon_type = bc.REPLICON_CHROMOSOME
                elif(replicon_type == 'p' or 'plasmid' in replicon_type):
                    replicon_type = bc.REPLICON_PLASMID
                else:
                    replicon_type = bc.REPLICON_CONTIG
                topology = topology.lower()
                if(topology == 'c' or 'circ' in topology):
                    topology = bc.TOPOLOGY_CIRCULAR
                else:
                    topology = bc.TOPOLOGY_LINEAR
                if(replicon_type == bc.REPLICON_CONTIG):
                    topology == bc.TOPOLOGY_LINEAR
                if(name == '' or name == '-'):
                    name = None
                replicon = {
                    'original_locus_id': original_locus_id,
                    'new_locus_id': new_locus_id,
                    'replicon_type': replicon_type,
                    'topology': topology,
                    'name': name
                }
                log.debug(
                    'parse replicon info: orig-id=%s, new-id=%s, type=%s, topology=%s, name=%s',
                    replicon['original_locus_id'], replicon['new_locus_id'], replicon['replicon_type'], replicon['topology'], replicon['name']
                )
                replicons[original_locus_id] = replicon
    except:
        log.error('wrong replicon table format!')
        sys.exit('ERROR: wrong replicon table file format!')
    return replicons


def qc_contigs(contigs, replicons):
    valid_contigs = []
    contig_counter = 1
    contig_prefix = cfg.locus if cfg.locus else 'contig'
    organism_definition = []
    if(cfg.genus):
        organism_definition.append(cfg.genus)
    if(cfg.species):
        organism_definition.append(cfg.species)
    if(len(organism_definition) > 0):
        organism_definition = ' '.join(organism_definition)
        organism_definition = f"[organism={organism_definition}]"
    else:
        organism_definition = None
    
    complete_genome = True
    for contig in contigs:
        if(contig['length'] >= cfg.min_contig_length):
            contig_name = f'{contig_prefix}_{contig_counter}'
            contig['simple_id'] = contig_name
            contig_counter += 1
            if(not cfg.keep_contig_headers):
                contig['orig_id'] = contig['id']
                contig['id'] = contig_name
                contig['orig_desc'] = contig['desc']
                contig_desc = []
                if(organism_definition):
                    contig_desc.append(organism_definition)
                if(cfg.strain):
                    contig_desc.append(f'[strain={cfg.strain}]')
                if(cfg.complete):
                    contig_desc.append('[completeness=complete]')
                contig['desc'] = ' '.join(contig_desc)
            if(cfg.complete):
                contig['complete'] = True
                contig['topology'] = bc.TOPOLOGY_CIRCULAR
            valid_contigs.append(contig)
            
            if(replicons):
                contig_id = contig['orig_id'] if 'orig_id' in contig else contig['id']
                replicon = replicons.get(contig_id, None)
                if(replicon):
                    contig['type'] = replicon['replicon_type']
                    contig['topology'] = replicon['topology']
                    contig['name'] = replicon['name']
                    if(replicon['replicon_type'] != bc.REPLICON_CONTIG):
                        contig['complete'] = True
                    if(not cfg.keep_contig_headers):
                        contig['id'] = replicon['new_locus_id'] if replicon['new_locus_id'] else contig['simple_id']
                        if(replicon['replicon_type'] != bc.REPLICON_CONTIG and 'completeness' not in contig['desc']):
                            contig['desc'] += ' [completeness=complete]'
                        contig['desc'] += f" [topology={replicon['topology']}]"
                        if(replicon['replicon_type'] == bc.REPLICON_PLASMID and replicon['name']):
                            contig['desc'] += f" [plasmid-name={replicon['name']}]"
                    contig.pop('simple_id')

            if(contig['type'] == bc.REPLICON_CONTIG):
                complete_genome = False

            log.info(
                "revised contig: id=%s, orig-id=%s, type=%s, topology=%s, name=%s, desc='%s', orig-desc='%s'",
                contig['id'], contig.get('orig_id', ''), contig['type'], contig['topology'], contig.get('name', ''), contig['desc'], contig.get('orig_desc', '')
            )
    return valid_contigs, complete_genome
