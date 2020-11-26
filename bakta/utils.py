
import argparse
import hashlib
import logging
import multiprocessing as mp
import os
import sys
import subprocess as sp

import bakta
import bakta.constants as bc
import bakta.config as cfg

log = logging.getLogger('UTILS')


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog='bakta',
        description='Rapid & standardized annotation of bacterial genomes & plasmids.',
        add_help=False
    )
    parser.add_argument('genome', metavar='<genome>', help='(Draft) genome in fasta format')

    arg_group_io = parser.add_argument_group('Input / Output')
    arg_group_io.add_argument('--db', '-d', action='store', default=None, help='Database path (default = <bakta_path>/db)')
    arg_group_io.add_argument('--min-contig-length', '-m', action='store', type=int, default=1, dest='min_contig_length', help='Minimum contig size (default = 1)')
    arg_group_io.add_argument('--prefix', '-p', action='store', default=None, help='Prefix for output files')
    arg_group_io.add_argument('--output', '-o', action='store', default=os.getcwd(), help='Output directory (default = current working directory)')

    arg_group_organism = parser.add_argument_group('Organism')
    arg_group_organism.add_argument('--genus', action='store', default=None, help='Genus name')
    arg_group_organism.add_argument('--species', action='store', default=None, help='Species name')
    arg_group_organism.add_argument('--strain', action='store', default=None, help='Strain name')
    arg_group_organism.add_argument('--plasmid', action='store', default=None, help='Plasmid name')
    
    arg_group_annotation = parser.add_argument_group('Annotation')
    arg_group_annotation.add_argument('--complete', action='store_true', help="All sequences are complete replicons (chromosome/plasmid[s])")
    arg_group_annotation.add_argument('--prodigal-tf', action='store', default=None, dest='prodigal_tf', help='Path to existing Prodigal training file to use for CDS prediction')
    arg_group_annotation.add_argument('--translation-table', action='store', type=int, default=11, choices=[11, 4], dest='translation_table', help='Translation table to use: 11/4 (default = 11)')
    arg_group_annotation.add_argument('--gram', action='store', default='?', choices=['+', '-', '?'], help="Gram type: +/-/? (default = '?')")
    arg_group_annotation.add_argument('--locus', action='store', default=None, help="Locus prefix (instead of 'contig')")
    arg_group_annotation.add_argument('--locus-tag', action='store', default=None, dest='locus_tag', help='Locus tag prefix')
    arg_group_annotation.add_argument('--keep-contig-headers', action='store_true', dest='keep_contig_headers', help='Keep original contig headers')
    arg_group_annotation.add_argument('--replicons', '-r', action='store', default=None, dest='replicons', help="Replicon information table (TSV)")

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
            sys.exit('ERROR: tRNAscan-SE not executable! Please either install tRNAscan-SE (>=v2.0.6) or skip the annotation of tRNAs via \'--skip-trna\'.')
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
            sys.exit('ERROR: aragorn not executable! Please either install Aragorn (>=1.2.38) or skip the annotation of tmRNAs via \'--skip-tmrna\'.')
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
            sys.exit('ERROR: \'cmscan\' not executable! Please either install Infernal (>=v1.1.2) or skip the annotation of rRNAs, ncRNAs & ncRNA regions via \'--skip-rrna --skip-ncrna --skip-ncrna-region\'.')
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
            sys.exit('ERROR: \'prodigal\' not executable! Please either install Prodigal (>=v2.6.3) or skip the annotation of CDSs via \'--skip-cds\'.')
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
            sys.exit('ERROR: \'hmmsearch\' not executable! Please either install HMMER (>=3.3.1) or skip the annotation of CDSs & sORFs via \'--skip-cds --skip-sorf\'.')
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
            sys.exit('ERROR: \'diamond\' not executable! Please either install Diamond (>=v2.0.4) or skip the annotation of CDSs & sORFs via \'--skip-cds --skip-sorf\'.')
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
            sys.exit('ERROR: pilercr not executable! Please either install PILER-CR or skip the annotation of CRISPRs via \'--skip-crispr\'.')
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

    contigs_by_id = {c['id']: c for c in genome['contigs']}
    coding_nts = 0
    for feat in features:
        if(feat.get('edge', False)):
            sequence_length = contigs_by_id[feat['contig']]['length']
            coding_nts += feat['stop'] + (sequence_length - feat['start'] + 1)  # feature coding nucleotides
        else:
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
                # TODO: add locus id checks
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
                log.info(
                    'replicon info: orig-id=%s, new-id=%s, type=%s, topology=%s, name=%s',
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
    organism_definition = f"[organism={cfg.taxon}]" if cfg.taxon else None
    
    complete_genome = True
    plasmid_number = 1
    for contig in contigs:
        if(contig['length'] >= cfg.min_contig_length):
            contig_id_generated = f'{contig_prefix}_{contig_counter}'
            contig['simple_id'] = contig_id_generated
            contig_counter += 1

            if('circular=true' in contig['description'].lower()):  # detection of Unicycler circularized sequences
                contig['complete'] = True
                contig['topology'] = bc.TOPOLOGY_CIRCULAR
                log.debug('qc: detected Unicycler circular topology via description: id=%s, description=%s', contig['id'], contig['description'])
            if('complete' in contig['description'].lower()):
                contig['complete'] = True
                contig['topology'] = bc.TOPOLOGY_CIRCULAR
                log.debug('qc: detected complete replicon via description: id=%s, description=%s', contig['id'], contig['description'])
            if('chromosome' in contig['description'].lower()):
                contig['type'] = bc.REPLICON_CHROMOSOME
                log.debug('qc: detected chromosome replicon type via description: id=%s, description=%s', contig['id'], contig['description'])
            elif('plasmid' in contig['description'].lower()):
                contig['type'] = bc.REPLICON_PLASMID
                log.debug('qc: detected plasmid replicon type via description: id=%s, description=%s', contig['id'], contig['description'])
            
            if(not cfg.keep_contig_headers):
                contig['orig_id'] = contig['id']
                contig['id'] = contig_id_generated
                contig['orig_description'] = contig['description']
                contig_desc = []
                if(organism_definition):
                    contig_desc.append(organism_definition)
                if(cfg.strain):
                    contig_desc.append(f'[strain={cfg.strain}]')
                if(cfg.complete or contig['complete']):
                    contig_desc.append('[completeness=complete]')
                    if(contig['topology'] != bc.REPLICON_CONTIG):
                        contig_desc.append(f"[topology={contig['topology']}]")
                contig['description'] = ' '.join(contig_desc)
            if(cfg.complete):
                contig['complete'] = True
                contig['topology'] = bc.TOPOLOGY_CIRCULAR
                if(contig['length'] >= bc.REPLICON_LENGTH_THRESHOLD_CHROMOSOME):
                    contig['type'] = bc.REPLICON_CHROMOSOME
                elif(contig['length'] < bc.REPLICON_LENGTH_THRESHOLD_PLASMID):
                    contig['type'] = bc.REPLICON_PLASMID
            
            valid_contigs.append(contig)
            
            if(replicons):  # use user provided replicon table
                contig_id = contig['orig_id'] if 'orig_id' in contig else contig['id']
                replicon = replicons.get(contig_id, None)
                if(replicon):
                    contig['type'] = replicon['replicon_type']
                    contig['topology'] = replicon['topology']
                    if(replicon['name']):
                        contig['name'] = replicon['name']
                    if(replicon['replicon_type'] != bc.REPLICON_CONTIG):
                        contig['complete'] = True
                    if(not cfg.keep_contig_headers):
                        contig['id'] = replicon['new_locus_id'] if replicon['new_locus_id'] else contig['simple_id']
                        if(replicon['replicon_type'] != bc.REPLICON_CONTIG and 'completeness' not in contig['description']):
                            contig['description'] += ' [completeness=complete]'
                        if('topology' not in contig['description']):
                            contig['description'] += f" [topology={replicon['topology']}]"
                        if(replicon['replicon_type'] == bc.REPLICON_PLASMID and replicon['name']):
                            contig['description'] += f" [plasmid-name={replicon['name']}]"
                    contig.pop('simple_id')

            if(contig['complete'] and contig['type'] == bc.REPLICON_PLASMID and not contig.get('name', None)):
                contig['name'] = f'unnamed{plasmid_number}'
                plasmid_number += 1
            
            if(contig['type'] == bc.REPLICON_CONTIG):
                complete_genome = False
            
            log.info(
                "qc: revised sequence: id=%s, orig-id=%s, type=%s, complete=%s, topology=%s, name=%s, description='%s', orig-description='%s'",
                contig['id'], contig.get('orig_id', ''), contig['type'], contig['complete'], contig['topology'], contig.get('name', ''), contig['description'], contig.get('orig_description', '')
            )
    return valid_contigs, complete_genome
