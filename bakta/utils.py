import argparse
import collections
import csv
import hashlib
import logging
import os
import platform as pf
import re
import shutil
import sys
import subprocess as sp

from argparse import Namespace
from datetime import datetime
from typing import Dict, Sequence, Tuple
from pathlib import Path

from Bio.Seq import Seq

import bakta
import bakta.constants as bc
import bakta.config as cfg


log = logging.getLogger('UTILS')


def print_version(self):
    return f'v{self.major}.{self.minor}.{self.patch}'


Version = collections.namedtuple('Version', ['major', 'minor', 'patch'], defaults=[0, 0])  # named tuple for version checking, defaults are zero for missing minor/patch
Version.__str__ = print_version


VERSION_MIN_DIGIT = -1
VERSION_MAX_DIGIT = 1000000000000
VERSION_REGEX = re.compile(r'(\d+)\.(\d+)(?:[\.-](\d+))?')  # regex to search for version number in tool output. Takes missing patch version into consideration.

# List of dependencies: tuples for: min version, max version, tool name & command line parameter, dependency check exclusion options
DEPENDENCY_TRNASCAN = (Version(2, 0, 11), Version(VERSION_MAX_DIGIT, VERSION_MAX_DIGIT, VERSION_MAX_DIGIT), VERSION_REGEX, 'tRNAscan-SE', ('tRNAscan-SE', '-h'), ['--skip-trna'])
DEPENDENCY_ARAGORN = (Version(1, 2, 41), Version(VERSION_MAX_DIGIT, VERSION_MAX_DIGIT, VERSION_MAX_DIGIT), VERSION_REGEX, 'Aragorn', ('aragorn', '-h'), ['skip-tmrna'])
DEPENDENCY_CMSCAN = (Version(1, 1, 4), Version(VERSION_MAX_DIGIT, VERSION_MAX_DIGIT, VERSION_MAX_DIGIT), VERSION_REGEX, 'CMscan', ('cmscan', '-h'), ['--skip-rrna', '--skip-ncrna', '--skip-ncrna-region'])
DEPENDENCY_PILERCR = (Version(1, 6), Version(VERSION_MAX_DIGIT, VERSION_MAX_DIGIT, VERSION_MAX_DIGIT), VERSION_REGEX, 'PilerCR', ('pilercr', '-options'), ['--skip-crispr'])
DEPENDENCY_PYRODIGAL = (Version(3, 5, 0), Version(VERSION_MAX_DIGIT, VERSION_MAX_DIGIT, VERSION_MAX_DIGIT), VERSION_REGEX, 'Pyrodigal', (sys.executable, '-c', 'import pyrodigal; print(pyrodigal.__version__)'), ['--skip-cds'])
DEPENDENCY_PYHMMER = (Version(0, 10, 15), Version(VERSION_MAX_DIGIT, VERSION_MAX_DIGIT, VERSION_MAX_DIGIT), VERSION_REGEX, 'Pyhmmer', (sys.executable, '-c', 'import pyhmmer; print(pyhmmer.__version__)'), ['--skip-cds', '--skip-sorf'])
DEPENDENCY_DIAMOND = (Version(2, 1, 10), Version(VERSION_MAX_DIGIT, VERSION_MAX_DIGIT, VERSION_MAX_DIGIT), VERSION_REGEX, 'Diamond', ('diamond', 'help'), ['--skip-cds', '--skip-sorf'])
DEPENDENCY_BLASTN = (Version(2, 14, 0), Version(VERSION_MAX_DIGIT, VERSION_MAX_DIGIT, VERSION_MAX_DIGIT), VERSION_REGEX, 'Blastn', ('blastn', '-version'), ['--skip-ori'])
DEPENDENCY_AMRFINDERPLUS = (Version(4, 0, 3), Version(VERSION_MAX_DIGIT, VERSION_MAX_DIGIT, VERSION_MAX_DIGIT), VERSION_REGEX, 'AMRFinderPlus', ('amrfinder', '--version'), ['--skip-cds'])
DEPENDENCY_PYCIRCLIZE = (Version(1, 7, 0), Version(VERSION_MAX_DIGIT, VERSION_MAX_DIGIT, VERSION_MAX_DIGIT), VERSION_REGEX, 'pyCirclize', (sys.executable, '-c', 'import pycirclize; print(pycirclize.__version__)'), ['--skip-plot'])


def init_parser(sub_command: str=''):
    parser = argparse.ArgumentParser(
        prog=f'bakta{sub_command}',
        description='Rapid & standardized annotation of bacterial genomes, MAGs & plasmids',
        epilog=f'Version: {cfg.version}\nDOI: {bc.BAKTA_DOI}\nURL: github.com/oschwengers/bakta\n\nCitation:\n{bc.BAKTA_CITATION}',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False
    )
    return parser


def parse_arguments():
    parser = init_parser()
    parser.add_argument('genome', metavar='<genome>', help='Genome sequences in (zipped) fasta format')

    arg_group_io = parser.add_argument_group('Input / Output')
    arg_group_io.add_argument('--db', '-d', action='store', default=None, help='Database path (default = <bakta_path>/db). Can also be provided as BAKTA_DB environment variable.')
    arg_group_io.add_argument('--min-contig-length', '-m', action='store', type=int, default=1, dest='min_contig_length', help='Minimum contig/sequence size (default = 1; 200 in compliant mode)')
    arg_group_io.add_argument('--prefix', '-p', action='store', default=None, help='Prefix for output files')
    arg_group_io.add_argument('--output', '-o', action='store', default=os.getcwd(), help='Output directory (default = current working directory)')
    arg_group_io.add_argument('--force', '-f', action='store_true', help='Force overwriting existing output folder (except for current working directory)')

    arg_group_organism = parser.add_argument_group('Organism')
    arg_group_organism.add_argument('--genus', action='store', default=None, help='Genus name')
    arg_group_organism.add_argument('--species', action='store', default=None, help='Species name')
    arg_group_organism.add_argument('--strain', action='store', default=None, help='Strain name')
    arg_group_organism.add_argument('--plasmid', action='store', default=None, help='Plasmid name')

    arg_group_annotation = parser.add_argument_group('Annotation')
    arg_group_annotation.add_argument('--complete', action='store_true', help='All sequences are complete replicons (chromosome/plasmid[s])')
    arg_group_annotation.add_argument('--prodigal-tf', action='store', default=None, dest='prodigal_tf', help='Path to existing Prodigal training file to use for CDS prediction')
    arg_group_annotation.add_argument('--translation-table', action='store', type=int, default=11, choices=[11, 4, 25], dest='translation_table', help='Translation table: 11/4/25 (default = 11)')
    arg_group_annotation.add_argument('--gram', action='store', default=bc.GRAM_UNKNOWN, choices=[bc.GRAM_POSITIVE, bc.GRAM_NEGATIVE, bc.GRAM_UNKNOWN], help=f'Gram type for signal peptide predictions: {bc.GRAM_POSITIVE}/{bc.GRAM_NEGATIVE}/{bc.GRAM_UNKNOWN} (default = {bc.GRAM_UNKNOWN})')
    arg_group_annotation.add_argument('--locus', action='store', default=None, help="Locus prefix (default = 'contig')")
    arg_group_annotation.add_argument('--locus-tag', action='store', default=None, dest='locus_tag', help='Locus tag prefix (default = autogenerated)')
    arg_group_annotation.add_argument('--locus-tag-increment', action='store', type=int, default=1, choices=[1, 5, 10], dest='locus_tag_increment', help='Locus tag increment: 1/5/10 (default = 1)')
    arg_group_annotation.add_argument('--keep-contig-headers', action='store_true', dest='keep_contig_headers', help='Keep original contig/sequence headers')
    arg_group_annotation.add_argument('--compliant', action='store_true', help='Force Genbank/ENA/DDJB compliance')
    arg_group_annotation.add_argument('--replicons', '-r', action='store', default=None, dest='replicons', help='Replicon information table (tsv/csv)')
    arg_group_annotation.add_argument('--regions', action='store', default=None, help='Path to pre-annotated regions in GFF3 or Genbank format (regions only, no functional annotations).')
    arg_group_annotation.add_argument('--proteins', action='store', default=None, dest='proteins', help='Fasta file of trusted protein sequences for CDS annotation')
    arg_group_annotation.add_argument('--hmms', action='store', default=None, dest='hmms', help='HMM file of trusted hidden markov models in HMMER format for CDS annotation')
    arg_group_annotation.add_argument('--meta', action='store_true', help='Run in metagenome mode. This only affects CDS prediction.')

    arg_group_workflow = parser.add_argument_group('Workflow')
    arg_group_workflow.add_argument('--skip-trna', action='store_true', dest='skip_trna', help='Skip tRNA detection & annotation')
    arg_group_workflow.add_argument('--skip-tmrna', action='store_true', dest='skip_tmrna', help='Skip tmRNA detection & annotation')
    arg_group_workflow.add_argument('--skip-rrna', action='store_true', dest='skip_rrna', help='Skip rRNA detection & annotation')
    arg_group_workflow.add_argument('--skip-ncrna', action='store_true', dest='skip_ncrna', help='Skip ncRNA detection & annotation')
    arg_group_workflow.add_argument('--skip-ncrna-region', action='store_true', dest='skip_ncrna_region', help='Skip ncRNA region detection & annotation')
    arg_group_workflow.add_argument('--skip-crispr', action='store_true', dest='skip_crispr', help='Skip CRISPR array detection & annotation')
    arg_group_workflow.add_argument('--skip-cds', action='store_true', dest='skip_cds', help='Skip CDS detection & annotation')
    arg_group_workflow.add_argument('--skip-pseudo', action='store_true', dest='skip_pseudo', help='Skip pseudogene detection & annotation')
    arg_group_workflow.add_argument('--skip-sorf', action='store_true', dest='skip_sorf', help='Skip sORF detection & annotation')
    arg_group_workflow.add_argument('--skip-gap', action='store_true', dest='skip_gap', help='Skip gap detection & annotation')
    arg_group_workflow.add_argument('--skip-ori', action='store_true', dest='skip_ori', help='Skip oriC/oriT detection & annotation')
    arg_group_workflow.add_argument('--skip-filter', action='store_true', dest='skip_filter', help='Skip feature overlap filters')
    arg_group_workflow.add_argument('--skip-plot', action='store_true', dest='skip_plot', help='Skip generation of circular genome plots')

    arg_group_general = parser.add_argument_group('General')
    arg_group_general.add_argument('--help', '-h', action='help', help='Show this help message and exit')
    arg_group_general.add_argument('--verbose', '-v', action='store_true', help='Print verbose information')
    arg_group_general.add_argument('--debug', action='store_true', help='Run Bakta in debug mode. Temp data will not be removed.')
    arg_group_general.add_argument('--threads', '-t', action='store', type=int, default=0, help='Number of threads to use (default = number of available CPUs)')
    arg_group_general.add_argument('--tmp-dir', action='store', default=None, dest='tmp_dir', help='Location for temporary files (default = system dependent auto detection)')
    arg_group_general.add_argument('--version', action='version', version=f'%(prog)s {cfg.version}')
    return parser.parse_args()


def setup_logger(output_path: Path, prefix: str, args: Namespace):
    logging.basicConfig(
        filename=str(output_path.joinpath(f'{prefix}.log')),
        filemode='w',
        format='%(asctime)s.%(msecs)03d - %(levelname)s - %(name)s - %(message)s',
        datefmt='%H:%M:%S',
        level=logging.DEBUG if args.debug else logging.INFO
    )
    log.info('version=%s', cfg.version)
    log.info('developer: Oliver Schwengers, github.com/oschwengers')
    log.info('command: %s', ' '.join(sys.argv))
    log.info('local time: %s', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    log.info('machine: type=%s, cores=%s', pf.processor(), os.cpu_count())
    log.info('system: type=%s, release=%s', pf.system(), pf.release())
    log.info('python: version=%s, implementation=%s', pf.python_version(), pf.python_implementation())


def cleanup(log, tmp_path: Path):
    shutil.rmtree(str(tmp_path))  # remove tmp dir
    log.info('removed tmp dir: %s', tmp_path)


def read_tool_output(dependency):
    """Method for reading tool version with regex. Input: regex expression, tool command. Retursn: version number."""
    version_regex = dependency[2]
    tool = dependency[3]
    command = dependency[4]
    skip_options = dependency[5]
    try:
        tool_output = str(sp.check_output(command, stderr=sp.STDOUT))  # stderr must be added in case the tool output is not piped into stdout
        log.debug('tool output: tool=%s, output=%s', tool, tool_output)
    except FileNotFoundError:
        log.error('dependency not found! tool=%s', tool)
        sys.exit(f"ERROR: {tool} not found or not executable! Please make sure {tool} is installed and executable or skip requiring workflow steps via via '{' '.join(skip_options)}'.")
    except sp.CalledProcessError:
        log.error('dependency check failed! tool=%s', tool)
        sys.exit(f"ERROR: {tool} could not be executed! Please make sure {tool} is installed and executable or skip requiring workflow steps via via '{' '.join(skip_options)}'.")
    version_match = version_regex.search(tool_output)
    try:
        if version_match is None:
            log.error('no dependency version detected! No regex hit in dependency output: regex=%s, command=%s', version_regex, command)
            sys.exit(f'ERROR: Could not detect/read {tool} version! No regex hit in dependency output.')
        major = version_match.group(1)
        minor = version_match.group(2)
        patch = version_match.group(3)
        if major is None:
            log.error('no dependency version detected! Major hit is None in dependency output: regex=%s, command=%s', version_regex, command)
            sys.exit(f'ERROR: Could not detect/read {tool} version! Major hit is None in dependency output.')
        elif minor is None:
            version_output = Version(int(major))
        elif patch is None:
            version_output = Version(int(major), int(minor))
        else:
            version_output = Version(int(major), int(minor), int(patch))
        return version_output
    except:
        log.error('no dependency version detected! Error while reading dependency output: regex=%s, command=%s', version_regex, command)
        sys.exit(f'ERROR: Could not detect/read {tool} version! Error while reading dependency output.')


def check_version(tool, min: int, max:int ) -> bool:
    """Method for checking tool versions with required version. Input: tool version, minimum and maximum version. Returns: boolean value for positive or negative check."""
    if tool.major < min.major or tool.major > max.major:
        return False
    else:
        if tool.major == min.major or tool.major == max.major:
            if tool.minor < min.minor and tool.major == min.major:
                return False
            elif tool.minor > max.minor and tool.major == max.major:
                return False
            else:
                if tool.minor == min.minor or tool.minor == max.minor:
                    if tool.patch < min.patch and tool.minor == min.minor:
                        return False
                    elif tool.patch > max.patch and tool.minor == max.minor and tool.major == max.major:
                        return False
                    else:
                        return True
                else:
                    return True
        else:
            return True


def test_dependency(dependency):
    """Test the proper installation of the required 3rd party executable."""
    version = read_tool_output(dependency)
    check_result = check_version(version, dependency[0], dependency[1])
    tool_name = dependency[3]
    if (not check_result):
        log.error('wrong dependency version for %s: installed=%s, minimum=%s', tool_name, version, dependency[0])
        sys.exit(f'ERROR: Wrong {tool_name} version installed. Please, either install {tool_name} version {dependency[0]} or use {dependency[4]}!')
    else:
        tool_path = shutil.which(dependency[4][0])
        log.info('dependency: tool=%s, version=%s, path=%s', tool_name, version, tool_path)


def test_dependencies():
    """Test the proper installation of all required 3rd party executables."""
    if(cfg.skip_trna is not None and cfg.skip_trna is False):
        test_dependency(DEPENDENCY_TRNASCAN)

    if(cfg.skip_tmrna is not None and cfg.skip_tmrna is False):
        test_dependency(DEPENDENCY_ARAGORN)

    if((cfg.skip_rrna is not None and cfg.skip_rrna is False) or (cfg.skip_ncrna is not None and cfg.skip_ncrna is False) or (cfg.skip_ncrna_region is not None and cfg.skip_ncrna_region is False)):
        test_dependency(DEPENDENCY_CMSCAN)

    if(cfg.skip_crispr is not None and cfg.skip_crispr is False):
        test_dependency(DEPENDENCY_PILERCR)

    if(cfg.skip_cds is not None and cfg.skip_cds is False):
        test_dependency(DEPENDENCY_PYRODIGAL)
        test_dependency(DEPENDENCY_AMRFINDERPLUS)

        # test if AMRFinderPlus db is installed
        amrfinderplus_db_path = cfg.db_path.joinpath('amrfinderplus-db')
        amrfinderplus_db_latest_path = amrfinderplus_db_path.joinpath('latest')
        process = sp.run(
            [
                'amrfinder',
                '--debug',
                '--database', str(amrfinderplus_db_latest_path)
            ], capture_output=True)
        if('No valid AMRFinder database found' in process.stderr.decode()):
            log.error('AMRFinderPlus database not installed')
            sys.exit(f"ERROR: AMRFinderPlus database not installed! Please, install AMRFinderPlus's internal database by executing: 'amrfinder_update --database {amrfinderplus_db_path}'. This must be done only once.")

    if((cfg.skip_cds is not None and cfg.skip_cds is False) or (cfg.skip_sorf is not None and cfg.skip_sorf is False)):
        test_dependency(DEPENDENCY_PYHMMER)
        test_dependency(DEPENDENCY_DIAMOND)

    if(cfg.skip_ori is not None and cfg.skip_ori is False):
        test_dependency(DEPENDENCY_BLASTN)
    
    if(cfg.skip_plot is not None and cfg.skip_plot is False):
        test_dependency(DEPENDENCY_PYCIRCLIZE)


def create_locus_tag_prefix(sequences: Sequence[dict], length: int=6) -> str:
    """Create either genus/species or sequence MD5 hex based locus tag prefix."""
    hash = hashlib.md5()
    for seq in sequences:
        hash.update(str.encode(seq['nt']))
    hexdigest = hash.hexdigest().upper()
    locus_prefix_chars = []
    i = 0
    while i < length:
        c = hexdigest[i]
        if(c >= '0' and c <= '9'):
            c = chr(ord('F') + int(c) + 1)
        locus_prefix_chars.append(c)
        i += 1
    locus_prefix = ''.join(locus_prefix_chars)
    log.info('generated sequence tag prefix: prefix=%s, length=%i, MD5=%s', locus_prefix, length, hexdigest)
    return locus_prefix


def calc_aa_hash(seq: str) -> Tuple[bytes, str]:
    aa_hash = hashlib.md5(seq.encode('utf-8'))
    return (aa_hash.digest(), aa_hash.hexdigest())


def has_annotation(feature: dict, attribute: str) -> bool:
    value = feature.get(attribute, None)
    if(value is not None and value != ''):
        return True
    else:
        return False


def calc_genome_stats(data: dict):
    genome_size = data['stats']['size']
    log.info('genome-size=%i', genome_size)

    # N50
    gc_sum = 0
    n_sum = 0
    for seq in data['sequences']:
        nt = seq['nt']
        gc_sum += nt.count('G') + nt.count('C')
        n_sum += nt.count('N')
    gc_ratio = gc_sum / (genome_size - n_sum)
    data['stats']['gc'] = gc_ratio
    log.info('GC=%0.3f', gc_ratio)

    n_ratio = n_sum / genome_size
    data['stats']['n_ratio'] = n_ratio
    log.info('N=%0.3f', n_ratio)

    sequence_length_sum = 0
    for seq in sorted(data['sequences'], key=lambda x: x['length'], reverse=True):
        nt_length = len(seq['nt'])
        sequence_length_sum += nt_length
        if(sequence_length_sum >= genome_size * 0.5):
            if 'n50' not in data['stats']:
                data['stats']['n50'] = nt_length
                log.info('N50=%i', nt_length)
        if(sequence_length_sum >= genome_size * 0.9):
            if 'n90' not in data['stats']:
                data['stats']['n90'] = nt_length
                log.info('N90=%i', nt_length)
    

    sequence_by_id = {seq['id']: seq for seq in data['sequences']}
    coding_nts = 0
    for feat in data['features']:
        if(feat.get('edge', False)):
            sequence_length = sequence_by_id[feat['sequence']]['length']
            coding_nts += feat['stop'] + (sequence_length - feat['start'] + 1)  # feature coding nucleotides
        else:
            coding_nts += feat['stop'] - feat['start'] + 1  # feature coding nucleotides
    coding_ratio = coding_nts / (genome_size - n_sum)
    data['stats']['coding_ratio'] = coding_ratio
    log.info('coding-ratio=%0.3f', coding_ratio)


def parse_replicon_table(replicon_table_path: Path) -> Dict[str, dict]:
    replicons = {}
    try:
        with replicon_table_path.open() as fh:
            dialect = csv.Sniffer().sniff(fh.read(1024), delimiters=",\t")
            fh.seek(0)
            reader = csv.reader(fh, dialect)
            for row in reader:
                (original_locus_id, new_locus_id, replicon_type, topology, name) = row
                original_locus_id = original_locus_id.strip()
                if(' ' in original_locus_id):
                    original_locus_id = original_locus_id.split(' ')[0]  # remove description
                new_locus_id = new_locus_id.strip()
                if(' ' in new_locus_id):
                    new_locus_id = new_locus_id.split(' ')[0]  # remove description

                # TODO: add locus id checks
                if(new_locus_id == '' or new_locus_id == '-'):
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


def qc_sequences(sequences: Sequence[dict], replicons: Dict[str, dict]) -> Tuple[Sequence[dict], bool]:
    valid_sequences = []
    sequence_counter = 1
    sequence_prefix = cfg.locus if cfg.locus else 'contig'
    complete_genome = True
    plasmid_number = 1
    sequence_ids = set()
    for seq in sequences:
        if(seq['length'] >= cfg.min_sequence_length):
            sequence_id_generated = f'{sequence_prefix}_{sequence_counter}'
            seq['simple_id'] = sequence_id_generated
            sequence_counter += 1

            sequence_description = seq['description'].lower()
            if(cfg.complete):
                seq['complete'] = True
                seq['topology'] = bc.TOPOLOGY_CIRCULAR
            elif('circular=true' in sequence_description):  # detection of Unicycler circularized sequences
                seq['complete'] = True
                seq['topology'] = bc.TOPOLOGY_CIRCULAR
                log.debug('qc: detected Unicycler circular topology via description: id=%s, description=%s', seq['id'], seq['description'])
            elif('complete' in sequence_description and 'complete=false' not in sequence_description):  # detection of public/described sequences
                seq['complete'] = True
                seq['topology'] = bc.TOPOLOGY_CIRCULAR
                log.debug('qc: detected complete replicon via description: id=%s, description=%s', seq['id'], seq['description'])
            
            if('chromosome' in sequence_description):
                seq['type'] = bc.REPLICON_CHROMOSOME
                log.debug('qc: detected chromosome replicon type via description: id=%s, description=%s', seq['id'], seq['description'])
            elif('plasmid' in sequence_description):
                seq['type'] = bc.REPLICON_PLASMID
                log.debug('qc: detected plasmid replicon type via description: id=%s, description=%s', seq['id'], seq['description'])

            sequence_desc = []
            if(cfg.keep_sequence_headers):
                if(seq['id'] in sequence_ids):
                    log.error('Fasta import: duplicated seq id! seq-id=%s', seq['id'])
                    sys.exit(f"ERROR: Detected duplicated sequence id! Sequence ID ({seq['id']}) occures multiple times!")
                else:
                    sequence_ids.add(seq['id'])
            else:
                seq['orig_id'] = seq['id']
                seq['id'] = sequence_id_generated
                seq['orig_description'] = seq['description']
                if(cfg.genus is not None or cfg.species is not None):
                    organism = ' '.join([t for t in [cfg.genus, cfg.species] if t is not None])
                    sequence_desc.append(f"[organism={organism}]")
                if(cfg.strain):
                    sequence_desc.append(f'[strain={cfg.strain}]')
                sequence_desc.append(f'[gcode={cfg.translation_table}]')

            if(seq['complete'] and seq['topology'] == bc.TOPOLOGY_CIRCULAR):  # detection of chromosomes/plasmids via sequence length thresholds
                if(seq['length'] >= bc.REPLICON_LENGTH_THRESHOLD_CHROMOSOME):
                    seq['type'] = bc.REPLICON_CHROMOSOME
                    log.debug('qc: detected replicon type via length: id=%s, type=%s, length=%i, description=%s', seq['id'], seq['type'], seq['length'], seq['description'])
                elif(seq['length'] < bc.REPLICON_LENGTH_THRESHOLD_PLASMID):
                    seq['type'] = bc.REPLICON_PLASMID
                    log.debug('qc: detected replicon type via length: id=%s, type=%s, length=%i, description=%s', seq['id'], seq['type'], seq['length'], seq['description'])
            valid_sequences.append(seq)

            if(len(sequences) == 1 and cfg.plasmid is not None):  # use plasmid mode
                seq['type'] = bc.REPLICON_PLASMID
                seq['topology'] = bc.TOPOLOGY_CIRCULAR
                seq['name'] = cfg.plasmid
            elif(replicons):  # use user provided replicon table
                sequence_id = seq['orig_id'] if 'orig_id' in seq else seq['id']
                replicon = replicons.get(sequence_id, None)
                if(replicon):
                    seq['type'] = replicon['replicon_type']
                    seq['topology'] = replicon['topology']
                    seq['complete'] = replicon['replicon_type'] != bc.REPLICON_CONTIG
                    if(replicon['name']):
                        seq['name'] = replicon['name']
                    if(not cfg.keep_sequence_headers):
                        seq['id'] = replicon['new_locus_id'] if replicon['new_locus_id'] else seq['simple_id']
                    seq.pop('simple_id')
            
            if(not cfg.keep_sequence_headers):
                if(seq['complete']):
                    sequence_desc.append('[completeness=complete]')
                if(seq['topology'] != bc.REPLICON_CONTIG):
                    sequence_desc.append(f"[topology={seq['topology']}]")
                if(seq['type'] == bc.REPLICON_CHROMOSOME):
                    sequence_desc.append('[location=chromosome]')
                elif(seq['type'] == bc.REPLICON_PLASMID):
                    if(not seq.get('name', None)):
                        seq['name'] = f'unnamed{plasmid_number}'
                        plasmid_number += 1
                    sequence_desc.append(f"[plasmid-name={seq['name']}]")
                seq['description'] = ' '.join(list(dict.fromkeys(sequence_desc)))  # remove duplicates remaining order

            if(seq['type'] == bc.REPLICON_CONTIG):
                complete_genome = False

            if(cfg.compliant):  # check INSDC compliance
                if(len(seq['id']) > 25):  # max 25 characters
                    log.error('INSDC compliance: seq id larger than 25! seq-id=%s', seq['id'])
                    sys.exit(f"ERROR: INSDC compliance failed! Sequence ID ({seq['id']}) larger than 25 characers!")
                if(bc.RE_INSDC_ID.fullmatch(seq['id']) is None):  # invalid characters
                    log.error('INSDC compliance: seq id contains invalid characters! seq-id=%s', seq['id'])
                    sys.exit(f"ERROR: INSDC compliance failed! Sequence ID ({seq['id']}) contains invalid characters!")

            log.info(
                "qc: revised sequence: id=%s, orig-id=%s, type=%s, complete=%s, topology=%s, name=%s, description='%s', orig-description='%s'",
                seq['id'], seq.get('orig_id', ''), seq['type'], seq['complete'], seq['topology'], seq.get('name', ''), seq['description'], seq.get('orig_description', '')
            )
    return valid_sequences, complete_genome


def extract_feature_sequence(feature: dict, sequence: dict) -> str:
    if(feature.get('edge', False)):
        nt = sequence['nt'][feature['start']-1:] + sequence['nt'][:feature['stop']]
    else:
        nt = sequence['nt'][feature['start']-1:feature['stop']]
    if(feature['strand'] == bc.STRAND_REVERSE):
        nt = str(Seq(nt).reverse_complement())
    return nt
