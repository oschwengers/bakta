import argparse
import collections
import csv
import hashlib
import logging
import os
import shutil
import platform as pf
import sys
import subprocess as sp
import re

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
DEPENDENCY_TRNASCAN = (Version(2, 0, 8), Version(VERSION_MAX_DIGIT, VERSION_MAX_DIGIT, VERSION_MAX_DIGIT), VERSION_REGEX, ('tRNAscan-SE', '-h'), ['--skip-trna'])
DEPENDENCY_ARAGORN = (Version(1, 2, 38), Version(VERSION_MAX_DIGIT, VERSION_MAX_DIGIT, VERSION_MAX_DIGIT), VERSION_REGEX, ('aragorn', '-h'), ['skip-tmrna'])
DEPENDENCY_CMSCAN = (Version(1, 1, 4), Version(VERSION_MAX_DIGIT, VERSION_MAX_DIGIT, VERSION_MAX_DIGIT), VERSION_REGEX, ('cmscan', '-h'), ['--skip-rrna', '--skip-ncrna', '--skip-ncrna-region'])
DEPENDENCY_PILERCR = (Version(1, 6), Version(VERSION_MAX_DIGIT, VERSION_MAX_DIGIT, VERSION_MAX_DIGIT), VERSION_REGEX, ('pilercr', '-options'), ['--skip-crispr'])
DEPENDENCY_HMMSEARCH = (Version(3, 3, 2), Version(VERSION_MAX_DIGIT, VERSION_MAX_DIGIT, VERSION_MAX_DIGIT), VERSION_REGEX, ('hmmsearch', '-h'), ['--skip-cds', '--skip-sorf'])
DEPENDENCY_DIAMOND = (Version(2, 0, 14), Version(VERSION_MAX_DIGIT, VERSION_MAX_DIGIT, VERSION_MAX_DIGIT), VERSION_REGEX, ('diamond', 'help'), ['--skip-cds', '--skip-sorf'])
DEPENDENCY_DEEPSIG = (Version(1, 2, 5), Version(VERSION_MAX_DIGIT, VERSION_MAX_DIGIT, VERSION_MAX_DIGIT), VERSION_REGEX, ('deepsig', '--version'), ['--gram ?'])
DEPENDENCY_BLASTN = (Version(2, 12, 0), Version(VERSION_MAX_DIGIT, VERSION_MAX_DIGIT, VERSION_MAX_DIGIT), VERSION_REGEX, ('blastn', '-version'), ['--skip-ori'])
DEPENDENCY_AMRFINDERPLUS = (Version(3, 10, 23), Version(VERSION_MAX_DIGIT, VERSION_MAX_DIGIT, VERSION_MAX_DIGIT), VERSION_REGEX, ('amrfinder', '--version'), ['--skip-cds'])
DEPENDENCY_CIRCOS = (Version(0, 69, 8), Version(VERSION_MAX_DIGIT, VERSION_MAX_DIGIT, VERSION_MAX_DIGIT), VERSION_REGEX, ('circos', '--version'), ['--skip-plot'])


def init_parser(sub_command: str=''):
    parser = argparse.ArgumentParser(
        prog=f'bakta{sub_command}',
        description='Rapid & standardized annotation of bacterial genomes, MAGs & plasmids',
        epilog=f'Version: {bakta.__version__}\nDOI: {bc.BAKTA_DOI}\nURL: github.com/oschwengers/bakta\n\nCitation:\n{bc.BAKTA_CITATION}',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False
    )
    return parser


def parse_arguments():
    parser = init_parser()
    parser.add_argument('genome', metavar='<genome>', help='Genome sequences in (zipped) fasta format')

    arg_group_io = parser.add_argument_group('Input / Output')
    arg_group_io.add_argument('--db', '-d', action='store', default=None, help='Database path (default = <bakta_path>/db). Can also be provided as BAKTA_DB environment variable.')
    arg_group_io.add_argument('--min-contig-length', '-m', action='store', type=int, default=1, dest='min_contig_length', help='Minimum contig size (default = 1; 200 in compliant mode)')
    arg_group_io.add_argument('--prefix', '-p', action='store', default=None, help='Prefix for output files')
    arg_group_io.add_argument('--output', '-o', action='store', default=os.getcwd(), help='Output directory (default = current working directory)')

    arg_group_organism = parser.add_argument_group('Organism')
    arg_group_organism.add_argument('--genus', action='store', default=None, help='Genus name')
    arg_group_organism.add_argument('--species', action='store', default=None, help='Species name')
    arg_group_organism.add_argument('--strain', action='store', default=None, help='Strain name')
    arg_group_organism.add_argument('--plasmid', action='store', default=None, help='Plasmid name')

    arg_group_annotation = parser.add_argument_group('Annotation')
    arg_group_annotation.add_argument('--complete', action='store_true', help='All sequences are complete replicons (chromosome/plasmid[s])')
    arg_group_annotation.add_argument('--prodigal-tf', action='store', default=None, dest='prodigal_tf', help='Path to existing Prodigal training file to use for CDS prediction')
    arg_group_annotation.add_argument('--translation-table', action='store', type=int, default=11, choices=[11, 4], dest='translation_table', help='Translation table: 11/4 (default = 11)')
    arg_group_annotation.add_argument('--gram', action='store', default=bc.GRAM_UNKNOWN, choices=[bc.GRAM_POSITIVE, bc.GRAM_NEGATIVE, bc.GRAM_UNKNOWN], help=f'Gram type for signal peptide predictions: {bc.GRAM_POSITIVE}/{bc.GRAM_NEGATIVE}/{bc.GRAM_UNKNOWN} (default = {bc.GRAM_UNKNOWN})')
    arg_group_annotation.add_argument('--locus', action='store', default=None, help="Locus prefix (default = 'contig')")
    arg_group_annotation.add_argument('--locus-tag', action='store', default=None, dest='locus_tag', help='Locus tag prefix (default = autogenerated)')
    arg_group_annotation.add_argument('--keep-contig-headers', action='store_true', dest='keep_contig_headers', help='Keep original contig headers')
    arg_group_annotation.add_argument('--replicons', '-r', action='store', default=None, dest='replicons', help='Replicon information table (tsv/csv)')
    arg_group_annotation.add_argument('--compliant', action='store_true', help='Force Genbank/ENA/DDJB compliance')
    arg_group_annotation.add_argument('--proteins', action='store', default=None, dest='proteins', help='Fasta file of trusted protein sequences for CDS annotation')

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
    arg_group_workflow.add_argument('--skip-plot', action='store_true', dest='skip_plot', help='Skip generation of circular genome plots')

    arg_group_general = parser.add_argument_group('General')
    arg_group_general.add_argument('--help', '-h', action='help', help='Show this help message and exit')
    arg_group_general.add_argument('--verbose', '-v', action='store_true', help='Print verbose information')
    arg_group_general.add_argument('--debug', action='store_true', help='Run Bakta in debug mode. Temp data will not be removed.')
    arg_group_general.add_argument('--threads', '-t', action='store', type=int, default=0, help='Number of threads to use (default = number of available CPUs)')
    arg_group_general.add_argument('--tmp-dir', action='store', default=None, dest='tmp_dir', help='Location for temporary files (default = system dependent auto detection)')
    arg_group_general.add_argument('--version', action='version', version=f'%(prog)s {bakta.__version__}')
    return parser.parse_args()


def setup_logger(output_path: Path, prefix: str, args: Namespace):
    logging.basicConfig(
        filename=str(output_path.joinpath(f'{prefix}.log')),
        filemode='w',
        format='%(asctime)s.%(msecs)03d - %(levelname)s - %(name)s - %(message)s',
        datefmt='%H:%M:%S',
        level=logging.DEBUG if args.debug else logging.INFO
    )
    log.info('version=%s', bakta.__version__)
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
    command = dependency[3]
    skip_options = dependency[4]
    try:
        tool_output = str(sp.check_output(command, stderr=sp.STDOUT))  # stderr must be added in case the tool output is not piped into stdout
    except FileNotFoundError:
        log.error('dependency not found! tool=%s', command[0])
        sys.exit(f"ERROR: {command[0]} not found or not executable! Please make sure {command[0]} is installed and executable or skip requiring workflow steps via via '{' '.join(skip_options)}'.")
    except sp.CalledProcessError:
        log.error('dependency check failed! tool=%s', command[0])
        sys.exit(f"ERROR: {command[0]} could not be executed! Please make sure {command[0]} is installed and executable or skip requiring workflow steps via via '{' '.join(skip_options)}'.")
    version_match = version_regex.search(tool_output)
    try:
        if version_match is None:
            log.error('no dependency version detected! no regex hit in dependency output: regex=%s, command=%s', version_regex, command)
            sys.exit('ERROR: Could not detect/read %s version!', command[0])
        major = version_match.group(1)
        minor = version_match.group(2)
        patch = version_match.group(3)
        if major is None:
            log.error('no dependency version detected! no regex hit in dependency output: regex=%s, command=%s', version_regex, command)
            sys.exit('ERROR: Could not detect/read %s version!', command[0])
        elif minor is None:
            version_output = Version(int(major))
        elif patch is None:
            version_output = Version(int(major), int(minor))
        else:
            version_output = Version(int(major), int(minor), int(patch))
        return version_output
    except:
        log.error('no dependency version detected! no regex hit in dependency output: regex=%s, command=%s', version_regex, command)
        sys.exit('ERROR: Could not detect/read %s version!', command[0])


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
    if (not check_result):
        log.error('wrong dependency version for %s: installed=%s, minimum=%s', dependency[3][0], version, dependency[0])
        sys.exit(f'ERROR: Wrong {dependency[3][0]} version installed. Please, either install {dependency[3][0]} version {dependency[0]} or use {dependency[4]}!')
    else:
        log.info('dependency: tool=%s, version=%s', dependency[3][0], version)


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
            sys.exit(f"ERROR: AMRFinderPlus database not installed! Please, install AMRFinderPlus's internal database by executing: 'amrfinder_update --database ${amrfinderplus_db_path}'. This must be done only once.")

    if((cfg.skip_cds is not None and cfg.skip_cds is False) or (cfg.skip_sorf is not None and cfg.skip_sorf is False)):
        test_dependency(DEPENDENCY_HMMSEARCH)
        test_dependency(DEPENDENCY_DIAMOND)
        if(cfg.gram is not None and cfg.gram != '?'):
            test_dependency(DEPENDENCY_DEEPSIG)

    if(cfg.skip_ori is not None and cfg.skip_ori is False):
        test_dependency(DEPENDENCY_BLASTN)
    
    if(cfg.skip_plot is not None and cfg.skip_plot is False):
        test_dependency(DEPENDENCY_CIRCOS)


def create_locus_tag_prefix(contigs: Sequence[dict], length: int=6) -> str:
    """Create either genus/species or sequence MD5 hex based locus tag prefix."""
    hash = hashlib.md5()
    for contig in contigs:
        hash.update(str.encode(contig['sequence']))
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


def calc_genome_stats(genome: dict, features: Sequence[dict]):
    genome_size = genome['size']
    log.info('genome-size=%i', genome_size)

    # N50
    gc_sum = 0
    n_sum = 0
    for contig in genome['contigs']:
        seq = contig['sequence']
        gc_sum += seq.count('G') + seq.count('C')
        n_sum += seq.count('N')
    gc_ratio = gc_sum / (genome_size - n_sum)
    genome['gc'] = gc_ratio
    log.info('GC=%0.3f', gc_ratio)

    n_ratio = n_sum / genome_size
    genome['n_ratio'] = n_ratio
    log.info('N=%0.3f', n_ratio)

    n50 = 0
    contig_length_sum = 0
    for contig in sorted(genome['contigs'], key=lambda x: x['length'], reverse=True):
        contig_length = len(contig['sequence'])
        contig_length_sum += contig_length
        if(contig_length_sum >= genome_size / 2):
            n50 = contig_length
            break
    genome['n50'] = n50
    log.info('N50=%i', n50)

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


def parse_replicon_table(replicon_table_path: Path) -> Dict[str, dict]:
    replicons = {}
    try:
        with replicon_table_path.open() as fh:
            dialect = csv.Sniffer().sniff(fh.read(1024), delimiters=",\t")
            fh.seek(0)
            reader = csv.reader(fh, dialect)
            for row in reader:
                (original_locus_id, new_locus_id, replicon_type, topology, name) = row
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


def qc_contigs(contigs: Sequence[dict], replicons: Dict[str, dict]) -> Tuple[Sequence[dict], bool]:
    valid_contigs = []
    contig_counter = 1
    contig_prefix = cfg.locus if cfg.locus else 'contig'
    organism_definition = f"[organism={cfg.taxon}]" if cfg.taxon else None

    complete_genome = True
    plasmid_number = 1
    contig_ids = set()
    for contig in contigs:
        if(contig['length'] >= cfg.min_contig_length):
            contig_id_generated = f'{contig_prefix}_{contig_counter}'
            contig['simple_id'] = contig_id_generated
            contig_counter += 1

            contig_description = contig['description'].lower()
            if(cfg.complete):
                contig['complete'] = True
                contig['topology'] = bc.TOPOLOGY_CIRCULAR
            elif('circular=true' in contig_description):  # detection of Unicycler circularized sequences
                contig['complete'] = True
                contig['topology'] = bc.TOPOLOGY_CIRCULAR
                log.debug('qc: detected Unicycler circular topology via description: id=%s, description=%s', contig['id'], contig['description'])
            elif('complete' in contig_description):  # detection of public/described sequences
                contig['complete'] = True
                contig['topology'] = bc.TOPOLOGY_CIRCULAR
                log.debug('qc: detected complete replicon via description: id=%s, description=%s', contig['id'], contig['description'])
            
            if('chromosome' in contig_description):
                contig['type'] = bc.REPLICON_CHROMOSOME
                log.debug('qc: detected chromosome replicon type via description: id=%s, description=%s', contig['id'], contig['description'])
            elif('plasmid' in contig_description):
                contig['type'] = bc.REPLICON_PLASMID
                log.debug('qc: detected plasmid replicon type via description: id=%s, description=%s', contig['id'], contig['description'])

            if(cfg.keep_contig_headers):
                if(contig['id'] in contig_ids):
                    log.error('Fasta import: duplicated contig id! contig-id=%s', contig['id'])
                    sys.exit(f"ERROR: Detected duplicated contig id! Contig ID ({contig['id']}) occures multiple times!")
                else:
                    contig_ids.add(contig['id'])
            else:
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
                    if(contig['type'] == bc.REPLICON_CHROMOSOME):
                        contig_desc.append('[location=chromosome]')
                contig_desc.append(f'[gcode={cfg.translation_table}]')
                contig['description'] = ' '.join(contig_desc)

            if(contig['complete'] and contig['topology'] == bc.TOPOLOGY_CIRCULAR):  # detection of chromosomes/plasmids via sequence length thresholds
                if(contig['length'] >= bc.REPLICON_LENGTH_THRESHOLD_CHROMOSOME):
                    contig['type'] = bc.REPLICON_CHROMOSOME
                    log.debug('qc: detected replicon type via length: id=%s, type=%s, length=%i, description=%s', contig['id'], contig['type'], contig['length'], contig['description'])
                elif(contig['length'] < bc.REPLICON_LENGTH_THRESHOLD_PLASMID):
                    contig['type'] = bc.REPLICON_PLASMID
                    log.debug('qc: detected replicon type via length: id=%s, type=%s, length=%i, description=%s', contig['id'], contig['type'], contig['length'], contig['description'])
            valid_contigs.append(contig)

            if(len(contigs) == 1 and contig['type'] == bc.REPLICON_PLASMID and cfg.plasmid is not None):
                contig['name'] = cfg.plasmid

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

            if(cfg.compliant):  # check INSDC compliance
                if(len(contig['id']) > 25):  # max 25 characters
                    log.error('INSDC compliance: contig id larger than 25! contig-id=%s', contig['id'])
                    sys.exit(f"ERROR: INSDC compliance failed! Contig ID ({contig['id']}) larger than 25 characers!")
                if(bc.RE_INSDC_ID.fullmatch(contig['id']) is None):  # invalid characters
                    log.error('INSDC compliance: contig id contains invalid characters! contig-id=%s', contig['id'])
                    sys.exit(f"ERROR: INSDC compliance failed! Contig ID ({contig['id']}) contains invalid characters!")

            log.info(
                "qc: revised sequence: id=%s, orig-id=%s, type=%s, complete=%s, topology=%s, name=%s, description='%s', orig-description='%s'",
                contig['id'], contig.get('orig_id', ''), contig['type'], contig['complete'], contig['topology'], contig.get('name', ''), contig['description'], contig.get('orig_description', '')
            )
    return valid_contigs, complete_genome


def extract_feature_sequence(feature: dict, contig: dict) -> str:
    if(feature.get('edge', False)):
        seq = contig['sequence'][feature['start']-1:] + contig['sequence'][:feature['stop']]
    else:
        seq = contig['sequence'][feature['start']-1:feature['stop']]
    if(feature['strand'] == bc.STRAND_REVERSE):
        seq = str(Seq(seq).reverse_complement())
    return seq
