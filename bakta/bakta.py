
import argparse
import logging
import multiprocessing as mp
import os
import sys
# import shutil
# import subprocess as sp
# from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import bakta
import bakta.constants as bc
import bakta.io as io
import bakta.utils as bu
import bakta.predictions as bp
import bakta.ups as ups
import bakta.psc as psc


def main():
    # parse arguments
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
    parser.add_argument('--pretty-json', action='store_true', dest='pretty_json', help='Write GFF3 annotation file')
    parser.add_argument('--gff3', action='store_true', help='Write GFF3 annotation file')
    parser.add_argument('--genbank', action='store_true', help='Write GenBank annotation file')
    parser.add_argument('--embl', action='store_true', help='Write EMBL annotation file')

    parser.add_argument('--locus', '-l', action='store', default='', help='Locus tag prefix')
    parser.add_argument('--genus', action='store', default='', help='Genus name')
    parser.add_argument('--species', action='store', default='', help='Species name')
    parser.add_argument('--strain', action='store', default='', help='Strain name')
    parser.add_argument('--plasmid', action='store', default='', help='Plasmid name')
    parser.add_argument('--gram', action='store', default='', choices=['+', '-'], help="Gram type: ''/+/- (default = '')")
    parser.add_argument('--complete', action='store_true', help="Replicons (chromosome/plasmid[s]) are complete")

    parser.add_argument('--version', action='version', version='%(prog)s ' + bakta.__version__)
    parser.add_argument('--citation', action='store_true', help='Print citation')
    args = parser.parse_args()

    # print citation
    if(args.citation):
        print(bc.CITATION)
        sys.exit()

    ############################################################################
    # Setup logging
    ############################################################################
    logging.basicConfig(
        filename='bakta.log',
        filemode='w',
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        level=logging.DEBUG if args.verbose else logging.INFO
    )
    log = logging.getLogger('main')
    log.info('version %s', bakta.__version__)

    ############################################################################
    # Checks and configurations
    # - setup global configuration
    # - test database
    # - test binary dependencies
    # - check parameters
    ############################################################################
    config = bu.setup_configuration()
    config['threads'] = args.threads

    if(args.db):
        db_path = Path(args.db).resolve()
        config['db'] = db_path
    bu.test_database(config)

    if('bundled-binaries' not in config):
        bu.test_dependencies()

    genome_path = Path(args.genome).resolve()
    if(not os.access(str(genome_path), os.R_OK)):
        log.error('genome file not readable! path=%s', genome_path)
        sys.exit('ERROR: genome file (%s) not readable!' % genome_path)
    if(genome_path.stat().st_size == 0):
        log.error('empty genome file! path=%s', genome_path)
        sys.exit('ERROR: genome file (%s) is empty!' % genome_path)

    output_path = Path(args.output) if args.output else Path.cwd()
    if(not output_path.exists()):
        output_path.mkdir(parents=True, exist_ok=True)
        log.info('create output dir: path=%s', output_path)
    output_path = output_path.resolve()

    log.info('configuration: db-path=%s', config['db'])
    log.info('configuration: bundled binaries=%s', config['bundled-binaries'])
    log.info('configuration: tmp-path=%s', config['tmp'])
    log.info('parameters: genome=%s', genome_path)
    log.info('parameters: output=%s', output_path)
    log.info('parameters: prefix=%s', args.prefix)
    log.info('parameters: locus=%s', args.locus)
    log.info('parameters: genus=%s', args.genus)
    log.info('parameters: species=%s', args.species)
    log.info('parameters: strain=%s', args.strain)
    log.info('parameters: complete=%s', args.complete)
    log.info('options: gff3=%s', args.gff3)
    log.info('options: genbank=%s', args.genbank)
    log.info('options: embl=%s', args.embl)
    log.info('options: threads=%d', args.threads)
    log.info('options: pretty-json=%d', args.pretty_json)
    if(args.verbose):
        print("Bakta v%s" % bakta.__version__)
        print('Options and arguments:')
        print("\tuse bundled binaries: %s" % str(config['bundled-binaries']))
        print("\tdb path: %s" % str(config['db']))
        print("\tgenome path: %s" % str(genome_path))
        print("\toutput path: %s" % str(output_path))
        print("\ttmp path: %s" % str(config['tmp']))
        print("\t# threads: %d" % args.threads)
        print("\tcomplete replicons: %s" % args.complete)

    ############################################################################
    # Parse input genome
    # - parse contigs in Fasta file
    # - apply contig length filter
    # - rename contigs
    ############################################################################
    print('parse genome...')
    try:
        contigs, discarded_contigs = io.import_contigs(genome_path, args.min_contig_length)
    except:
        log.error('wrong genome file format!', exc_info=True)
        sys.exit('ERROR: wrong genome file format!')
    log.info('imported contigs: # valid=%d, # discarded=%d', len(contigs), len(discarded_contigs))
    print("\timported %i valid contig(s)\n\tdiscarded %i unvalid contig(s)" % (len(contigs), len(discarded_contigs)))
    if(len(contigs) == 0):
        log.warning('no valid contigs!')
        sys.exit('Error: input file contains no valid contigs.')
    contigs_path = config['tmp'].joinpath('contigs.fna')
    io.export_contigs(contigs, contigs_path)
    data = {
        'genome_size': sum(map(lambda k: k['length'], contigs)),
        'contigs': contigs
    }

    ############################################################################
    # tRNA prediction
    ############################################################################
    # print('predict tRNAs...')
    # log.debug('start tRNA prediction')
    # data['t_rnas'] = bp.predict_t_rnas(config, data, contigs_path)
    # print("\tfound %i tRNAs" % len(data['t_rnas']))

    ############################################################################
    # rRNA prediction
    ############################################################################
    # print('predict rRNAs...')
    # log.debug('start rRNA prediction')
    # data['r_rnas'] = bp.predict_r_rnas(config, data, contigs_path)
    # print("\tfound %i rRNAs" % len(data['r_rnas']))

    ############################################################################
    # ncRNA prediction
    ############################################################################
    # print('predict ncRNAs...')
    # log.debug('start ncRNA prediction')
    # data['nc_rnas'] = bp.predict_nc_rnas(config, data, contigs_path)
    # print("\tfound %i ncRNAs" % len(data['nc_rnas']))

    ############################################################################
    # CRISPR prediction
    ############################################################################
    # print('predict CRISPR cassettes...')
    # log.debug('start CRISPR prediction')
    # data['crisprs'] =  bp.predict_crispr(config, data, contigs_path)
    # print("\tfound %i CRISPR cassettes" % len(data['crisprs']))

    ############################################################################
    # CDS prediction
    # - Prodigal prediction
    # - lookup UPS matches for CDSs
    # - search PSC for unannotated CDS
    ############################################################################
    print('predict CDSs...')
    log.debug('start CDS prediction')
    data['cdss'] = bp.predict_cdss(config, data['contigs'], contigs_path)
    print("\tfound %i CDSs" % len(data['cdss']))
    cdss_found, cdss_not_found = ups.lookup_upss(config, data['cdss'])
    print("\tfound %i UPSs for CDSs" % len(cdss_found))
    cdss_found, cdss_not_found = psc.search_pscs(config, cdss_not_found)
    print("\tfound %i PSCs for CDSs" % len(cdss_found))
    psc.lookup_pscs(config, data['cdss'])

    ############################################################################
    # ORF prediction
    # - in-mem ORF extraction
    # - overlap filtering (tRNA, rRNA, CDS)
    # - lookup UPS matches for ORFs
    # - filter ORFs w/o UPS match
    ############################################################################
    print('predict ORFs...')
    log.debug('start ORF prediction')
    orfs = bp.extract_orfs(data['contigs'])
    print("\tfound %i ORFs" % len(orfs))

    print('filter ORFs...')
    log.debug('start ORF filtering')
    orfs, discarded_orfs = bp.overlap_filter_orfs(data, orfs)
    print("\tdiscarded %i ORFs, %i remaining" % (len(discarded_orfs), len(orfs)))

    orfs_found, orfs_not_found = ups.lookup_upss(config, orfs)
    print("\tfound %i UPSs for ORFs, %i discarded" % (len(orfs_found), len(orfs_not_found)))
    for orf in orfs_found:
        orf['type'] = 'cds'  # change seq type from ORF to CDS
        orf['inference'] = 'UniProtKB'
    data['cdss'] += orfs_found  # add ORFs identified by an UPS to CDSs

    ############################################################################
    # Create annotations
    # - filter features based on precedence and overlaps
    # - sort features
    # - create locus tags for features
    ############################################################################
    print('select features and create locus tags...')
    log.debug('start feature selection and creation of locus tags')
    annotations = {}
    locus_tag_nr = 1
    locus_prefix = args.locus if args.locus != '' else bu.create_locus_prefix(contigs)

    # tRNAs
    for t_rna in data['t_rnas']:
        locus_tag = "%s%04d" % (locus_prefix, locus_tag_nr)
        locus_tag_nr += 1
        t_rna['locus'] = locus_tag
        annotations[locus_tag] = t_rna

    # rRNAs
    for r_rna in data['r_rnas']:
        locus_tag = "%s%04d" % (locus_prefix, locus_tag_nr)
        locus_tag_nr += 1
        r_rna['locus'] = locus_tag
        annotations[locus_tag] = r_rna

    # ncRNAs
    for nc_rna in data['nc_rnas']:
        locus_tag = "%s%04d" % (locus_prefix, locus_tag_nr)
        locus_tag_nr += 1
        nc_rna['locus'] = locus_tag
        annotations[locus_tag] = nc_rna

    ############################################################################
    # Write output files
    # - write comprehensive annotation results as JSON
    # - write optional output files in GFF3/GenBank/EMBL formats
    # - remove temp directory
    ############################################################################
    print('write JSON output...')
    log.debug('write JSON output')

    prefix = genome_path.stem if args.prefix == '' else args.prefix
    json_path = output_path.joinpath("%s.json" % prefix)
    io.write_json(annotations, json_path, args.pretty_json)

    if(args.gff3):
        print('write GFF3 output...')
        log.debug('write GFF3 output')
        gff3_path = output_path.resolve("%s.gff3" % prefix)
        io.write_gff3(annotations, gff3_path)

    if(args.genbank):
        print('write GenBank output...')
        log.debug('write GenBank output')
        genbank_path = output_path.resolve("%s.gbff" % prefix)
        io.write_genbank(annotations, genbank_path)

    if(args.embl):
        print('write EMBL output...')
        log.debug('write EMBL output')
        embl_path = output_path.resolve("%s.embl" % prefix)
        io.write_embl(annotations, embl_path)

    # remove tmp dir
    # shutil.rmtree(str(config['tmp']))
    log.debug('removed tmp dir: %s', config['tmp'])


if __name__ == '__main__':
    main()
