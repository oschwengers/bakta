
import logging
import sys
import shutil

import bakta
import bakta.constants as bc
import bakta.config as cfg
import bakta.io as io
import bakta.utils as bu
import bakta.predictions as bp
import bakta.ups as ups
import bakta.psc as psc


def main(args):

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
    log.info('command line: %s', ' '.join(sys.argv))

    ############################################################################
    # Checks and configurations
    # - check parameters and setup global configuration
    # - test database
    # - test binary dependencies
    ############################################################################
    cfg.setup(args)  # check parameters and prepare global configuration
    bu.test_database()
    bu.test_dependencies()
    if(cfg.verbose):
        print("Bakta v%s" % bakta.__version__)
        print('Options and arguments:')
        for label, value in [
            ('db path', cfg.db_path),
            ('genome path', cfg.genome_path),
            ('output path', cfg.output_path),
            ('tmp path', cfg.tmp_path),
            ('# threads', cfg.threads),
            ('complete replicons', cfg.complete)
        ]:
            print("\t%s: %s" % (label, str(value)))

    ############################################################################
    # Parse input genome
    # - parse contigs in Fasta file
    # - apply contig length filter
    # - rename contigs
    ############################################################################
    print('parse genome...')
    try:
        contigs, discarded_contigs = io.import_contigs(cfg.genome_path, cfg.min_contig_length)
    except:
        log.error('wrong genome file format!', exc_info=True)
        sys.exit('ERROR: wrong genome file format!')
    log.info('imported contigs: # valid=%d, # discarded=%d', len(contigs), len(discarded_contigs))
    print("\timported %i valid contig(s)\n\tdiscarded %i unvalid contig(s)" % (len(contigs), len(discarded_contigs)))
    if(len(contigs) == 0):
        log.warning('no valid contigs!')
        sys.exit('Error: input file contains no valid contigs.')
    contigs_path = cfg.tmp_path.joinpath('contigs.fna')
    io.export_contigs(contigs, contigs_path)
    data = {
        'genome_size': sum(map(lambda k: k['length'], contigs)),
        'contigs': contigs
    }

    ############################################################################
    # tRNA prediction
    ############################################################################
    print('predict tRNAs...')
    log.debug('start tRNA prediction')
    data[bc.FEATURE_T_RNA] = bp.predict_t_rnas(data, contigs_path)
    print("\tfound %i tRNAs" % len(data[bc.FEATURE_T_RNA]))

    ############################################################################
    # tmRNA prediction
    ############################################################################
    print('predict tmRNAs...')
    log.debug('start tmRNA prediction')
    data[bc.FEATURE_TM_RNA] = bp.predict_tm_rnas(data, contigs_path)
    print("\tfound %i tmRNAs" % len(data[bc.FEATURE_TM_RNA]))

    ############################################################################
    # rRNA prediction
    ############################################################################
    print('predict rRNAs...')
    log.debug('start rRNA prediction')
    data[bc.FEATURE_R_RNA] = bp.predict_r_rnas(data, contigs_path)
    print("\tfound %i rRNAs" % len(data[bc.FEATURE_R_RNA]))

    ############################################################################
    # ncRNA gene prediction
    ############################################################################
    print('predict ncRNA genes...')
    log.debug('start ncRNA gene prediction')
    data[bc.FEATURE_NC_RNA_GENE] = bp.predict_nc_rna_genes(data, contigs_path)
    print("\tfound %i ncRNA genes" % len(data[bc.FEATURE_NC_RNA_GENE]))

    ############################################################################
    # ncRNA region prediction
    ############################################################################
    print('predict ncRNA regions...')
    log.debug('start ncRNA region prediction')
    data[bc.FEATURE_NC_RNA_REGION] = bp.predict_nc_rna_regions(data, contigs_path)
    print("\tfound %i ncRNA regions" % len(data[bc.FEATURE_NC_RNA_REGION]))

    ############################################################################
    # CRISPR prediction
    ############################################################################
    # print('predict CRISPR cassettes...')
    # log.debug('start CRISPR prediction')
    # data[bc.FEATURE_CRISPR] = bp.predict_crispr(data, contigs_path)
    # print("\tfound %i CRISPR cassettes" % len(data[bc.FEATURE_CRISPR]))

    ############################################################################
    # CDS prediction
    # - Prodigal prediction
    # - lookup UPS matches for CDSs
    # - search PSC for unannotated CDS
    ############################################################################
    print('predict CDSs...')
    log.debug('start CDS prediction')
    data[bc.FEATURE_CDS] = bp.predict_cdss(data['contigs'], contigs_path)
    print("\tfound %i CDSs" % len(data[bc.FEATURE_CDS]))
    upss_found, cdss_not_found = ups.lookup_upss(data[bc.FEATURE_CDS])
    print("\tfound %i UPSs for CDSs" % len(upss_found))
    pscs_found, cdss_not_found = psc.search_pscs(cdss_not_found)
    print("\tfound %i PSCs for CDSs" % len(pscs_found))
    
    ############################################################################
    # sORF prediction
    # - in-mem sORF extraction
    # - overlap filtering (tRNA, rRNA, CDS)
    # - lookup UPS matches for sORFs
    # - filter sORFs w/o UPS match
    ############################################################################
    print('predict short open reading frames (sORFs)...')
    log.debug('start sORF prediction')
    orfs = bp.extract_sorfs(data['contigs'])
    print("\tfound %i potential sORFs" % len(orfs))

    print('filter potential sORFs by overlaps...')
    log.debug('start sORF filtering')
    orfs, discarded_orfs = bp.overlap_filter_sorfs(data, orfs)
    print("\tdiscarded %i sORFs, %i remaining" % (len(discarded_orfs), len(orfs)))

    orfs_found, orfs_not_found = ups.lookup_upss(orfs)
    print("\tfound %i UPSs for sORFs, %i discarded" % (len(orfs_found), len(orfs_not_found)))
    for orf in orfs_found:
        orf['type'] = bc.FEATURE_CDS  # change feature type from ORF to CDS
        orf['inference'] = 'UniRef100'
    data[bc.FEATURE_CDS] += orfs_found  # add ORFs identified by an UPS to CDSs

    print("lookup PSC annotations for PSCs, UPSs and sORFs...")
    psc.lookup_pscs(data[bc.FEATURE_CDS])  # lookup PSC info

    ############################################################################
    # Create annotations
    # - filter features based on precedence and overlaps
    # - sort features
    # - create locus tags for features
    ############################################################################
    print('select features and create locus tags...')
    log.debug('start feature selection and creation of locus tags')
    features_by_contig = {}
    for feature_list in [
            data['t_rnas'],
            data['r_rnas'],
            data['nc_rna_genes'],
            # data['crisprs'],
            data['cdss']
        ]:
        for feature in feature_list:
            contig_features = features_by_contig.get(feature['contig'])
            if(contig_features is None):
                contig_features = []
                features_by_contig[feature['contig']] = contig_features
            contig_features.append(feature)
    features = []
    for contig in contigs:
        contig_features = features_by_contig[contig['id']]
        features.extend(sorted(contig_features, key=lambda k: k['start']))

    locus_tag_nr = 5
    locus_prefix = bu.create_locus_tag_prefix(contigs)
    for feature in features:
        locus_tag = "%s%04i" % (locus_prefix, locus_tag_nr)
        feature['locus'] = locus_tag
        locus_tag_nr += 5

    ############################################################################
    # Write output files
    # - write comprehensive annotation results as JSON
    # - write optional output files in GFF3/GenBank/EMBL formats
    # - remove temp directory
    ############################################################################
    print('write JSON output...')
    log.debug('write JSON output')

    prefix = cfg.genome_path.stem if cfg.prefix is None else cfg.prefix
    json_path = cfg.output_path.joinpath("%s.json" % prefix)
    io.write_json(features, json_path)

    if(cfg.gff3):
        print('write GFF3 output...')
        log.debug('write GFF3 output')
        gff3_path = cfg.output_path.resolve("%s.gff3" % prefix)
        io.write_gff3(features, gff3_path)

    if(cfg.genbank):
        print('write GenBank output...')
        log.debug('write GenBank output')
        genbank_path = cfg.output_path.resolve("%s.gbff" % prefix)
        io.write_genbank(features, genbank_path)

    if(cfg.embl):
        print('write EMBL output...')
        log.debug('write EMBL output')
        embl_path = cfg.output_path.resolve("%s.embl" % prefix)
        io.write_embl(features, embl_path)

    # remove tmp dir
    shutil.rmtree(str(cfg.tmp_path))
    log.debug('removed tmp dir: %s', cfg.tmp_path)


if __name__ == '__main__':
    # parse arguments
    args = bu.parse_arguments()

    if(args.citation):  # print citation
        print(bc.CITATION)
        sys.exit()
    else:  # start bakta
        main(args)
