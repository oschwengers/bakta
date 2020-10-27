import logging
import sys
import shutil
from pathlib import Path

import bakta
import bakta.constants as bc
import bakta.config as cfg
import bakta.io.fasta as fasta
import bakta.io.json as json
import bakta.io.tsv as tsv
import bakta.io.gff as gff
import bakta.io.embl as embl
import bakta.io.genbank as genbank
import bakta.features.annotation as anno
import bakta.features.t_rna as t_rna
import bakta.features.tm_rna as tm_rna
import bakta.features.r_rna as r_rna
import bakta.features.nc_rna as nc_rna
import bakta.features.nc_rna_region as nc_rna_region
import bakta.features.crispr as crispr
import bakta.features.orf as orf
import bakta.features.cds as cds
import bakta.features.s_orf as s_orf
import bakta.features.gaps as gaps
import bakta.features.ori as ori
import bakta.utils as bu
import bakta.ups as ups
import bakta.ips as ips
import bakta.psc as psc


def main(args):

    ############################################################################
    # Setup logging
    ############################################################################
    prefix = args.prefix if args.prefix else Path(args.genome).stem
    logging.basicConfig(
        filename='%s.log' % prefix,
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
            ('complete replicons', cfg.complete),
            ('prodigal training file', cfg.prodigal_tf)
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
        contigs, discarded_contigs = fasta.import_contigs(cfg.genome_path, cfg.min_contig_length)
    except:
        log.error('wrong genome file format!', exc_info=True)
        sys.exit('ERROR: wrong genome file format!')
    log.info('imported contigs: # valid=%d, # discarded=%d', len(contigs), len(discarded_contigs))
    print("\timported %i valid contig(s)\n\tdiscarded %i unvalid contig(s)" % (len(contigs), len(discarded_contigs)))
    if(len(contigs) == 0):
        log.warning('no valid contigs!')
        sys.exit('Error: input file contains no valid contigs.')
    contigs_path = cfg.tmp_path.joinpath('contigs.fna')
    fasta.export_contigs(contigs, contigs_path)
    data = {
        'genome_size': sum(map(lambda k: k['length'], contigs)),
        'contigs': contigs
    }

    ############################################################################
    # tRNA prediction
    ############################################################################
    if(cfg.skip_trna):
        print('skip tRNA prediction...')
    else:
        print('predict tRNAs...')
        log.debug('start tRNA prediction')
        data[bc.FEATURE_T_RNA] = t_rna.predict_t_rnas(data, contigs_path)
        print("\tfound %i tRNAs" % len(data[bc.FEATURE_T_RNA]))

    ############################################################################
    # tmRNA prediction
    ############################################################################
    if(cfg.skip_tmrna):
        print('skip tmRNA prediction...')
    else:
        print('predict tmRNAs...')
        log.debug('start tmRNA prediction')
        data[bc.FEATURE_TM_RNA] = tm_rna.predict_tm_rnas(data, contigs_path)
        print("\tfound %i tmRNAs" % len(data[bc.FEATURE_TM_RNA]))

    ############################################################################
    # rRNA prediction
    ############################################################################
    if(cfg.skip_rrna):
        print('skip rRNA prediction...')
    else:
        print('predict rRNAs...')
        log.debug('start rRNA prediction')
        data[bc.FEATURE_R_RNA] = r_rna.predict_r_rnas(data, contigs_path)
        print("\tfound %i rRNAs" % len(data[bc.FEATURE_R_RNA]))

    ############################################################################
    # ncRNA gene prediction
    ############################################################################
    if(cfg.skip_ncrna):
        print('skip ncRNA prediction...')
    else:
        print('predict ncRNAs...')
        log.debug('start ncRNA prediction')
        data[bc.FEATURE_NC_RNA] = nc_rna.predict_nc_rnas(data, contigs_path)
        print("\tfound %i ncRNAs" % len(data[bc.FEATURE_NC_RNA]))

    ############################################################################
    # ncRNA region prediction
    ############################################################################
    if(cfg.skip_ncrna_region):
        print('skip ncRNA region prediction...')
    else:
        print('predict ncRNA regions...')
        log.debug('start ncRNA region prediction')
        data[bc.FEATURE_NC_RNA_REGION] = nc_rna_region.predict_nc_rna_regions(data, contigs_path)
        print("\tfound %i ncRNA regions" % len(data[bc.FEATURE_NC_RNA_REGION]))

    ############################################################################
    # CRISPR prediction
    ############################################################################
    if(cfg.skip_crispr):
        print('skip CRISPR array prediction...')
    else:
        print('predict CRISPR arrays...')
        log.debug('start CRISPR prediction')
        data[bc.FEATURE_CRISPR] = crispr.predict_crispr(data, contigs_path)
        print("\tfound CRISPR arrays: %i" % len(data[bc.FEATURE_CRISPR]))

    ############################################################################
    # CDS prediction
    # - Prodigal prediction
    # - lookup UPS matches
    # - lookup IPS matches
    # - search PSC for unannotated CDSs
    ############################################################################
    if(cfg.skip_cds):
        print('skip CDS prediction...')
    else:
        print('predict CDSs...')
        log.debug('start CDS prediction')
        data[bc.FEATURE_CDS] = cds.predict(data['contigs'], contigs_path)
        print("\tpredicted CDSs: %i " % len(data[bc.FEATURE_CDS]))
        discarded_cds = orf.detect_spurious(data[bc.FEATURE_CDS])
        print("\tdiscarded spurious CDSs: %i " % len(discarded_cds))
        cdss_ups, cdss_not_found = ups.lookup(data[bc.FEATURE_CDS])
        print("\tfound UPSs: %i " % len(cdss_ups))
        cdss_ips, tmp = ips.lookup(cdss_ups)
        cdss_not_found.extend(tmp)
        print("\tfound IPSs: %i " % len(cdss_ips))
        if(len(cdss_not_found) > 0):
            cdss_psc, cdss_not_found = psc.search(cdss_not_found)
            print("\tfound PSCs: %i " % len(cdss_psc))
        print("lookup PSC annotations for PSCs and IPSs and mark hypotheticals...")
        psc.lookup(data[bc.FEATURE_CDS])  # lookup PSC info
        cds.mark_hypotheticals(data[bc.FEATURE_CDS])  # mark hypotheticals
        for feat in data[bc.FEATURE_CDS]:
            anno.combine_ips_psc_annotation(feat) # combine IPS and PSC annotations
    
    ############################################################################
    # sORF prediction
    # - in-mem sORF extraction
    # - overlap filtering (tRNA, tmRNA, rRNA, CDS)
    # - lookup UPS matches
    # - lookup IPS matches
    # - filter sORFs w/o IPS match
    ############################################################################
    if(cfg.skip_sorf):
        print('skip sORF prediction...')
    else:
        print('extract sORF...')
        log.debug('start sORF prediction')
        sorfs = s_orf.extract(data['contigs'])
        print("\textracted potential sORFs: %i" % len(sorfs))

        print('apply sORF overlap filter...')
        log.debug('apply sORF overlap filter')
        sorfs, discarded_sorfs = s_orf.overlap_filter(data, sorfs)
        print("\tdiscarded overlapping sORFs: %i, %i remaining" % (len(discarded_sorfs), len(sorfs)))

        discarded_sorfs = orf.detect_spurious(sorfs)
        print("\tdiscarded spurious sORFs: %i " % len(discarded_sorfs))

        print('lookup sORFs UPSs/IPSs...')
        sorf_upss, sorfs_not_found = ups.lookup(sorfs)
        sorf_ipss, tmp = ips.lookup(sorf_upss)
        sorfs_not_found.extend(tmp)
        print("\tfound UPSs/IPSs: %i" % len(sorf_ipss))

        sorf_pscs = []
        if(len(sorfs_not_found) > 0):
            print('search sORFs PSCs...')
            tmp, sorfs_not_found = s_orf.search_pscs(sorfs_not_found)
            sorf_pscs.extend(tmp)
            print("\tfound PSCs: %i" % len(sorf_pscs))

        print("lookup PSC annotations for sORFs...")
        sorf_pscs.extend(sorf_ipss)
        psc.lookup(sorf_pscs)  # lookup PSC info
        sorfs_filtered = s_orf.annotation_filter(sorfs)
        data[bc.FEATURE_SORF] = sorfs_filtered
        print("\tfiltered sORFs: %i" % len(sorfs_filtered))
        for feat in data[bc.FEATURE_SORF]:
            anno.combine_ips_psc_annotation(feat) # combine IPS and PSC annotations
    
    ############################################################################
    # gap annotation
    # - in-mem gap detection
    # - gap annotation
    ############################################################################
    if(cfg.skip_gap):
        print('skip gap annotation...')
    else:
        print('detect gaps...')
        log.debug('start gap detection')
        assembly_gaps = gaps.detect_assembly_gaps(data['contigs'])
        print("\tfound gaps: %i" % len(assembly_gaps))
        data[bc.FEATURE_GAP] = assembly_gaps
    
    ############################################################################
    # oriC/T prediction
    ############################################################################
    if(cfg.skip_ori):
        print('skip oriC/T annotation...')
    else:
        print('detect oris...')
        log.debug('start oriC detection')
        oriCs = ori.predict_oris(data, contigs_path, bc.FEATURE_ORIC)
        print("\tfound oriCs: %i" % len(oriCs))
        data[bc.FEATURE_ORIC] = oriCs

        log.debug('start oriT detection')
        oriTs = ori.predict_oris(data, contigs_path, bc.FEATURE_ORIT)
        print("\tfound oriTs: %i" % len(oriTs))
        data[bc.FEATURE_ORIT] = oriTs

    ############################################################################
    # Filter overlapping features
    ############################################################################
    print('apply feature overlap filters...')
    anno.detect_feature_overlaps(data)

    ############################################################################
    # Create annotations
    # - filter features based on precedence and overlaps
    # - sort features
    # - create locus tags for features
    ############################################################################
    print('select features and create locus tags...')
    log.debug('start feature selection and creation of locus tags')
    features_by_contig = {k['id']: [] for k in contigs}
    for feature_type in [
            bc.FEATURE_T_RNA,
            bc.FEATURE_TM_RNA,
            bc.FEATURE_R_RNA,
            bc.FEATURE_NC_RNA,
            bc.FEATURE_NC_RNA_REGION,
            bc.FEATURE_CRISPR,
            bc.FEATURE_CDS,
            bc.FEATURE_SORF,
            bc.FEATURE_GAP,
            bc.FEATURE_ORIC,
            bc.FEATURE_ORIT
        ]:
        feature_list = data.get(feature_type, [])
        for feature in feature_list:
            if('discarded' not in feature):
                contig_features = features_by_contig.get(feature['contig'])
                contig_features.append(feature)
    features = []
    for contig in contigs:
        contig_features = features_by_contig[contig['id']]
        contig_features.sort(key=lambda k: k['start'])
        features.extend(contig_features)

    locus_tag_nr = 5
    # use user provided locus tag if not None/non-empty or generate a sequence based locus prefix
    locus_tag_prefix = cfg.locus_tag if cfg.locus_tag else bu.create_locus_tag_prefix(contigs)
    log.info('locus tag prefix: %s', locus_tag_prefix)
    for feature in features:
        locus_tag = "%s_%05d" % (locus_tag_prefix, locus_tag_nr)
        if(feature['type'] != bc.FEATURE_GAP and feature['type'] != bc.FEATURE_CRISPR and feature['type'] != bc.FEATURE_ORIC and feature['type'] != bc.FEATURE_ORIT):
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

    log.info('file prefix: %s', prefix)
    json_path = cfg.output_path.joinpath("%s.json" % prefix)
    json.write_json(features, json_path)

    if(cfg.tsv):
        print('write TSV output...')
        log.debug('write tsv output')
        tsv_path = cfg.output_path.joinpath("%s.tsv" % prefix)
        tsv.write_tsv(contigs, features_by_contig, tsv_path)

    if(cfg.gff3):
        print('write GFF3 output...')
        log.debug('write GFF3 output')
        gff3_path = cfg.output_path.joinpath("%s.gff3" % prefix)
        gff.write_gff3(contigs, features_by_contig, gff3_path)

    if(cfg.genbank):
        print('write GenBank output...')
        log.debug('write GenBank output')
        genbank_path = cfg.output_path.joinpath("%s.gbff" % prefix)
        genbank.write_genbank(features, genbank_path)

    if(cfg.embl):
        print('write EMBL output...')
        log.debug('write EMBL output')
        embl_path = cfg.output_path.joinpath("%s.embl" % prefix)
        embl.write_embl(features, embl_path)

    if(cfg.faa):
        print('write translated CDS sequences...')
        log.debug('write translated CDS output')
        faa_path = cfg.output_path.joinpath("%s.faa" % prefix)
        fasta.write_faa(features, faa_path)

    ############################################################################
    # Print statistics
    # - genome stats
    # - annotation stats
    ############################################################################
    print('\ngenome statistics:')
    genome_stats = bu.calc_genome_stats(data, features)
    print('\tGenome size: %i bp' % genome_stats['genome_size'])
    print('\tContigs/replicons: %i' % genome_stats['no_contigs'])
    print('\tGC: %2.1f%%' % (100 * genome_stats['gc_ratio']))
    print('\tN50: %i' % genome_stats['n50'])
    print('\tN ratio: %2.1f%%' % (100 * genome_stats['n_ratio']))
    print('\tcoding density: %2.1f%%' % (100 * genome_stats['coding_ratio']))

    print('\nannotation statistics:')
    print('\ttRNAs: %i' % len([f for f in features if f['type'] == bc.FEATURE_T_RNA]))
    print('\ttmRNAs: %i' % len([f for f in features if f['type'] == bc.FEATURE_TM_RNA]))
    print('\trRNAs: %i' % len([f for f in features if f['type'] == bc.FEATURE_R_RNA]))
    print('\tncRNAs: %i' % len([f for f in features if f['type'] == bc.FEATURE_NC_RNA]))
    print('\tncRNA-regions: %i' % len([f for f in features if f['type'] == bc.FEATURE_NC_RNA_REGION]))
    print('\tCRISPR arrays: %i' % len([f for f in features if f['type'] == bc.FEATURE_CRISPR]))
    print('\tCDSs: %i' % len([f for f in features if f['type'] == bc.FEATURE_CDS]))
    print('\tsORFs: %i' % len([f for f in features if f['type'] == bc.FEATURE_SORF]))
    print('\tgaps: %i' % len([f for f in features if f['type'] == bc.FEATURE_GAP]))
    print('\toriCs: %i' % len([f for f in features if f['type'] == bc.FEATURE_ORIC]))
    print('\toriTs: %i' % len([f for f in features if f['type'] == bc.FEATURE_ORIT]))

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
