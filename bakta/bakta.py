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


def main():

    # parse arguments
    args = bu.parse_arguments()
    
    if(args.citation):  # print citation
        print(bc.CITATION)
        sys.exit()

    ############################################################################
    # Setup logging
    ############################################################################
    prefix = args.prefix if args.prefix else Path(args.genome).stem
    logging.basicConfig(
        filename=f'{prefix}.log',
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
    # bu.test_database()
    bu.test_dependencies()
    if(cfg.verbose):
        print(f'Bakta v{bakta.__version__}')
        print('Options and arguments:')
        print(f'\tinput: {cfg.genome_path}')
        print(f'\tdb: {cfg.db_path}')
        print(f'\toutput: {cfg.output_path}')
        print(f'\ttmp directory: {cfg.tmp_path}')
        print(f'\t# threads: {cfg.threads}')
        if(cfg.keep_contig_headers): print(f'\tkeep contig headers: {cfg.keep_contig_headers}')
        if(cfg.locus): print(f'\tlocus prefix: {cfg.locus}')
        if(cfg.locus_tag): print(f'\tlocus tag prefix: {cfg.locus_tag}')
        if(cfg.complete): print(f'\tcomplete replicons: {cfg.complete}')
        if(cfg.prodigal_tf): print(f'\tprodigal training file: {cfg.prodigal_tf}')
        if(cfg.replicons): print(f'\treplicon table: {cfg.replicons}')
        print(f'\ttranslation table: {cfg.translation_table}')
        taxon = f"{cfg.genus if cfg.genus else ''} {cfg.species if cfg.species else ''} {cfg.strain if cfg.strain else ''}"
        if(cfg.plasmid): taxon = f'{taxon} plasmid {cfg.plasmid}'
        taxon = ' '.join(taxon.split())
        if(taxon): print(f'\ttaxon: {taxon}')
        if(cfg.gram != '?'): print(f'\tgram: {cfg.gram}')
        if(cfg.skip_trna): print(f'\tskip tRNA: {cfg.skip_trna}')
        if(cfg.skip_tmrna): print(f'\tskip tmRNA: {cfg.skip_tmrna}')
        if(cfg.skip_rrna): print(f'\tskip rRNA: {cfg.skip_rrna}')
        if(cfg.skip_ncrna): print(f'\tskip ncRNA: {cfg.skip_ncrna}')
        if(cfg.skip_ncrna_region): print(f'\tskip ncRNA region: {cfg.skip_ncrna_region}')
        if(cfg.skip_crispr): print(f'\tskip CRISPR: {cfg.skip_crispr}')
        if(cfg.skip_cds): print(f'\tskip CDS: {cfg.skip_cds}')
        if(cfg.skip_sorf): print(f'\tskip sORF: {cfg.skip_sorf}')
        if(cfg.skip_gap): print(f'\tskip gap: {cfg.skip_gap}')
        if(cfg.skip_ori): print(f'\tskip oriC/V/T: {cfg.skip_ori}')

    ############################################################################
    # Import genome
    # - parse contigs in Fasta file
    # - apply contig length filter
    # - rename contigs
    ############################################################################
    print('parse genome sequences...')
    try:
        contigs = fasta.import_contigs(cfg.genome_path)
        log.info('imported sequences=%i', len(contigs))
        print(f'\timported: {len(contigs)}')
    except:
        log.error('wrong genome file format!', exc_info=True)
        sys.exit('ERROR: wrong genome file format!')
    replicons = bu.parse_replicon_table(cfg.replicons) if cfg.replicons else None
    contigs = bu.qc_contigs(contigs, replicons)
    print(f'\tfiltered & revised: {len(contigs)}')
    no_chromosomes = len([c for c in contigs if c['type'] == bc.REPLICON_CHROMOSOME])
    if(no_chromosomes > 0):
        print(f"\tchromosomes: {no_chromosomes}")
    no_plasmids = len([c for c in contigs if c['type'] == bc.REPLICON_PLASMID])
    if(no_plasmids > 0):
        print(f"\tplasmids: {no_plasmids}")
    no_contigs = len([c for c in contigs if c['type'] == bc.REPLICON_CONTIG])
    if(no_contigs > 0):
        print(f"\tcontigs: {no_contigs}")
    if(len(contigs) == 0):
        log.warning('no valid contigs!')
        sys.exit('Error: input file contains no valid contigs.')
    contigs_path = cfg.tmp_path.joinpath('contigs.fna')
    fasta.export_contigs(contigs, contigs_path)
    genome = {
        'genus': cfg.genus,
        'species': cfg.species,
        'strain': cfg.strain,
        'gram': cfg.gram,
        'translation_table': cfg.translation_table,
        'size': sum([c['length'] for c in contigs]),
        'features': {},
        'contigs': contigs
    }
    if(cfg.plasmid):
        genome['plasmid'] = cfg.plasmid
    
    ############################################################################
    # tRNA prediction
    ############################################################################
    if(cfg.skip_trna):
        print('skip tRNA prediction...')
    else:
        print('predict tRNAs...')
        log.debug('start tRNA prediction')
        genome['features'][bc.FEATURE_T_RNA] = t_rna.predict_t_rnas(genome, contigs_path)
        print(f"\tfound: {len(genome['features'][bc.FEATURE_T_RNA])}")

    ############################################################################
    # tmRNA prediction
    ############################################################################
    if(cfg.skip_tmrna):
        print('skip tmRNA prediction...')
    else:
        print('predict tmRNAs...')
        log.debug('start tmRNA prediction')
        genome['features'][bc.FEATURE_TM_RNA] = tm_rna.predict_tm_rnas(genome, contigs_path)
        print(f"\tfound: {len(genome['features'][bc.FEATURE_TM_RNA])}")

    ############################################################################
    # rRNA prediction
    ############################################################################
    if(cfg.skip_rrna):
        print('skip rRNA prediction...')
    else:
        print('predict rRNAs...')
        log.debug('start rRNA prediction')
        genome['features'][bc.FEATURE_R_RNA] = r_rna.predict_r_rnas(genome, contigs_path)
        print(f"\tfound: {len(genome['features'][bc.FEATURE_R_RNA])}")

    ############################################################################
    # ncRNA gene prediction
    ############################################################################
    if(cfg.skip_ncrna):
        print('skip ncRNA prediction...')
    else:
        print('predict ncRNAs...')
        log.debug('start ncRNA prediction')
        genome['features'][bc.FEATURE_NC_RNA] = nc_rna.predict_nc_rnas(genome, contigs_path)
        print(f"\tfound: {len(genome['features'][bc.FEATURE_NC_RNA])}")

    ############################################################################
    # ncRNA region prediction
    ############################################################################
    if(cfg.skip_ncrna_region):
        print('skip ncRNA region prediction...')
    else:
        print('predict ncRNA regions...')
        log.debug('start ncRNA region prediction')
        genome['features'][bc.FEATURE_NC_RNA_REGION] = nc_rna_region.predict_nc_rna_regions(genome, contigs_path)
        print(f"\tfound: {len(genome['features'][bc.FEATURE_NC_RNA_REGION])}")

    ############################################################################
    # CRISPR prediction
    ############################################################################
    if(cfg.skip_crispr):
        print('skip CRISPR array prediction...')
    else:
        print('predict CRISPR arrays...')
        log.debug('start CRISPR prediction')
        genome['features'][bc.FEATURE_CRISPR] = crispr.predict_crispr(genome, contigs_path)
        print(f"\tfound: {len(genome['features'][bc.FEATURE_CRISPR])}")

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
        print('predict & annotate CDSs...')
        log.debug('predict CDS')
        genome['features'][bc.FEATURE_CDS] = cds.predict(genome, contigs_path)
        print(f"\tpredicted: {len(genome['features'][bc.FEATURE_CDS])} ")
       
        log.debug('detect spurious CDS')
        discarded_cds = orf.detect_spurious(genome['features'][bc.FEATURE_CDS])
        print(f'\tdiscarded spurious: {len(discarded_cds)}')
        
        log.debug('lookup CDS UPS/IPS')
        cdss_ups, cdss_not_found = ups.lookup(genome['features'][bc.FEATURE_CDS])
        cdss_ips, tmp = ips.lookup(cdss_ups)
        cdss_not_found.extend(tmp)
        print(f'\tdetected IPSs: {len(cdss_ips)}')
        
        if(len(cdss_not_found) > 0):
            log.debug('search CDS PSC')
            cdss_psc, cdss_not_found = psc.search(cdss_not_found)
            print(f'\tfound PSCs: {len(cdss_psc)}')
        
        print(f"\tlookup annotations...")
        log.debug('lookup CDS PSCs')
        psc.lookup(genome['features'][bc.FEATURE_CDS])  # lookup PSC info
        
        print(f"\tmark hypotheticals and combine annotations...")
        log.debug('mark CDS hypotheticals')
        cds.mark_hypotheticals(genome['features'][bc.FEATURE_CDS])  # mark hypotheticals
        for feat in genome['features'][bc.FEATURE_CDS]:
            log.debug('combine CDS annotations')
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
        log.debug('predict sORF')
        sorfs = s_orf.extract(genome)
        print(f'\tpotential: {len(sorfs)}')

        log.debug('apply sORF overlap filter')
        sorfs, discarded_sorfs = s_orf.overlap_filter(genome, sorfs)
        print(f'\tdiscarded due to overlaps: {len(discarded_sorfs)}')

        log.debug('detect spurious sORF')
        discarded_sorfs = orf.detect_spurious(sorfs)
        print(f'\tdiscarded spurious: {len(discarded_sorfs)}')

        log.debug('lookup sORF UPS/IPS')
        sorf_upss, sorfs_not_found = ups.lookup(sorfs)
        sorf_ipss, tmp = ips.lookup(sorf_upss)
        sorfs_not_found.extend(tmp)
        print(f'\tdetected IPSs: {len(sorf_ipss)}')

        sorf_pscs = []
        if(len(sorfs_not_found) > 0):
            log.debug('search sORF PSC')
            tmp, sorfs_not_found = s_orf.search_pscs(sorfs_not_found)
            sorf_pscs.extend(tmp)
            print(f'\tfound PSCs: {len(sorf_pscs)}')

        print("lookup PSC annotations for sORFs...")
        log.debug('lookup sORF PSCs')
        sorf_pscs.extend(sorf_ipss)
        psc.lookup(sorf_pscs)  # lookup PSC info
        log.debug('filter sORF by annotations')
        sorfs_filtered = s_orf.annotation_filter(sorfs)
        for feat in sorfs_filtered:
            log.debug('combine sORF annotations')
            anno.combine_ips_psc_annotation(feat) # combine IPS and PSC annotations
        genome['features'][bc.FEATURE_SORF] = sorfs_filtered
        print(f'\tfiltered sORFs: {len(sorfs_filtered)}')
    
    ############################################################################
    # gap annotation
    # - in-mem gap detection
    # - gap annotation
    ############################################################################
    if(cfg.skip_gap):
        print('skip gap annotation...')
    else:
        print('detect gaps...')
        log.debug('detect gaps')
        assembly_gaps = gaps.detect_assembly_gaps(genome)
        genome['features'][bc.FEATURE_GAP] = assembly_gaps
        print(f'\tfound: {len(assembly_gaps)}')
    
    ############################################################################
    # oriC/T prediction
    ############################################################################
    if(cfg.skip_ori):
        print('skip oriC/T annotation...')
    else:
        print('detect oriCs/oriVs...')
        log.debug('detect oriC/V')
        oriCs = ori.predict_oris(genome, contigs_path, bc.FEATURE_ORIC)
        genome['features'][bc.FEATURE_ORIC] = oriCs
        print(f'\tfound: {len(oriCs)}')

        print('detect oriTs...')
        log.debug('detect oriT')
        oriTs = ori.predict_oris(genome, contigs_path, bc.FEATURE_ORIT)
        genome['features'][bc.FEATURE_ORIT] = oriTs
        print(f'\tfound: {len(oriTs)}')

    ############################################################################
    # Filter overlapping features
    ############################################################################
    print('apply feature overlap filters...')
    anno.detect_feature_overlaps(genome)

    ############################################################################
    # Create annotations
    # - filter features based on precedence and overlaps
    # - sort features
    # - create locus tags for features
    ############################################################################
    print('select features and create locus tags...')
    log.debug('start feature selection and creation of locus tags')
    features_by_contig = {k['id']: [] for k in genome['contigs']}
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
            bc.FEATURE_ORIV,
            bc.FEATURE_ORIT
        ]:
        feature_list = genome['features'].get(feature_type, [])
        for feature in feature_list:
            if('discarded' not in feature):
                contig_features = features_by_contig.get(feature['contig'])
                contig_features.append(feature)
    features = []
    for contig in genome['contigs']:
        contig_features = features_by_contig[contig['id']]
        contig_features.sort(key=lambda k: k['start'])
        features.extend(contig_features)
    log.info('selected features: %i', len(features))
    print(f'selected: {len(features)}')

    locus_tag_nr = 5
    # use user provided locus tag if not None/non-empty or generate a sequence based locus prefix
    locus_tag_prefix = cfg.locus_tag if cfg.locus_tag else bu.create_locus_tag_prefix(contigs)
    log.info('locus tag prefix: %s', locus_tag_prefix)
    for feature in features:
        locus_tag = f'{locus_tag_prefix}_{locus_tag_nr:05}'
        if(feature['type'] != bc.FEATURE_GAP and feature['type'] != bc.FEATURE_CRISPR and feature['type'] != bc.FEATURE_ORIC and feature['type'] != bc.FEATURE_ORIT):
            feature['locus'] = locus_tag
            locus_tag_nr += 5

    ############################################################################
    # Print statistics
    # - genome stats
    # - annotation stats
    ############################################################################
    print('\ngenome statistics:')
    genome_stats = bu.calc_genome_stats(genome, features)
    print(f"\tGenome size: {genome['size']:,} bp")
    print(f"\tContigs/replicons: {len(genome['contigs'])}")
    print(f"\tGC: {100 * genome_stats['gc']:.1f} %")
    print(f"\tN50: {genome_stats['n50']:,}")
    print(f"\tN ratio: {100 * genome_stats['n_ratio']:.1f} %")
    print(f"\tcoding density: {100 * genome_stats['coding_ratio']:.1f} %")

    print('\nannotation statistics:')
    print(f"\ttRNAs: {len([f for f in features if f['type'] == bc.FEATURE_T_RNA])}")
    print(f"\ttmRNAs: {len([f for f in features if f['type'] == bc.FEATURE_TM_RNA])}")
    print(f"\trRNAs: {len([f for f in features if f['type'] == bc.FEATURE_R_RNA])}")
    print(f"\tncRNAs: {len([f for f in features if f['type'] == bc.FEATURE_NC_RNA])}")
    print(f"\tncRNA regions: {len([f for f in features if f['type'] == bc.FEATURE_NC_RNA_REGION])}")
    print(f"\tCRISPR arrays: {len([f for f in features if f['type'] == bc.FEATURE_CRISPR])}")
    print(f"\tCDSs: {len([f for f in features if f['type'] == bc.FEATURE_CDS])}")
    print(f"\tsORFs: {len([f for f in features if f['type'] == bc.FEATURE_SORF])}")
    print(f"\tgaps: {len([f for f in features if f['type'] == bc.FEATURE_GAP])}")
    print(f"\toriCs/oriVs: {len([f for f in features if f['type'] == bc.FEATURE_ORIC])}")
    print(f"\toriTs: {len([f for f in features if f['type'] == bc.FEATURE_ORIT])}")

    ############################################################################
    # Write output files
    # - write comprehensive annotation results as JSON
    # - write optional output files in GFF3/GenBank/EMBL formats
    # - remove temp directory
    ############################################################################
    print('\nwrite JSON output...')
    log.debug('write JSON output')

    log.info(f'file prefix: {prefix}')
    json_path = cfg.output_path.joinpath(f'{prefix}.json')
    json.write_json(genome, features, json_path)

    if(cfg.tsv):
        print('write TSV output...')
        log.debug('write tsv output')
        tsv_path = cfg.output_path.joinpath(f'{prefix}.tsv')
        tsv.write_tsv(genome['contigs'], features_by_contig, tsv_path)

    if(cfg.gff3):
        print('write GFF3 output...')
        log.debug('write GFF3 output')
        gff3_path = cfg.output_path.joinpath(f'{prefix}.gff3')
        gff.write_gff3(genome['contigs'], features_by_contig, gff3_path)

    if(cfg.genbank):
        print('write GenBank output...')
        log.debug('write GenBank output')
        genbank_path = cfg.output_path.joinpath(f'{prefix}.gbff')
        genbank.write_genbank(features, genbank_path)

    if(cfg.fna):
        print('write genome sequences...')
        log.debug('write genome sequence output')
        fna_path = cfg.output_path.joinpath(f'{prefix}.fna')
        fasta.export_contigs(genome['contigs'], fna_path, description=True, wrap=True)

    if(cfg.faa):
        print('write translated CDS sequences...')
        log.debug('write translated CDS output')
        faa_path = cfg.output_path.joinpath(f'{prefix}.faa')
        fasta.write_faa(features, faa_path)

    # remove tmp dir
    shutil.rmtree(str(cfg.tmp_path))
    log.debug('removed tmp dir: %s', cfg.tmp_path)


if __name__ == '__main__':
    main()
