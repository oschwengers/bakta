import atexit
import logging
import os
import shutil
import sys
import platform as pf

from datetime import datetime
from pathlib import Path

import bakta
import bakta.constants as bc
import bakta.config as cfg
import bakta.io.fasta as fasta
import bakta.io.json as json
import bakta.io.tsv as tsv
import bakta.io.gff as gff
import bakta.io.insdc as insdc
import bakta.expert.amrfinder as exp_amr
import bakta.expert.protein_sequences as exp_aa_seq
import bakta.features.annotation as anno
import bakta.features.t_rna as t_rna
import bakta.features.tm_rna as tm_rna
import bakta.features.r_rna as r_rna
import bakta.features.nc_rna as nc_rna
import bakta.features.nc_rna_region as nc_rna_region
import bakta.features.crispr as crispr
import bakta.features.orf as orf
import bakta.features.cds as feat_cds
import bakta.features.s_orf as s_orf
import bakta.features.gaps as gaps
import bakta.features.ori as ori
import bakta.db as db
import bakta.utils as bu
import bakta.ups as ups
import bakta.ips as ips
import bakta.psc as psc
import bakta.pscc as pscc


def main():
    args = bu.parse_arguments()  # parse arguments

    ############################################################################
    # Setup logging
    ############################################################################
    cfg.prefix = args.prefix if args.prefix else Path(args.genome).stem
    try:
        output_path = Path(args.output) if args.output else Path.cwd()
        if(not output_path.exists()):
            output_path.mkdir(parents=True, exist_ok=True)
        elif(not os.access(str(output_path), os.X_OK)):
            sys.exit(f'ERROR: output path ({output_path}) not accessible!')
        elif(not os.access(str(output_path), os.W_OK)):
            sys.exit(f'ERROR: output path ({output_path}) not writable!')
        output_path = output_path.resolve()
        cfg.output_path = output_path
    except:
        sys.exit(f'ERROR: could not resolve or create output directory ({args.output})!')
    logging.basicConfig(
        filename=str(output_path.joinpath(f'{cfg.prefix}.log')),
        filemode='w',
        format='%(asctime)s.%(msecs)03d - %(levelname)s - %(name)s - %(message)s',
        datefmt='%H:%M:%S',
        level=logging.DEBUG if args.verbose else logging.INFO
    )
    log = logging.getLogger('MAIN')
    log.info('version=%s', bakta.__version__)
    log.info('developer: Oliver Schwengers, https://github.com/oschwengers')
    log.info('command: %s', ' '.join(sys.argv))
    log.info('local time: %s', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    log.info('machine: type=%s, cores=%s', pf.processor(), os.cpu_count())
    log.info('system: type=%s, release=%s', pf.system(), pf.release())
    log.info('python: version=%s, implementation=%s', pf.python_version(), pf.python_implementation())

    ############################################################################
    # Checks and configurations
    # - check parameters and setup global configuration
    # - test database
    # - test binary dependencies
    ############################################################################
    cfg.setup(args)  # check parameters and prepare global configuration
    atexit.register(cleanup, log, cfg.tmp_path)  # register cleanup exit hook
    cfg.db_info = db.check(cfg.db_path)
    bu.test_dependencies()
    if(cfg.verbose):
        print(f'Bakta v{bakta.__version__}')
        print('Options and arguments:')
        print(f'\tinput: {cfg.genome_path}')
        print(f"\tdb: {cfg.db_path}, version {cfg.db_info['major']}.{cfg.db_info['minor']}")
        print(f'\toutput: {cfg.output_path}')
        print(f'\tprefix: {cfg.prefix}')
        print(f'\ttmp directory: {cfg.tmp_path}')
        print(f'\t# threads: {cfg.threads}')
        if(cfg.compliant): print(f'\tINSDC compliant: {cfg.compliant}')
        if(cfg.keep_contig_headers): print(f'\tkeep contig headers: {cfg.keep_contig_headers}')
        if(cfg.locus): print(f'\tlocus prefix: {cfg.locus}')
        if(cfg.locus_tag): print(f'\tlocus tag prefix: {cfg.locus_tag}')
        if(cfg.complete): print(f'\tcomplete replicons: {cfg.complete}')
        if(cfg.prodigal_tf): print(f'\tprodigal training file: {cfg.prodigal_tf}')
        if(cfg.replicons): print(f'\treplicon table: {cfg.replicons}')
        print(f'\ttranslation table: {cfg.translation_table}')
        if(cfg.taxon): print(f'\ttaxon: {cfg.taxon}')
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
    contigs, complete_genome = bu.qc_contigs(contigs, replicons)
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
        'taxon': cfg.taxon,
        'gram': cfg.gram,
        'translation_table': cfg.translation_table,
        'size': sum([c['length'] for c in contigs]),
        'complete': cfg.complete or complete_genome,
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
    # - conduct expert systems analysis
    # - lookup & combine annotations
    # - analyze hypotheticals
    ############################################################################
    if(cfg.skip_cds):
        print('skip CDS prediction...')
    else:
        print('predict & annotate CDSs...')
        log.debug('predict CDS')
        cdss = feat_cds.predict(genome, contigs_path)
        print(f"\tpredicted: {len(cdss)} ")

        log.debug('detect spurious CDS')
        discarded_cdss = orf.detect_spurious(cdss) if len(cdss) > 0 else []
        print(f'\tdiscarded spurious: {len(discarded_cdss)}')
        cdss = [cds for cds in cdss if 'discarded' not in cds]

        log.debug('lookup CDS UPS/IPS')
        cdss_ups, cdss_not_found = ups.lookup(cdss)
        cdss_ips, tmp = ips.lookup(cdss_ups)
        cdss_not_found.extend(tmp)
        print(f'\tdetected IPSs: {len(cdss_ips)}')

        if(len(cdss_not_found) > 0):
            cds_fasta_path = cfg.tmp_path.joinpath('cds.unidentified.faa')
            with cds_fasta_path.open(mode='w') as fh:
                for cds in cdss_not_found:
                    fh.write(f">{cds['aa_hexdigest']}-{cds['contig']}-{cds['start']}\n{cds['aa']}\n")
            log.debug('search CDS PSC')
            cdss_psc, cdss_not_found = psc.search(cdss_not_found, cds_fasta_path)
            print(f'\tfound PSCs: {len(cdss_psc)}')
        print('\tlookup annotations...')
        log.debug('lookup CDS PSCs')
        psc.lookup(cdss)  # lookup PSC info
        pscc.lookup(cdss)  # lookup PSCC info

        print('\tconduct expert systems...')  # conduct expert systems annotation
        cds_fasta_path = cfg.tmp_path.joinpath('cds.faa')
        with cds_fasta_path.open(mode='w') as fh:
            for cds in cdss:
                fh.write(f">{cds['aa_hexdigest']}-{cds['contig']}-{cds['start']}\n{cds['aa']}\n")
        log.debug('conduct expert system: amrfinder')
        expert_amr_found = exp_amr.search(cdss, cds_fasta_path)
        print(f'\t\tamrfinder: {len(expert_amr_found)}')
        log.debug('conduct expert system: aa seqs')
        expert_aa_found = exp_aa_seq.search(cdss, cds_fasta_path)
        print(f'\t\tprotein sequences: {len(expert_aa_found)}')

        print('\tcombine annotations and mark hypotheticals...')
        log.debug('combine CDS annotations')
        for cds in cdss:
            anno.combine_annotation(cds)  # combine IPS & PSC annotations and mark hypotheticals

        log.debug('analyze hypotheticals')
        hypotheticals = [cds for cds in cdss if 'hypothetical' in cds]
        if(len(hypotheticals) > 0):
            print(f'\tanalyze hypothetical proteins: {len(hypotheticals)}')
            pfam_hits = feat_cds.predict_pfam(hypotheticals)
            print(f"\tdetected Pfam hits: {len(pfam_hits)} ")
            feat_cds.analyze_proteins(hypotheticals)
            print('\tcalculated proteins statistics')

        genome['features'][bc.FEATURE_CDS] = cdss

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
        discarded_sorfs = orf.detect_spurious(sorfs) if len(sorfs) > 0 else []
        print(f'\tdiscarded spurious: {len(discarded_sorfs)}')
        sorfs = [sorf for sorf in sorfs if 'discarded' not in sorf]

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

        print("\tlookup annotations...")
        log.debug('lookup sORF PSCs')
        sorf_pscs.extend(sorf_ipss)
        psc.lookup(sorf_pscs)  # lookup PSC info
        print('\tfilter and combine annotations...')
        log.debug('filter sORF by annotations')
        sorfs_filtered = s_orf.annotation_filter(sorfs)
        log.debug('combine sORF annotations')
        for feat in sorfs_filtered:
            anno.combine_annotation(feat)  # combine IPS and PSC annotations
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
    log.info('selected features=%i', len(features))
    print(f'selected: {len(features)}')

    locus_tag_nr = 5
    # use user provided locus tag if not None/non-empty or generate a sequence based locus prefix
    locus_tag_prefix = cfg.locus_tag if cfg.locus_tag else bu.create_locus_tag_prefix(contigs)
    log.info('locus tag prefix=%s', locus_tag_prefix)
    for feature in features:
        locus_tag = f'{locus_tag_prefix}_{locus_tag_nr:05}'
        if(feature['type'] not in [bc.FEATURE_NC_RNA_REGION, bc.FEATURE_GAP]):
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
    cdss = [f for f in features if f['type'] == bc.FEATURE_CDS]
    print(f"\tCDSs: {len(cdss)}, hypotheticals: {len([cds for cds in cdss if 'hypothetical' in cds])}")
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
    json_path = cfg.output_path.joinpath(f'{cfg.prefix}.json')
    json.write_json(genome, features, json_path)

    print('write TSV output...')
    tsv_path = cfg.output_path.joinpath(f'{cfg.prefix}.tsv')
    tsv.write_tsv(genome['contigs'], features_by_contig, tsv_path)

    print('write GFF3 output...')
    gff3_path = cfg.output_path.joinpath(f'{cfg.prefix}.gff3')
    gff.write_gff3(genome, features_by_contig, gff3_path)

    print('write INSDC (GenBank/EMBL) output...')
    genbank_path = cfg.output_path.joinpath(f'{cfg.prefix}.gbff')
    embl_path = cfg.output_path.joinpath(f'{cfg.prefix}.embl')
    insdc.write_insdc(genome, features, genbank_path, embl_path)

    print('write genome sequences...')
    fna_path = cfg.output_path.joinpath(f'{cfg.prefix}.fna')
    fasta.export_contigs(genome['contigs'], fna_path, description=True, wrap=True)

    print('write feature nucleotide sequences...')
    ffn_path = cfg.output_path.joinpath(f'{cfg.prefix}.ffn')
    fasta.write_ffn(features, ffn_path)

    print('write translated CDS sequences...')
    faa_path = cfg.output_path.joinpath(f'{cfg.prefix}.faa')
    fasta.write_faa(features, faa_path)

    if(cfg.skip_cds is False):
        hypotheticals = [feat for feat in features if feat['type'] == bc.FEATURE_CDS and 'hypothetical' in feat]
        print('write hypothetical TSV output...')
        tsv_path = cfg.output_path.joinpath(f'{cfg.prefix}.hypotheticals.tsv')
        tsv.write_hypothetical_tsv(hypotheticals, tsv_path)

        print('write translated hypothetical CDS sequences...')
        faa_path = cfg.output_path.joinpath(f'{cfg.prefix}.hypotheticals.faa')
        fasta.write_faa(hypotheticals, faa_path)


def cleanup(log, tmp_path):
    shutil.rmtree(str(tmp_path))  # remove tmp dir
    log.info('removed tmp dir: %s', tmp_path)


if __name__ == '__main__':
    main()
