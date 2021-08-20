
import logging
import re

import bakta.config as cfg
import bakta.constants as bc
import bakta.io.insdc as insdc


log = logging.getLogger('ANNOTATION')


RE_CONTIG = re.compile(r'\s+contig\s*', flags=re.IGNORECASE)
RE_HOMOLOG = re.compile(r'\shomolog(?: (\d+))?', flags=re.IGNORECASE)
RE_PUTATIVE = re.compile(r'(potential|possible|probable|predicted)', flags=re.IGNORECASE)
RE_MULTIWHITESPACE = re.compile(r'\s{2,}')
RE_NODE = re.compile(r'NODE_', flags=re.IGNORECASE)
RE_POTENTIAL_CONTIG_NAME = re.compile(r'(genome|shotgun)', flags=re.IGNORECASE)
RE_NO_LETTERS = re.compile(r'[^A-Za-z]')


def combine_annotation(feature):
    ups = feature.get('ups', None)
    ips = feature.get('ips', None)
    psc = feature.get('psc', None)
    pscc = feature.get('pscc', None)
    expert = feature.get('expert', None)

    gene = None
    product = None
    db_xrefs = set()
    if(pscc):
        pscc_product = pscc.get('product', None)
        if(pscc_product):
            product = pscc_product
        for db_xref in pscc['db_xrefs']:
            db_xrefs.add(db_xref)
    if(psc and psc.get('valid', True)):
        psc_gene = psc.get('gene', None)
        if(psc_gene):
            gene = psc_gene
        psc_product = psc.get('product', None)
        if(psc_product):
            product = psc_product
        for db_xref in psc['db_xrefs']:
            db_xrefs.add(db_xref)
    if(ups):
        for db_xref in ups['db_xrefs']:
            db_xrefs.add(db_xref)
    if(ips):
        ips_gene = ips.get('gene', None)
        if(ips_gene):
            gene = ips_gene
        ips_product = ips.get('product', None)
        if(ips_product):
            product = ips_product
        for db_xref in ips['db_xrefs']:
            db_xrefs.add(db_xref)
    if(expert):
        rank = 0
        for expert_system, expert_hit in expert.items():
            expert_rank = expert_hit['rank']
            if(expert_rank > rank):
                expert_gene = expert_hit.get('gene', None)
                if(expert_gene):
                    gene = expert_gene
                expert_product = expert_hit.get('product', None)
                if(expert_product):
                    product = expert_product
                expert_db_xrefs = expert_hit.get('db_xrefs', None)
                if(expert_db_xrefs):
                    for expert_db_xref in expert_db_xrefs:
                        db_xrefs.add(expert_db_xref)
                rank = expert_rank

    if(gene):
        feature['gene'] = gene

    if(product):
        feature['product'] = product
        revise_cds_product(feature)
        if(cfg.compliant and not feature.get('hypothetical', False)):
            insdc.revise_product_insdc(feature)
    else:
        feature['product'] = bc.HYPOTHETICAL_PROTEIN
        feature['hypothetical'] = True

    if(gene and feature.get('hypothetical', False)):
        feature['gene'] = None
        log.info(
            'remove gene symbol of hypothetical protein: contig=%s, start=%i, stop=%i, strand=%s, gene=%s',
            feature['contig'], feature['start'], feature['stop'], feature['strand'], gene
        )

    feature['db_xrefs'] = sorted(list(db_xrefs))


def detect_feature_overlaps(genome):
    """Apply feature type specific hierarchical feature overlap filters.
    tRNA < tmRNA
    CDS < tmRNA, tRNA, rRNA, CRISPR
    sORF < mRNA, tRNA, rRNA, CRISPR, CDS (in-frame & entirely overlapping), sORF (shorter, weaker annotations)
    """
    contig_t_rnas = {k['id']: [] for k in genome['contigs']}
    for t_rna in genome['features'].get(bc.FEATURE_T_RNA, []):
        t_rnas = contig_t_rnas[t_rna['contig']]
        t_rnas.append(t_rna)
    contig_tm_rnas = {k['id']: [] for k in genome['contigs']}
    for tm_rna in genome['features'].get(bc.FEATURE_TM_RNA, []):
        tm_rnas = contig_tm_rnas[tm_rna['contig']]
        tm_rnas.append(tm_rna)
    contig_r_rnas = {k['id']: [] for k in genome['contigs']}
    for r_rna in genome['features'].get(bc.FEATURE_R_RNA, []):
        r_rnas = contig_r_rnas[r_rna['contig']]
        r_rnas.append(r_rna)
    contig_ncrna_regions = {k['id']: [] for k in genome['contigs']}
    for ncRNA_region in genome['features'].get(bc.FEATURE_NC_RNA_REGION, []):
        ncRNA_regions = contig_ncrna_regions[ncRNA_region['contig']]
        ncRNA_regions.append(ncRNA_region)
    contig_crispr_arrays = {k['id']: [] for k in genome['contigs']}
    for crispr_array in genome['features'].get(bc.FEATURE_CRISPR, []):
        crispr_arrays = contig_crispr_arrays[crispr_array['contig']]
        crispr_arrays.append(crispr_array)
    contig_cdss = {k['id']: [] for k in genome['contigs']}
    for cds in genome['features'].get(bc.FEATURE_CDS, []):
        cdss = contig_cdss[cds['contig']]
        cdss.append(cds)
    contig_sorfs = {k['id']: [] for k in genome['contigs']}
    for sorf in genome['features'].get(bc.FEATURE_SORF, []):
        sorfs = contig_sorfs[sorf['contig']]
        sorfs.append(sorf)

    for contig in genome['contigs']:  # find feature overlaps contig-wise to increase the performance
        log.debug('filter features on contig: %s', contig['id'])

        # mark tRNAs overlapping with tmRNAs
        for tRNA in contig_t_rnas[contig['id']]:
            for tmRNA in contig_tm_rnas[contig['id']]:
                if(tRNA['stop'] < tmRNA['start'] or tRNA['start'] > tmRNA['stop']):
                    continue
                else:  # overlap -> remove tRNA
                    overlap = f"[{max(tRNA['start'], tmRNA['start'])},{min(tRNA['stop'], tmRNA['stop'])}]"
                    tRNA['discarded'] = {
                        'type': bc.DISCARD_TYPE_OVERLAP,
                        'feature_type': bc.FEATURE_TM_RNA,
                        'description': f"{bc.FEATURE_TM_RNA} overlap with ({tmRNA['product']}) at {overlap}"
                    }
                    log.info(
                        "overlap: tRNA (%s) [%i, %i] overlapping with tmRNA (%s) [%i, %i] at %s on contig=%s",
                        tRNA['product'], tRNA['start'], tRNA['stop'], tmRNA['product'], tmRNA['start'], tmRNA['stop'], overlap, tRNA['contig']
                    )

        # mark ncRNA-regions overlapping with ncRNA-regions
        for ncRNA_region in contig_ncrna_regions[contig['id']]:
            for ncRNA_region_overlap in contig_ncrna_regions[contig['id']]:
                if(ncRNA_region['stop'] < ncRNA_region_overlap['start'] or ncRNA_region['start'] > ncRNA_region_overlap['stop']):
                    continue
                if(ncRNA_region['db_xrefs'][0] == ncRNA_region_overlap['db_xrefs'][0]):
                    continue  # same
                else:  # overlap -> select ncRNA based on bitscore, discard the other
                    if(ncRNA_region['score'] < ncRNA_region_overlap['score']):
                        overlap = f"[{max(ncRNA_region['start'], ncRNA_region_overlap['start'])},{min(ncRNA_region['stop'], ncRNA_region_overlap['stop'])}]"
                        ncRNA_region['discarded'] = {
                            'type': bc.DISCARD_TYPE_OVERLAP,
                            'feature_type': bc.FEATURE_NC_RNA_REGION,
                            'description': f"{bc.FEATURE_NC_RNA_REGION} overlap with ({ncRNA_region_overlap['product']}) at {overlap}"
                        }
                        log.info(
                            "overlap: ncRNA-region (%s) [%i, %i] overlapping with ncRNA-region (%s) [%i, %i] at %s on contig=%s, lower bitscore (%f/%f)",
                            ncRNA_region['product'], ncRNA_region['start'], ncRNA_region['stop'], ncRNA_region_overlap['product'], ncRNA_region_overlap['start'], ncRNA_region_overlap['stop'], overlap, ncRNA_region['contig'], ncRNA_region['score'], ncRNA_region_overlap['score']
                        )

        # mark CDS overlapping with tRNAs, tmRNAs, rRNAs, CRISPRs
        for cds in contig_cdss[contig['id']]:
            # tmRNA overlaps
            for tmRNA in contig_tm_rnas[contig['id']]:
                if(cds['stop'] < tmRNA['start'] or cds['start'] > tmRNA['stop']):
                    continue
                else:  # overlap -> remove cds
                    overlap = f"[{max(cds['start'], tmRNA['start'])},{min(cds['stop'], tmRNA['stop'])}]"
                    cds['discarded'] = {
                        'type': bc.DISCARD_TYPE_OVERLAP,
                        'feature_type': bc.FEATURE_TM_RNA,
                        'description': f"{bc.FEATURE_TM_RNA} overlap with ({tmRNA['product']}) at {overlap}"
                    }
                    log.info(
                        "overlap: CDS (%s/%s) [%i, %i] overlapping tmRNA (%s) [%i, %i], %s, contig=%s",
                        cds.get('gene', '-'), cds.get('product', '-'), cds['start'], cds['stop'], tmRNA['gene'], tmRNA['start'], tmRNA['stop'], overlap, cds['contig']
                    )
            # tRNA overlaps
            for tRNA in contig_t_rnas[contig['id']]:
                if(cds['stop'] < tRNA['start'] or cds['start'] > tRNA['stop']):
                    continue
                else:  # overlap -> remove cds
                    overlap = f"[{max(cds['start'], tRNA['start'])},{min(cds['stop'], tRNA['stop'])}]"
                    cds['discarded'] = {
                        'type': bc.DISCARD_TYPE_OVERLAP,
                        'feature_type': bc.FEATURE_T_RNA,
                        'description': f"{bc.FEATURE_T_RNA} overlap with ({tRNA['product']}) at {overlap}"
                    }
                    log.info(
                        "overlap: CDS (%s/%s) [%i, %i] overlapping tRNA (%s) [%i, %i], %s, contig=%s",
                        cds.get('gene', '-'), cds.get('product', '-'), cds['start'], cds['stop'], tRNA['gene'], tRNA['start'], tRNA['stop'], overlap, cds['contig']
                    )
            # rRNA overlaps
            for rRNA in contig_r_rnas[contig['id']]:
                if(cds['stop'] < rRNA['start'] or cds['start'] > rRNA['stop']):
                    continue
                else:  # overlap -> remove cds
                    overlap = f"[{max(cds['start'], rRNA['start'])},{min(cds['stop'], rRNA['stop'])}]"
                    cds['discarded'] = {
                        'type': bc.DISCARD_TYPE_OVERLAP,
                        'feature_type': bc.FEATURE_R_RNA,
                        'description': f"{bc.FEATURE_R_RNA} overlap with ({rRNA['product']}) at {overlap}"
                    }
                    log.info(
                        "overlap: CDS (%s/%s) [%i, %i] overlapping rRNA (%s) [%i, %i], %s, contig=%s",
                        cds.get('gene', '-'), cds.get('product', '-'), cds['start'], cds['stop'], rRNA['gene'], rRNA['start'], rRNA['stop'], overlap, cds['contig']
                    )
            # CRISPR overlaps
            for crispr in contig_crispr_arrays[contig['id']]:
                if(cds['stop'] < crispr['start'] or cds['start'] > crispr['stop']):
                    continue
                else:  # overlap -> remove cds
                    overlap = f"[{max(cds['start'], crispr['start'])},{min(cds['stop'], crispr['stop'])}]"
                    cds['discarded'] = {
                        'type': bc.DISCARD_TYPE_OVERLAP,
                        'feature_type': bc.FEATURE_CRISPR,
                        'description': f'overlaps {bc.FEATURE_CRISPR} at {overlap}'
                    }
                    log.info(
                        "overlap: CDS (%s/%s) [%i, %i] overlapping CRISPR [%i, %i], %s, contig=%s",
                        cds.get('gene', '-'), cds.get('product', '-'), cds['start'], cds['stop'], crispr['start'], crispr['stop'], overlap, cds['contig']
                    )

        # remove sORF overlapping with tRNAs, tmRNAs, rRNAs, CRISPRs, inframe CDSs, shorter inframe sORFs
        for sorf in contig_sorfs[contig['id']]:
            # tmRNA overlaps
            for tmRNA in contig_tm_rnas[contig['id']]:
                if(sorf['stop'] < tmRNA['start'] or sorf['start'] > tmRNA['stop']):
                    continue
                else:  # overlap -> remove sorf
                    overlap = f"[{max(sorf['start'], tmRNA['start'])},{min(sorf['stop'], tmRNA['stop'])}]"
                    sorf['discarded'] = {
                        'type': bc.DISCARD_TYPE_OVERLAP,
                        'feature_type': bc.FEATURE_TM_RNA,
                        'description': f"{bc.FEATURE_TM_RNA} overlap with ({tmRNA['product']}) at {overlap}"
                    }
                    log.info(
                        "overlap: sORF (%s/%s) [%i, %i] overlapping tmRNA (%s) [%i, %i], %s, contig=%s",
                        sorf.get('gene', '-'), sorf.get('product', '-'), sorf['start'], sorf['stop'], tmRNA['gene'], tmRNA['start'], tmRNA['stop'], overlap, sorf['contig']
                    )
            # tRNA overlaps
            for tRNA in contig_t_rnas[contig['id']]:
                if(sorf['stop'] < tRNA['start'] or sorf['start'] > tRNA['stop']):
                    continue
                else:  # overlap -> remove sorf
                    overlap = f"[{max(sorf['start'], tRNA['start'])},{min(sorf['stop'], tRNA['stop'])}]"
                    sorf['discarded'] = {
                        'type': bc.DISCARD_TYPE_OVERLAP,
                        'feature_type': bc.FEATURE_T_RNA,
                        'description': f"{bc.FEATURE_T_RNA} overlap with ({tRNA['product']}) at {overlap}"
                    }
                    log.info(
                        "overlap: sORF (%s/%s) [%i, %i] overlapping tRNA (%s) [%i, %i], %s, contig=%s",
                        sorf.get('gene', '-'), sorf.get('product', '-'), sorf['start'], sorf['stop'], tRNA['gene'], tRNA['start'], tRNA['stop'], overlap, sorf['contig']
                    )
            # rRNA overlaps
            for rRNA in contig_r_rnas[contig['id']]:
                if(sorf['stop'] < rRNA['start'] or sorf['start'] > rRNA['stop']):
                    continue
                else:  # overlap -> remove sorf
                    overlap = f"[{max(sorf['start'], rRNA['start'])},{min(sorf['stop'], rRNA['stop'])}]"
                    sorf['discarded'] = {
                        'type': bc.DISCARD_TYPE_OVERLAP,
                        'feature_type': bc.FEATURE_R_RNA,
                        'description': f"{bc.FEATURE_R_RNA} overlap with ({rRNA['product']}) at {overlap}"
                    }
                    log.info(
                        "overlap: sORF (%s/%s) [%i, %i] overlapping rRNA (%s) [%i, %i], %s, contig=%s",
                        sorf.get('gene', '-'), sorf.get('product', '-'), sorf['start'], sorf['stop'], rRNA['gene'], rRNA['start'], rRNA['stop'], overlap, sorf['contig']
                    )
            # CRISPR overlaps
            for crispr in contig_crispr_arrays[contig['id']]:
                if(sorf['stop'] < crispr['start'] or sorf['start'] > crispr['stop']):
                    continue
                else:  # overlap -> remove sorf
                    overlap = f"[{max(sorf['start'], crispr['start'])},{min(sorf['stop'], crispr['stop'])}]"
                    sorf['discarded'] = {
                        'type': bc.DISCARD_TYPE_OVERLAP,
                        'feature_type': bc.FEATURE_CRISPR,
                        'description': f'overlaps {bc.FEATURE_CRISPR} at {overlap}'
                    }
                    log.info(
                        "overlap: sORF (%s/%s) [%i, %i] overlapping CRISPR [%i, %i], %s, contig=%s",
                        sorf.get('gene', '-'), sorf.get('product', '-'), sorf['start'], sorf['stop'], crispr['start'], crispr['stop'], overlap, sorf['contig']
                    )
            # CDS overlaps (skipped as most overlaps are filtered in sORF detection)

            # sORF overlaps
            for overlap_sorf in contig_sorfs[contig['id']]:
                if(sorf['stop'] < overlap_sorf['start'] or sorf['start'] > overlap_sorf['stop']):
                    continue  # no overlap
                elif(sorf['start'] == overlap_sorf['start'] and sorf['stop'] == overlap_sorf['stop']):
                    continue  # same
                else:  # overlap -> remove sorf
                    score_sorf = calc_sorf_annotation_score(sorf)
                    score_overlap_sorf = calc_sorf_annotation_score(overlap_sorf)

                    if(score_sorf < score_overlap_sorf):  # lower annotation score
                        overlap = f"[{max(sorf['start'], overlap_sorf['start'])},{min(sorf['stop'], overlap_sorf['stop'])}]"
                        sorf['discarded'] = {
                            'type': bc.DISCARD_TYPE_OVERLAP,
                            'feature_type': bc.FEATURE_SORF,
                            'description': f"overlaps {bc.FEATURE_SORF} ({overlap_sorf.get('gene', '-')}/{overlap_sorf.get('product', '-')}) at {overlap} with lower score ({score_sorf}/{score_overlap_sorf})"
                        }
                        log.info(
                            "overlap: sORF (%s/%s) [%i, %i] overlapping sORF (%s/%s) [%i, %i], %s, contig=%s, lower annotation score (%i/%i)",
                            sorf.get('gene', '-'), sorf.get('product', '-'), sorf['start'], sorf['stop'], overlap_sorf.get('gene', '-'), overlap_sorf.get('product', '-'), overlap_sorf['start'], overlap_sorf['stop'], overlap, sorf['contig'], score_sorf, score_overlap_sorf
                        )
                    elif(score_sorf == score_overlap_sorf and len(sorf['aa']) < len(overlap_sorf['aa'])):  # equal annotation score but shorter sequence -> potential fragment or too short ORF prediction
                        overlap = f"[{max(sorf['start'], overlap_sorf['start'])},{min(sorf['stop'], overlap_sorf['stop'])}]"
                        sorf['discarded'] = {
                            'type': bc.DISCARD_TYPE_OVERLAP,
                            'feature_type': bc.FEATURE_SORF,
                            'description': f"overlaps {bc.FEATURE_SORF} ({overlap_sorf.get('gene', '-')}/{overlap_sorf.get('product', '-')}) at {overlap} with equal score ({score_sorf}) but lower length ({len(sorf['aa'])}/{len(overlap_sorf['aa'])})"
                        }
                        log.info(
                            "overlap: sORF (%s/%s) [%i, %i] overlapping sORF (%s/%s) [%i, %i], %s, contig=%s, equal annotation score (%i), lower length (%i/%i)",
                            sorf.get('gene', '-'), sorf.get('product', '-'), sorf['start'], sorf['stop'], overlap_sorf.get('gene', '-'), overlap_sorf.get('product', '-'), overlap_sorf['start'], overlap_sorf['stop'], overlap, sorf['contig'], score_sorf, len(sorf['aa']), len(overlap_sorf['aa'])
                        )


def calc_sorf_annotation_score(sorf):
    """Calc an annotation score rewarding each identification & annotation"""
    score = 0

    if('ups' in sorf):
        score += 1

    ips = sorf.get('ips', None)
    if(ips):
        score += 1
        ips_gene = ips.get('gene', None)
        if(ips_gene):
            score += 1
        ips_product = ips.get('product', None)
        if(ips_product):
            score += 1

    psc = sorf.get('psc', None)
    if(psc):
        score += 1
        psc_gene = psc.get('gene', None)
        if(psc_gene):
            score += 1
        psc_product = psc.get('product', None)
        if(psc_product):
            score += 1
    log.debug(
        'sorf score: contig=%s, start=%i, stop=%i, gene=%s, product=%s, score=%i',
        sorf['contig'], sorf['start'], sorf['stop'], sorf.get('gene', '-'), sorf.get('product', '-'), score
    )
    return score


def revise_cds_product(feature):
    """Revise product name for INSDC compliant submissions"""
    product = feature['product']

    old_product = product
    if('.' in product):  # remove periods
        product = product.replace('.', '')
        log.info('fix product: replace periods. new=%s, old=%s', product, old_product)

    old_product = product
    if('=' in product):  # remove equal chars
        product = product.replace('=', '')
        log.info('fix product: remove = character. new=%s, old=%s', product, old_product)

    old_product = product
    if(RE_CONTIG.search(product)):
        product = bc.HYPOTHETICAL_PROTEIN
        feature['hypothetical'] = True
        log.info('fix product: remove product containing "contig". new=%s, old=%s', product, old_product)

    old_product = product
    if(RE_HOMOLOG.search(product)):  # replace Homologs
        product = RE_HOMOLOG.sub('-like protein', product)
        log.info('fix product: replace Homolog. new=%s, old=%s', product, old_product)

    old_product = product
    if(RE_PUTATIVE.search(product)):  # replace putative synonyms)
        product = RE_PUTATIVE.sub('putative', product)
        log.info('fix product: replace putative synonyms. new=%s, old=%s', product, old_product)

    old_product = product
    if(RE_MULTIWHITESPACE.search(product)):  # squeeze multiple whitespaces
        product = RE_MULTIWHITESPACE.sub(' ', product)
        log.info('fix product: squeeze multiple whitespaces. new=%s, old=%s', product, old_product)

    old_product = product
    product = product.strip()  # trim leading/trailing whitespaces
    if(product != old_product):
        log.info('fix product: trim leading/trailing whitespace. new=%s, old=%s', product, old_product)

    if(
        RE_NODE.search(product) or  # potential contig name (SPAdes)
        RE_POTENTIAL_CONTIG_NAME.search(product) or  # potential contig name (SPAdes)
        RE_NO_LETTERS.fullmatch(product)  # no letters -> set to Hypothetical
        ):  # remove suspect products and mark as hypothetical
        product = bc.HYPOTHETICAL_PROTEIN
        feature['hypothetical'] = True
        log.info('remove product: mark proteins with suspect products as hypothetical. old=%s', old_product)

    feature['product'] = product
