
import logging

import bakta.constants as bc

log = logging.getLogger('ANNOTATION')


def combine_ips_psc_annotation(feature):
    ups = feature.get('ups', None)
    ips = feature.get('ips', None)
    psc = feature.get('psc', None)

    gene = ''
    product = bc.HYPOTHETICAL_PROTEIN
    db_xrefs = set()
    if(psc):
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
    feature['gene'] = gene
    feature['product'] = product
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
                        "overlap: sORF (%s/%s) [%i, %i] overlapping tRNA (%s) [%i, %i], %s, contig=%s",
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
                            "discard overlapping sORF: (%s/%s) [%i, %i] overlapping sORF (%s/%s) [%i, %i], %s, contig=%s, lower annotation score (%i/%i)",
                            sorf.get('gene', '-'), sorf.get('product', '-'), sorf['start'], sorf['stop'], overlap_sorf.get('gene', '-'), overlap_sorf.get('product', '-'), overlap_sorf['start'], overlap_sorf['stop'], overlap, sorf['contig'], score_sorf, score_overlap_sorf
                        )
                    elif(score_sorf == score_overlap_sorf and len(sorf['sequence']) < len(overlap_sorf['sequence'])):  # equal annotation score but shorter sequence -> potential fragment or too short ORF prediction
                        overlap = f"[{max(sorf['start'], overlap_sorf['start'])},{min(sorf['stop'], overlap_sorf['stop'])}]"
                        sorf['discarded'] = {
                            'type': bc.DISCARD_TYPE_OVERLAP,
                            'feature_type': bc.FEATURE_SORF,
                            'description': f"overlaps {bc.FEATURE_SORF} ({overlap_sorf.get('gene', '-')}/{overlap_sorf.get('product', '-')}) at {overlap} with equal score ({score_sorf}) but lower length ({len(sorf['sequence'])}/{len(overlap_sorf['sequence'])})"
                        }
                        log.info(
                            "discard overlapping sORF: (%s/%s) [%i, %i] overlapping sORF (%s/%s) [%i, %i], %s, contig=%s, equal annotation score (%i), lower length (%i/%i)",
                            sorf.get('gene', '-'), sorf.get('product', '-'), sorf['start'], sorf['stop'], overlap_sorf.get('gene', '-'), overlap_sorf.get('product', '-'), overlap_sorf['start'], overlap_sorf['stop'], overlap, sorf['contig'], score_sorf, len(sorf['sequence']), len(overlap_sorf['sequence'])
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
    return score
