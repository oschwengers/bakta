
import logging
import re

import bakta.config as cfg
import bakta.constants as bc
import bakta.io.insdc as insdc

from typing import Sequence

log = logging.getLogger('ANNOTATION')


RE_MULTIWHITESPACE = re.compile(r'\s{2,}')

RE_PROTEIN_CONTIG = re.compile(r'\s+contig\s*', flags=re.IGNORECASE)
RE_PROTEIN_HOMOLOG = re.compile(r'\shomolog(?: (\d+))?', flags=re.IGNORECASE)
RE_PROTEIN_PUTATIVE = re.compile(r'(potential|possible|probable|predicted)', flags=re.IGNORECASE)
RE_PROTEIN_NODE = re.compile(r'NODE_', flags=re.IGNORECASE)
RE_PROTEIN_POTENTIAL_CONTIG_NAME = re.compile(r'(genome|shotgun)', flags=re.IGNORECASE)
RE_PROTEIN_DOMAIN_CONTAINING = re.compile(r'domain-containing protein', flags=re.IGNORECASE)
RE_PROTEIN_NO_LETTERS = re.compile(r'[^A-Za-z]')
RE_PROTEIN_SUSPECT_CHARS = re.compile(r'[.@=?%]')
RE_PROTEIN_PERIOD_SEPARATOR = re.compile(r'([a-zA-Z0-9]+)\.([a-zA-Z0-9]+)')
RE_PROTEIN_WRONG_PRIMES = re.compile(r'[\u2032\u0060\u00B4]')  # prime (′), grave accent (`), acute accent (´)
RE_PROTEIN_WEIGHT = re.compile(r' [0-9]+(?:\.[0-9]+)? k?da ', flags=re.IGNORECASE)
RE_PROTEIN_SYMBOL = re.compile(r'[A-Z][a-z]{2}[A-Z][0-9]?')
RE_DOMAIN_OF_UNKNOWN_FUCTION = re.compile(r'(DUF\d{3,4})', flags=re.IGNORECASE)
RE_UNCHARACTERIZED_PROTEIN_FAMILY = re.compile(r'(UPF\d{3,4})', flags=re.IGNORECASE)
RE_GENE_CAPITALIZED = re.compile(r'^[A-Z].+', flags=re.DOTALL)
RE_GENE_SUSPECT_CHARS = re.compile(r'[\?]', flags=re.DOTALL)
RE_GENE_SYMBOL = re.compile(r'[a-z]{3}[A-Z][0-9]?')


def combine_annotation(feature: dict):
    pseudogene_psc = None
    pseudogene_pscc = None
    if(bc.PSEUDOGENE in feature):
        pseudogene_psc = feature[bc.PSEUDOGENE].get('psc', None)
        pseudogene_pscc = feature[bc.PSEUDOGENE].get('pscc', None)
    ups = feature.get('ups', None)
    ips = feature.get('ips', None)
    psc = feature.get('psc', None)
    pscc = feature.get('pscc', None)
    expert_hits = feature.get('expert', [])

    gene = None
    genes = set()
    product = None
    db_xrefs = set()
    if(pseudogene_pscc):
        pseudogene_pscc_genes = pseudogene_pscc.get('gene', None)
        if(pseudogene_pscc_genes):
            pseudogene_pscc_genes = pseudogene_pscc_genes.replace('/', ',').split(',')
            genes.update(pseudogene_pscc_genes)
            gene = pseudogene_pscc_genes[0]
        pscc_product = pseudogene_pscc.get('product', None)
        if(pscc_product):
            product = pscc_product
        for db_xref in pseudogene_pscc['db_xrefs']:
            db_xrefs.add(db_xref)
    if(pseudogene_psc):
        pseudogene_psc_genes = pseudogene_psc.get('gene', None)
        if(pseudogene_psc_genes):
            pseudogene_psc_genes = pseudogene_psc_genes.replace('/', ',').split(',')
            genes.update(pseudogene_psc_genes)
            gene = pseudogene_psc_genes[0]
        psc_product = pseudogene_psc.get('product', None)
        if(psc_product):
            product = psc_product
        for db_xref in pseudogene_psc['db_xrefs']:
            db_xrefs.add(db_xref)
    if(pscc):
        pscc_product = pscc.get('product', None)
        if(pscc_product):
            product = pscc_product
        for db_xref in pscc['db_xrefs']:
            db_xrefs.add(db_xref)
    if(psc and psc.get('valid', True)):
        psc_genes = psc.get('gene', None)
        if(psc_genes):
            psc_genes = psc_genes.replace('/', ',').split(',')
            genes.update(psc_genes)
            gene = psc_genes[0]
        psc_product = psc.get('product', None)
        if(psc_product):
            product = psc_product
        for db_xref in psc['db_xrefs']:
            db_xrefs.add(db_xref)
    if(ups):
        for db_xref in ups['db_xrefs']:
            db_xrefs.add(db_xref)
    if(ips):
        ips_genes = ips.get('gene', None)
        if(ips_genes):
            ips_genes = ips_genes.replace('/', ',').split(',')
            genes.update(ips_genes)
            gene = ips_genes[0]
        ips_product = ips.get('product', None)
        if(ips_product):
            product = ips_product
        for db_xref in ips['db_xrefs']:
            db_xrefs.add(db_xref)
    rank = 0
    for hit in expert_hits:
        db_xrefs.update(hit.get('db_xrefs', []))
        expert_rank = hit['rank']
        if(expert_rank > rank):
            expert_genes = hit.get('gene', None)
            if(expert_genes):
                expert_genes = expert_genes.replace('/', ',').split(',')
                genes.update(expert_genes)
                gene = expert_genes[0]
            product = hit.get('product', None)
            rank = expert_rank

    if(product):
        product = revise_cds_product(product)
        if(product):
            if(cfg.compliant):
                product = insdc.revise_product_insdc(product)
            feature['product'] = product
        
            protein_gene_symbol = extract_protein_gene_symbol(product)
            if(protein_gene_symbol):
                genes.add(protein_gene_symbol)
            revised_genes = revise_cds_gene_symbols(genes)
            revised_gene = None
            if gene is not None:
                revised_gene = revise_cds_gene_symbols([gene])  # special treatment for selected gene symbol
                revised_gene = revised_gene[0] if len(revised_gene) > 0 else None
            if(revised_gene is None  and  len(revised_genes) >= 1):  # select first from gene symbol list if no symbol was selected before
                revised_gene = revised_genes[0]

            feature['gene'] = revised_gene
            feature['genes'] = sorted(revised_genes)
        else:
            mark_as_hypothetical(feature)
    else:
        mark_as_hypothetical(feature)

    feature['db_xrefs'] = sorted(list(db_xrefs))


def detect_feature_overlaps(genome: dict):
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


def calc_sorf_annotation_score(sorf: dict) -> int:
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


def extract_protein_gene_symbol(product: str) -> str:
    gene_symbols = []
    for part in product.split(' '):  # try to extract valid gene symbols
        m = RE_GENE_SYMBOL.fullmatch(part)
        if(m):
            symbol = m[0]
            log.info('fix gene: extract symbol from protein name. symbol=%s', symbol)
            gene_symbols.append(symbol)
        else:
            m = RE_PROTEIN_SYMBOL.fullmatch(part)  # extract protein names
            if(m):
                symbol = m[0]
                symbol = symbol[0].lower() + symbol[1:]
                log.info('fix gene: extract symbol from protein name. symbol=%s', symbol)
                gene_symbols.append(symbol)
    if(len(gene_symbols) == 0):  # None found
        return None
    elif(len(gene_symbols) == 1):  # found 1
        return gene_symbols[0]
    else:  # found more than one, take the 2nd as the 1st often describes a broader gene family like "xyz family trancsriptional regulator ..."
        return gene_symbols[1]


def revise_cds_gene_symbols(raw_genes: Sequence[str]):
    revised_genes = set()
    for gene in raw_genes:
        old_gene = gene
        if(RE_GENE_SUSPECT_CHARS.search(gene)):  # check for suspect characters -> remove gene symbol
            log.info('fix gene: remove gene symbol containing suspect chars. old=%s', old_gene)
            continue

        old_gene = gene
        gene = gene.replace('gene', '')
        if(gene != old_gene):  # remove gene literal
            log.info('fix gene: remove gene literal. new=%s, old=%s', gene, old_gene)

        old_gene = gene
        if(gene[-1] == '-'):  # remove orphan hyphen
            gene = gene[:-1]
            log.info('fix gene: remove orphan hypen. new=%s, old=%s', gene, old_gene)
        
        old_gene = gene
        gene = RE_MULTIWHITESPACE.sub(' ', gene).strip()  # revise whitespaces
        if(gene != old_gene):
            log.info('fix gene: revise whitespaces. new=%s, old=%s', gene, old_gene)

        old_gene = gene
        if(RE_GENE_CAPITALIZED.fullmatch(gene)):
            gene = gene[0].lower() + gene[1:]
            log.info('fix gene: lowercase first char. new=%s, old=%s', gene, old_gene)

        if(len(gene) >= 3):
            if(len(gene) <= 12):
                revised_genes.add(gene)
            else:
                old_gene = gene
                gene = extract_protein_gene_symbol(gene)
                if(gene):
                    revised_genes.add(gene)
    return list(revised_genes)


def revise_cds_product(product: str):
    """Revise product name for INSDC compliant submissions"""
    
    old_product = product
    product = RE_PROTEIN_WEIGHT.sub(' ', product)  # remove protein weight in (k)Da
    if(product != old_product):
        log.info('fix product: remove protein weight in (k)Da. new=%s, old=%s', product, old_product)

    old_product = product
    product = re.sub(RE_PROTEIN_PERIOD_SEPARATOR, r'\1-\2', product)  # replace separator periods
    if(product != old_product):
        log.info('fix product: replace separator periods. new=%s, old=%s', product, old_product)

    old_product = product
    product = RE_PROTEIN_SUSPECT_CHARS.sub('', product)  # remove suspect characters
    if(product != old_product):
        log.info('fix product: replace invalid characters. new=%s, old=%s', product, old_product)

    old_product = product
    product = RE_PROTEIN_WRONG_PRIMES.sub('\u0027', product)  # replace wrong prime characters with single quote (U+0027) (') according to https://www.ncbi.nlm.nih.gov/genome/doc/internatprot_nomenguide/
    if(product != old_product):
        log.info('fix product: replace wrong prime characters. new=%s, old=%s', product, old_product)

    old_product = product
    product = product.replace('FOG:', '')  # remove FOG ids
    if(product != old_product):
        log.info('fix product: replace FOG ids. new=%s, old=%s', product, old_product)

    old_product = product
    dufs = []  # replace DUF-containing products
    for m in RE_DOMAIN_OF_UNKNOWN_FUCTION.finditer(product):
        dufs.append(m.group(1).upper())
    if(len(dufs) >= 1):
        product = f"{' '.join(dufs)} domain{'s' if len(dufs) > 1 else ''}-containing protein"
        if(product != old_product):
            log.info('fix product: revise DUF. new=%s, old=%s', product, old_product)
    
    old_product = product
    if('conserved' in product.lower()):  # replace conserved UPF proteins
        upfs = []
        for m in RE_UNCHARACTERIZED_PROTEIN_FAMILY.finditer(product):
            upfs.append(m.group(1).upper())
        if(len(upfs) >= 1):
            product = f"{' '.join(upfs)} protein"
            if(product != old_product):
                log.info('fix product: revise UPF. new=%s, old=%s', product, old_product)

    old_product = product
    product = RE_PROTEIN_HOMOLOG.sub('-like protein', product)  # replace Homologs
    if(product != old_product):
        if(product.count('protein') == 2):
            product = product.replace('protein', '', 1)  # remove former protein term if existing
        log.info('fix product: replace Homolog. new=%s, old=%s', product, old_product)

    old_product = product
    product = RE_MULTIWHITESPACE.sub(' ', product).strip()  # revise whitespaces
    if(product != old_product):
        log.info('fix product: revise whitespaces. new=%s, old=%s', product, old_product)

    old_product = product
    product = RE_PROTEIN_PUTATIVE.sub('putative', product)  # replace putative synonyms)
    if(product != old_product):
        log.info('fix product: replace putative synonyms. new=%s, old=%s', product, old_product)

    old_product = product
    if(RE_PROTEIN_DOMAIN_CONTAINING.search(product)):  # replace domain name underscores in domain names
        product = product.replace('_', '-')
        if(product != old_product):
            log.info('fix product: replace domain name underscores. new=%s, old=%s', product, old_product)

    old_product = product
    if(
        RE_PROTEIN_CONTIG.search(product) or  # protein containing 'contig'
        RE_PROTEIN_NODE.search(product) or  # potential contig name (SPAdes)
        RE_PROTEIN_POTENTIAL_CONTIG_NAME.search(product) or  # potential contig name (SPAdes)
        RE_PROTEIN_NO_LETTERS.fullmatch(product)  # no letters -> set to Hypothetical
        ):  # remove suspect products and mark as hypothetical
        product = None
        log.info('remove product: mark proteins with suspect products as hypothetical. old=%s', old_product)

    return product


def mark_as_hypothetical(feature: dict):
    log.info(
        'marked as hypothetical: contig=%s, start=%i, stop=%i, strand=%s',
        feature['contig'], feature['start'], feature['stop'], feature['strand']
    )
    feature['hypothetical'] = True
    feature['gene'] = None
    feature['genes'] = []
    feature['product'] = bc.HYPOTHETICAL_PROTEIN


def get_adjacent_genes(feature: dict, features: Sequence[dict], neighbors=3):
    for idx, feat in enumerate(features):
        if feat['locus'] == feature['locus']:
            upstream_genes = []
            if(idx >= 1):
                start = idx - neighbors
                if(start < 0 ):
                    start = 0
                upstream_genes = features[start:idx]
            downstream_genes = []
            if(idx + 1 < len(features)):
                end = idx + 1 + neighbors
                if(end > len(features)):
                    end = len(features)
                downstream_genes = features[idx+1:end]
            upstream_genes.extend(downstream_genes)
            for gene in upstream_genes:
                log.debug(
                    'extracted neighbor genes: contig=%s, start=%i, stop=%i, gene=%s, product=%s',
                    gene['contig'], gene['start'], gene['stop'], gene.get('gene', '-'), gene.get('product', '-')
                )
            return upstream_genes
    return []


def select_gene_symbols(features: Sequence[dict]):
    improved_genes = []
    for feat in [f for f in features if len(f.get('genes', [])) > 1]:  # all CDS/sORF with multiple gene symbols
        old_gene_symbol = feat['gene']
        gene_symbol_prefixes = set([symbol[:3] for symbol in feat['genes'] if len(symbol) > 3])
        if(len(gene_symbol_prefixes) == 1 ):  # multiple gene symbols of the same prefix: 
            product_parts = feat.get('product', '').split()
            for gene_symbol in feat['genes']:
                protein_symbol = gene_symbol[0].upper() + gene_symbol[1:]
                if(protein_symbol == product_parts[-1]):  # gene symbol is last part of product description which often is a specific gene/protein name
                    if(gene_symbol != old_gene_symbol):
                        feat['gene'] = gene_symbol
                        log.info(
                            'gene product symbol selection: contig=%s, start=%i, stop=%i, new-gene=%s, old-gene=%s, genes=%s, product=%s',
                            feat['contig'], feat['start'], feat['stop'], gene_symbol, old_gene_symbol, ','.join(feat['genes']), feat.get('product', '-')
                        )
                        improved_genes.append(feat)
        else:  # multiple gene symbols of varying prefixes are available, e.g. acrS, envR
            log.debug(
                'select gene symbol: contig=%s, start=%i, stop=%i, gene=%s, genes=%s, product=%s',
                feat['contig'], feat['start'], feat['stop'], feat.get('gene', '-'), ','.join(feat['genes']), feat.get('product', '-')
            )
            adjacent_genes = get_adjacent_genes(feat, features, neighbors=3)
            adjacent_gene_symbol_lists = [gene.get('genes', []) for gene in adjacent_genes]
            adjacent_gene_symbols = [item for sublist in adjacent_gene_symbol_lists for item in sublist]  # flatten lists
            adjacent_gene_symbol_prefixes = [gene_symbol[:3] for gene_symbol in adjacent_gene_symbols if len(gene_symbol) > 3]  # extract gene symbol prefixes, e.g. tra for traI, traX, traM
            adjacent_gene_symbol_prefix_counts = {}
            for gene_symbol_prefix in adjacent_gene_symbol_prefixes:
                if gene_symbol_prefix in adjacent_gene_symbol_prefix_counts:
                    adjacent_gene_symbol_prefix_counts[gene_symbol_prefix] += 1
                else:
                    adjacent_gene_symbol_prefix_counts[gene_symbol_prefix] = 1
            log.debug('neighbor gene symbol prefix counts: %s', adjacent_gene_symbol_prefix_counts)
            count = 0
            selected_gene_symbol = old_gene_symbol
            for gene_symbol in feat['genes']:
                gene_symbol_prefix = gene_symbol[:3]
                gene_symbol_count = adjacent_gene_symbol_prefix_counts.get(gene_symbol_prefix, 0)
                if gene_symbol_count > count:  # select symbol if its prefix is dominant in the gene neighborhood (neihboorhood of 3 genes up-/downstream as operon proxy)
                    selected_gene_symbol = gene_symbol
                    count = gene_symbol_count
            if(selected_gene_symbol != old_gene_symbol):
                feat['gene'] = selected_gene_symbol
                log.info(
                    'gene neighborhood symbol selection: contig=%s, start=%i, stop=%i, new-gene=%s, old-gene=%s, genes=%s, product=%s',
                    feat['contig'], feat['start'], feat['stop'], selected_gene_symbol, old_gene_symbol, ','.join(feat['genes']), feat.get('product', '-')
                )
                improved_genes.append(feat)
    return improved_genes