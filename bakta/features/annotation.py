
import bakta.constants as bc

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
    feature['db_xrefs'] = list(db_xrefs)
