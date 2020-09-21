
import bakta.constants as bc

def combine_ups_psc_annotation(feature):
    ips = feature.get('ips', None)
    psc = feature.get('psc', None)

    gene = ''
    product = bc.HYPOTHETICAL_PROTEIN
    if(psc):
        psc_gene = psc.get('gene', None)
        if(psc_gene):
            gene = psc_gene
        psc_product = psc.get('product', None)
        if(psc_product):
            product = psc_product
    if(ips):
        ips_gene = ips.get('gene', None)
        if(ips_gene):
            gene = ips_gene
        ips_product = ips.get('product', None)
        if(ips_product):
            product = ips_product
    feature['gene'] = gene
    feature['product'] = product