
import bakta.constants as bc

def combine_ups_psc_annotation(feature):
    ups = feature.get('ups', None)
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
    if(ups):
        ups_gene = ups.get('gene', None)
        if(ups_gene):
            gene = ups_gene
        ups_product = ups.get('product', None)
        if(ups_product):
            product = ups_product
    feature['gene'] = gene
    feature['product'] = product