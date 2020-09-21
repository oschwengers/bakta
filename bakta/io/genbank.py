
import logging

from Bio import SeqIO


log = logging.getLogger('io:genbank')


############################################################################
# General terms
# 
# LOCUS       NC_002695            5498578 bp    DNA     circular CON 15-AUG-2018
# DEFINITION  Escherichia coli O157:H7 str. Sakai DNA, complete genome.
# SOURCE      Escherichia coli O157:H7 str. Sakai
# ORGANISM    Escherichia coli O157:H7 str. Sakai
#             Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacterales;
#             Enterobacteriaceae; Escherichia.
############################################################################

############################################################################
# CDS features
# 
# source          1..5498578
#                 /organism="Escherichia coli O157:H7 str. Sakai"
#                 /mol_type="genomic DNA"
#                 /strain="Sakai"
#                 /sub_strain="RIMD 0509952"
#                 /serovar="O157:H7"
#                 /db_xref="taxon:386585"
############################################################################

############################################################################
# CDS features
# 
# gene            232497..233300
#                 /gene="dkgB"
#                 /locus_tag="ECs_0203"
# CDS             232497..233300
#                 /gene="dkgB"
#                 /locus_tag="ECs_0203"
#                 /codon_start=1
#                 /transl_table=11
#                 /product="2,5-diketo-D-gluconate reductase B"
#                 /protein_id="NP_308230.1"
#                 /translation="MAIPAFG.....PEWD"
############################################################################

############################################################################
# tRNA features
# 
# tRNA            232258..232334
#                 /gene="aspU"
#                 /locus_tag="ECs_R0026"
#                 /product="tRNA-Asp"
#                 /note="anticodon: gtc"
############################################################################

############################################################################
# rRNA features
# 
# gene            232090..232200
#                 /gene="rrfH"
#                 /locus_tag="ECs_R0003"
# rRNA            232090..232200
#                 /gene="rrfH"
#                 /locus_tag="ECs_R0003"
#                 /product="5S ribosomal RNA"
############################################################################

############################################################################
# Feature inference terms
#
# for further information: https://www.ncbi.nlm.nih.gov/genbank/evidence/
# 
# tRNA: 'profile:tRNAscan:2.0'
# tmRNA: 'profile:aragorn:1.2'
# rRNA: 'profile:Rfam:%s' % subject_id
# ncRNA genes: 'profile:Rfam:%s' % subject_id
# ncRNA regions: 'profile:Rfam:%s' % subject_id
# CDS hyp: 'ab initio prediction:Prodigal:2.6'
# CDS PSC: 'similar to AA sequence:UniProtKB:%s' % psc[DB_PSC_COL_UNIREF90]
# CDS IPS: 'similar to AA sequence:UniProtKB:%s' % psc[DB_PSC_COL_UNIREF100]
# CDS sORF: 'similar to AA sequence:UniProtKB:%s' % psc[DB_PSC_COL_UNIREF100]
############################################################################


def write_genbank(annotations, genbank_path):
    """Export annotations in GenBank format."""

    return