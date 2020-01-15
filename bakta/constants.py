
############################################################################
# ORF prediction setttings
############################################################################
MIN_ORF_LENGTH = 7  # smallest non-fragmentated protein found in UniProtKB/Swiss-Prot
MAX_ORF_LENGTH = 30  # smallest ORF length detected by Prodigal


############################################################################
# Protein identification settings
############################################################################
MIN_PROTEIN_IDENTITY = 0.9  # min protein identity for PSC detection
MIN_PROTEIN_COVERAGE = 0.8  # min protein coverage for PSC detection
HYPOTHETICAL_PROTEIN = 'hypothetical protein'  # hypothetical protein product description


############################################################################
# DB identifiers
############################################################################
DB_XREF_UNIREF_100 = 'UniRef100'
DB_XREF_UNIREF_90 = 'UniRef90'
DB_XREF_UNIREF_50 = 'UniRef50'
DB_XREF_UNIPARC = 'UniParc'
DB_XREF_UNIPROTKB = 'UniProtKB'
DB_XREF_REFSEQ_NRP = 'RefSeq_NRP'
DB_XREF_EC = 'EC'
DB_XREF_COG = 'COG'
DB_XREF_PFAM = 'Pfam'
DB_XREF_IS = 'IS'
DB_XREF_GO = 'GO'


############################################################################
# DB identifier prefixes
############################################################################
DB_PREFIX_UNIREF_100 = 'UniRef100_'
DB_PREFIX_UNIREF_90 = 'UniRef90_'
DB_PREFIX_UNIREF_50 = 'UniRef50_'
DB_PREFIX_UNIPARC = 'UPI'
DB_PREFIX_REFSEQ_NRP = 'WP_'
DB_PREFIX_COG = 'COG'

############################################################################
# Inference tags
############################################################################
INFERENCE_T_RNA = 'COORDINATES:profile:tRNAscan-SE'


############################################################################
# INSDC feature tags
############################################################################
INSDC_FEATURE_T_RNA = 'tRNA'
INSDC_FEATURE_TM_RNA = 'tmRNA'
INSDC_FEATURE_R_RNA = 'rRNA'
# rRNA size via /product
INSDC_FEATURE_NC_RNA = 'ncRNA'
# mandatory: /ncRNA_class="TYPE"
INSDC_FEATURE_REPEAT_REGION = 'repeat_region'
# /rpt_family="text"
# /rpt_type=<repeat_type>
# /rpt_unit_range=<base_range>
# /rpt_unit_seq="text"
INSDC_FEATURE_CDS = 'CDS'

############################################################################
# Miscellaneous feature tags
############################################################################
CITATION = 'Schwengers O., Goesmann A. (2020)\nBakta: comprehensive and rapid annotation of bacterial genomes.\nGitHub https://github.com/oschwengers/bakta'
