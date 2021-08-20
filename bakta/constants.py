import re


############################################################################
# Common regex
############################################################################
RE_MULTIWHITESPACE = re.compile(r'\s+')


############################################################################
# sORF prediction setttings
############################################################################
MIN_SORF_LENGTH = 7  # smallest non-fragmentated protein found in UniProtKB/Swiss-Prot
MAX_SORF_LENGTH = 30  # smallest ORF length detected by Prodigal


############################################################################
# Protein identification settings
############################################################################
MIN_PSCC_IDENTITY = 0.5  # min protein identity for PSC detection
MIN_PSC_COVERAGE = 0.8  # min protein coverage for PSC detection
MIN_PSC_IDENTITY = 0.9  # min protein identity for PSC detection
MIN_SORF_COVERAGE = 0.9  # min sORF coverage for PSC detection
MIN_SORF_IDENTITY = 0.9  # min sORF identity for PSC detection
HYPOTHETICAL_PROTEIN = 'hypothetical protein'  # hypothetical protein product description


############################################################################
# DB identifiers
############################################################################
DB_XREF_UNIPROTKB = 'UniProtKB/TrEMBL'
DB_XREF_REFSEQ_NRP = 'RefSeq'
DB_XREF_EC = 'EC'
DB_XREF_COG = 'COG'
DB_XREF_PFAM = 'Pfam'
DB_XREF_RFAM = 'RFAM'
DB_XREF_IS = 'IS'
DB_XREF_GO = 'GO'


############################################################################
# DB identifier prefixes
############################################################################
DB_PREFIX_UNIREF_100 = 'UniRef100_'
DB_PREFIX_UNIREF_90 = 'UniRef90_'
DB_PREFIX_UNIREF_50 = 'UniRef50_'
DB_PREFIX_UNIPARC = 'UPI'  # https://registry.identifiers.org/registry/uniparc
DB_PREFIX_REFSEQ_NRP = 'WP_'
DB_PREFIX_GO = 'GO:'  # https://registry.identifiers.org/registry/go
DB_PREFIX_COG = 'COG'  # https://www.ebi.ac.uk/miriam/main/collections/MIR:00000296
DB_PREFIX_IS = 'IS'  # https://www.ebi.ac.uk/miriam/main/collections/MIR:00000173
DB_PREFIX_PFAM = 'PF'  # https://registry.identifiers.org/registry/pfam


############################################################################
# Bakta feature tags
############################################################################
FEATURE_T_RNA = 'tRNA'
FEATURE_TM_RNA = 'tmRNA'
FEATURE_R_RNA = 'rRNA'
FEATURE_NC_RNA = 'ncRNA'
FEATURE_NC_RNA_REGION = 'ncRNA-region'
FEATURE_CRISPR = 'crispr'
FEATURE_ORF = 'orf'
FEATURE_SORF = 'sorf'
FEATURE_CDS = 'cds'
FEATURE_GAP = 'gap'
FEATURE_ORIC = 'oriC'
FEATURE_ORIV = 'oriV'
FEATURE_ORIT = 'oriT'

############################################################################
# INSDC feature tags
############################################################################
INSDC_FEATURE_T_RNA = 'tRNA'
INSDC_FEATURE_TM_RNA = 'tmRNA'
INSDC_FEATURE_R_RNA = 'rRNA'
# rRNA size via /product
INSDC_FEATURE_NC_RNA = 'ncRNA'
INSDC_FEATURE_NC_RNA_CLASS = 'ncRNA_class'  # mandatory: /ncRNA_class="TYPE"
INSDC_FEATURE_NC_RNA_CLASS_RNASEP = 'RNase_P_RNA'
INSDC_FEATURE_NC_RNA_CLASS_RIBOZYME = 'ribozyme'
INSDC_FEATURE_NC_RNA_CLASS_ANTISENSE = 'antisense_RNA'
INSDC_FEATURE_NC_RNA_CLASS_OTHER = 'other'
INSDC_FEATURE_REPEAT_REGION = 'repeat_region'
INSDC_FEATURE_REPEAT_FAMILY = 'rpt_family'  # /rpt_family="text"  # CRISPR
INSDC_FEATURE_REPEAT_TYPE = 'rpt_type'  # /rpt_type=<repeat_type>  # 'direct'
INSDC_FEATURE_REPEAT_UNIT_RANGE = 'rpt_unit_range'  # /rpt_unit_range=<base_range>
INSDC_FEATURE_REPEAT_UNIT_SEQ = 'rpt_unit_seq'  # /rpt_unit_seq="text"
INSDC_FEATURE_CDS = 'CDS'
INSDC_FEATURE_GAP = 'gap'
INSDC_FEATURE_ASSEMBLY_GAP = 'assembly_gap'
INSDC_FEATURE_MISC_FEATURE = 'misc_feature'
INSDC_FEATURE_MISC_RNA = 'misc_RNA'
INSDC_FEATURE_ORIGIN_TRANSFER = 'oriT'
INSDC_FEATURE_ORIGIN_REPLICATION = 'rep_origin'
INSDC_FEATURE_REGULATORY = 'regulatory'
INSDC_FEATURE_REGULATORY_CLASS = 'regulatory_class'
INSDC_FEATURE_REGULATORY_CLASS_ATTENUATOR = 'attenuator'
INSDC_FEATURE_REGULATORY_CLASS_RIBOSWITCH = 'riboswitch'
INSDC_FEATURE_REGULATORY_CLASS_RESPONSE_ELEMENT = 'response_element'
INSDC_FEATURE_REGULATORY_CLASS_RECODING_STIMULATORY_REGION = 'recoding_stimulatory_region'
INSDC_FEATURE_REGULATORY_CLASS_RIBOSOME_BINDING_SITE = 'ribosome_binding_site'
INSDC_FEATURE_REGULATORY_CLASS_OTHER = 'other'

############################################################################
# Feature ends
############################################################################
FEATURE_END_5_PRIME = '5-prime'
FEATURE_END_3_PRIME = '3-prime'
FEATURE_END_UNKNOWN = '?'

############################################################################
# Strand types prefixes
############################################################################
STRAND_FORWARD = '+'
STRAND_REVERSE = '-'
STRAND_UNKNOWN = '?'
STRAND_NA = '.'

############################################################################
# Replicon types, length thresholds & topology
############################################################################
REPLICON_CHROMOSOME = 'chromosome'
REPLICON_PLASMID = 'plasmid'
REPLICON_CONTIG = 'contig'

REPLICON_LENGTH_THRESHOLD_PLASMID = 112_000  # Nasuia deltocephalinicola -> DOI: 10.1093/gbe/evt118
REPLICON_LENGTH_THRESHOLD_CHROMOSOME = 2_800_000  # max plasmid length (except 1 outlier-> https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/

TOPOLOGY_CIRCULAR = 'circular'
TOPOLOGY_LINEAR = 'linear'

############################################################################
# Miscellaneous constants
############################################################################
DISCARD_TYPE_SPURIOUS = 'spurious'
DISCARD_TYPE_OVERLAP = 'overlap'

DB_VERSIONS_URL = 'https://raw.githubusercontent.com/oschwengers/bakta/master/db-versions.json'
CITATION = 'Schwengers et al. (2021)\nBakta: comprehensive and rapid annotation of bacterial genomes.\nGitHub https://github.com/oschwengers/bakta'
