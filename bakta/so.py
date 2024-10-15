from collections import namedtuple


SO = namedtuple('SO', ['name', 'id'])


SO_REPLICON = SO('replicon', 'SO:0001235')
SO_REPLICON_CHROMOSOME = SO('chromosome', 'SO:0000340')
SO_REPLICON_PLASMID = SO('plasmid', 'SO:0000155')
SO_CONTIG = SO('sequence', 'SO:0000149')

SO_OPERON = SO('operon', 'SO:0000178')
SO_PROMOTER = SO('promoter', 'SO:0000167')
SO_OPERATOR = SO('operator', 'SO:0000057')

SO_REGULATORY_REGION = SO('regulatory_region', 'SO:0005836')
SO_CIS_REG = SO('transcriptional_cis_regulatory_region', 'SO:0001055')
SO_CIS_REG_ATTENUATOR = SO('attenuator', 'SO:0000140')
SO_CIS_REG_RESPONSE_ELEMENT = SO('response_element', 'SO:0002205')
SO_CIS_REG_RIBOSWITCH = SO('riboswitch', 'SO:0000035')
SO_CIS_REG_THERMOMETER = SO('RNA_thermometer', 'SO:0002168')
SO_CIS_REG_FRAMESHIFT = SO('cis_regulatory_frameshift_element', 'SO:0001427')
SO_CIS_REG_RECODING_STIMULATION_REGION = SO('recoding_stimulatory_region', 'SO:1001268')
SO_CIS_REG_RIBOSOME_BINDING_SITE = SO('ribosome_entry_site', 'SO:0000139')

SO_TERMINATOR = SO('terminator', 'SO:0000141')
SO_TERMINATOR_RHO_DEPENDENT = SO('rho_dependent_bacterial_terminator', 'SO:0000981')
SO_TERMINATOR_RHO_INDEPENDENT = SO('rho_independent_bacterial_terminator', 'SO:0000982')

SO_GENE = SO('gene', 'SO:0000704')
SO_PSEUDO_GENE = SO('pseudogene', 'SO:0000336')
SO_CDS = SO('CDS', 'SO:0000316')
SO_SORF = SO('sORF', 'SO:0002028')
SO_SIGNAL_PEPTIDE = SO('signal_peptide', 'SO:0000418')

SO_TRNA = SO('tRNA', 'SO:0000253')
SO_TRNA_ALA = SO('alanyl_tRNA', 'SO:0000254')
SO_TRNA_GLN = SO('glycyl_tRNA', 'SO:0000261')
SO_TRNA_GLU = SO('glutaminyl_tRNA', 'SO:0000259')
SO_TRNA_GLY = SO('glutamyl_tRNA', 'SO:0000260')
SO_TRNA_PRO = SO('prolyl_tRNA', 'SO:0000268')
SO_TRNA_MET = SO('methionyl_tRNA', 'SO:0000266')
SO_TRNA_ASP = SO('asparaginyl_tRNA', 'SO:0000256')
SO_TRNA_THR = SO('threonyl_tRNA', 'SO:0000270')
SO_TRNA_VAL = SO('valyl_tRNA', 'SO:0000273')
SO_TRNA_TYR = SO('tyrosyl_tRNA', 'SO:0000272')
SO_TRNA_CYS = SO('cysteinyl_tRNA', 'SO:0000258')
SO_TRNA_ILE = SO('isoleucyl_tRNA', 'SO:0000263')
SO_TRNA_SER = SO('seryl_tRNA', 'SO:0000269')
SO_TRNA_LEU = SO('leucyl_tRNA', 'SO:0000264')
SO_TRNA_TRP = SO('tryptophanyl_tRNA', 'SO:0000271')
SO_TRNA_LYS = SO('lysyl_tRNA', 'SO:0000265')
SO_TRNA_ASN = SO('aspartyl_tRNA', 'SO:0000257')
SO_TRNA_ARG = SO('arginyl_tRNA', 'SO:0001036')
SO_TRNA_HIS = SO('histidyl_tRNA', 'SO:0000262')
SO_TRNA_PHE = SO('phenylalanyl_tRNA', 'SO:0000267')
SO_TRNA_SELCYS = SO('selenocysteinyl_tRNA', 'SO:0005857')

SO_TMRNA = SO('tmRNA', 'SO:0000584')
SO_RRNA = SO('rRNA', 'SO:0000252')
SO_RRNA_5S = SO('rRNA_5S', 'SO:0000652')
SO_RRNA_16S = SO('rRNA_16S', 'SO:0001000')
SO_RRNA_23S = SO('rRNA_23S', 'SO:0001001')

SO_NCRNA_GENE = SO('ncRNA', 'SO:0000655')
SO_NCRNA_GENE_RIBOZYME = SO('ribozyme', 'SO:0000374')
SO_NCRNA_GENE_RNASEP = SO('RNase_P_RNA', 'SO:0000386')
SO_NCRNA_GENE_ANTISENSE = SO('antisense_RNA', 'SO:0000644')

SO_CRISPR = SO('CRISPR', 'SO:0001459')

SO_ORIC = SO('oriC', 'SO:0000953')
SO_ORIT = SO('oriT', 'SO:0000724')
SO_ORIV = SO('oriV', 'SO:0000952')

SO_PROPHAGE = SO('proviral_region', 'SO:0000113')
SO_IS = SO('insertion_sequence', 'SO:0000973')
SO_TRANSPOSABLE_ELEMENT = SO('transposable_element', 'SO:0000101')
SO_MOBILE_ELEMENT = SO('mobile_genetic_element', 'SO:0001037')

SO_GAP = SO('gap', 'SO:0000730')
SO_REPEAT = SO('repeat_region', 'SO:0000657')
