#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [bakta]
id: bakta
label: "Bakta: rapid & standardized annotation of bacterial genomes, MAGs & plasmids"

doc: |
      The software and documentation can be found here:
      https://github.com/oschwengers/bakta

      Necessary database files can be found here:
      https://doi.org/10.5281/zenodo.4247252

hints:
  SoftwareRequirement:
    packages:
      bakta:
      version: [ "1.6.1" ]

requirements:
  ResourceRequirement:
    ramMin: 4096
    coresMin: 1

inputs:
  - doc: Genome assembly in Fasta format
    id: fasta_file
    inputBinding: {position: 0}
    type: File
    format: edam:format_1929
  - doc: Database path (default = <bakta_path>/db)
    id: db
    inputBinding: {prefix: --db}
    type: Directory
  - doc: Output directory (default = current working directory)
    id: output
    inputBinding: {prefix: --output}
    type: Directory
  - doc: Output files prefix
    id: prefix
    inputBinding: {prefix: --prefix}
    type: string
  - doc: Min contig length
    id: min_contig_length
    inputBinding: {prefix: --min-contig-length}
    type: int
  - doc: Genus
    id: genus
    inputBinding: {prefix: --genus}
    type: string
  - doc: Species
    id: species
    inputBinding: {prefix: --species}
    type: string
  - doc: Strain
    id: strain
    inputBinding: {prefix: --strain}
    type: string
  - doc: All sequences are complete replicons (chromosome/plasmid[s])
    id: complete
    inputBinding: {prefix: --complete}
    type: boolean
  - doc: Prodigal training file for CDS prediction
    id: prodigal_tf_file
    inputBinding: {prefix: --prodigal-tf}
    type: File
  - doc: Translation table: 11/4 (default = 11)
    id: translation_table
    inputBinding: {prefix: --translation-table}
    type: int
  - doc: Locus tag
    id: locus_tag
    inputBinding: {prefix: --locus-tag}
    type: string
  - doc: Gram type (default = ?)
    id: gram
    inputBinding: {prefix: --gram}
    type: string
  - doc: Keep original contig headers
    id: keep_contig_headers
    inputBinding: {prefix: --keep-contig-headers}
    type: boolean
  - doc: Replicon information table (tsv/csv)
    id: replicons
    inputBinding: {prefix: --replicons}
    type: File
  - doc: Force Genbank/ENA/DDJB compliance
    id: compliant
    inputBinding: {prefix: --compliant}
    type: boolean
  - doc: Fasta file of trusted protein sequences for CDS annotation
    id: proteins
    inputBinding: {prefix: --proteins}
    type: File
  - doc: Skip tRNA detection & annotation
    id: skip_tRNA
    inputBinding: {prefix: --skip-trna}
    type: boolean
  - doc: Skip tmRNA detection & annotation
    id: skip_tmrna
    inputBinding: {prefix: --skip-tmrna}
    type: boolean
  - doc: Skip rRNA detection & annotation
    id: skip_rrna
    inputBinding: {prefix: --skip-rrna}
    type: boolean
  - doc: Skip ncRNA detection & annotation
    id: skip_ncrna
    inputBinding: {prefix: --skip-ncrna}
    type: boolean
  - doc: Skip ncRNA region detection & annotation
    id: skip_ncrna_region
    inputBinding: {prefix: --skip-ncrna-region}
    type: boolean
  - doc: Skip CRISPR detection & annotation
    id: skip_crispr
    inputBinding: {prefix: --skip-crispr}
    type: boolean
  - doc: Skip CDS detection & annotation
    id: skip_cds
    inputBinding: {prefix: --skip-cds}
    type: boolean
  - doc: Skip Pseudogene detection & annotation
    id: skip_pseudo
    inputBinding: {prefix: --skip-pseudo}
    type: boolean
  - doc: Skip sORF detection & annotation
    id: skip_sorf
    inputBinding: {prefix: --skip-sorf}
    type: boolean
  - doc: Skip gap detection & annotation
    id: skip_gap
    inputBinding: {prefix: --skip-gap}
    type: boolean
  - doc: Skip ori detection & annotation
    id: skip_ori
    inputBinding: {prefix: --skip-ori}
    type: boolean
  - doc: Skip genome plotting
    id: skip_plot
    inputBinding: {prefix: --skip-plot}
    type: boolean
  - doc: Directory for temporary files (default = system dependent auto detection)
    id: tmp_dir
    inputBinding: {prefix: --tmp-dir}
    type: Directory
  - doc: Threads
    id: threads
    inputBinding: {prefix: --threads}
    type: int

outputs:
  - doc: Hypothetical CDS AA sequences as Fasta
    id: hypo_sequences_cds
    type: File
    format: edam:format_2200
    outputBinding: {glob: '*.hypotheticals.faa'}
  - doc: Information on hypothetical CDS as TSV
    id: hypo_annotation_tsv
    type: File
    format: edam:format_3475
    outputBinding: {glob: '*.hypotheticals.tsv'}
  - doc: Annotation as TSV
    id: annotation_tsv
    type: File
    format: edam:format_3475
    outputBinding: {glob: '*.tsv'}
  - doc: Annotation summary as txt
    id: summary_txt
    type: File
    format: edam:format_2330
    outputBinding: {glob: '*.txt'}
  - doc: Annotation as JSON
    id: annotation_json
    type: File
    format: edam:format_3464
    outputBinding: {glob: '*.json'}
  - doc: Annotation as GFF3
    id: annotation_gff3
    type: File
    format: edam:format_1939
    outputBinding: {glob: '*.gff3'}
  - doc: Annotation as GenBank
    id: annotation_gbff
    type: File
    format: edam:format_1936
    outputBinding: {glob: '*.gbff'}
  - doc: Annotation as EMBL
    id: annotation_embl
    type: File
    format: edam:format_1927
    outputBinding: {glob: '*.embl'}
  - doc: Genome Sequences as Fasta
    id: sequences_fna
    type: File
    format: edam:format_2200
    outputBinding: {glob: '*.fna'}
  - doc: Gene DNA sequences as Fasta
    id: sequences_fna
    type: File
    format: edam:format_2200
    outputBinding: {glob: '*.ffn'}
  - doc: CDS AA sequences as Fasta
    id: sequences_cds
    type: File
    format: edam:format_2200
    outputBinding: {glob: '*.faa'}
  - doc: Circular genome plot as PNG
    id: plot_png
    type: File
    format: edam:format_3603
    outputBinding: {glob: '*.png'}
  - doc: Circular genome plot as SVG
    id: plot_svg
    type: File
    format: edam:format_3604
    outputBinding: {glob: '*.svg'}

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0003-4216-2721
    s:email: mailto:oliver.schwengers@computational.bio.uni-giessen.de
    s:name: Oliver Schwengers

s:citation: https://doi.org/10.1099/mgen.0.000685
s:codeRepository: https://github.com/oschwengers/bakta
s:license: https://spdx.org/licenses/GNU GPL3
s:programmingLanguage: Python

$namespaces:
  s: https://schema.org/
  edam: http://edamontology.org/

$schemas:
 - https://schema.org/version/latest/schema.rdf
 - http://edamontology.org/EDAM_1.18.owl
