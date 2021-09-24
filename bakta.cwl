#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [bakta]
id: bakta
label: "Rapid & standardized annotation of bacterial genomes & plasmids."

doc: |
      The software and documentation can be found here:
      https://github.com/oschwengers/bakta

      Necessary database files can be found here:
      https://doi.org/10.5281/zenodo.4247252

hints:
  SoftwareRequirement:
    packages:
      bakta:
      version: [ "1.1" ]

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
  - doc: Threads
    id: threads
    inputBinding: {prefix: --threads}
    type: int

outputs:
  - doc: Annotation as TSV
    id: annotation_tsv
    type: File
    format: edam:format_1929
    outputBinding: {glob: '*.tsv'}
  - doc: Annotation as JSON
    id: annotation_json
    type: File
    format: edam:format_1929
    outputBinding: {glob: '*.json'}
  - doc: Annotation as GFF3
    id: annotation_gff3
    type: File
    format: edam:format_1929
    outputBinding: {glob: '*.gff3'}
  - doc: Annotation as GenBank
    id: annotation_gbff
    type: File
    format: edam:format_1929
    outputBinding: {glob: '*.gbff'}
  - doc: Annotation as EMBL
    id: annotation_embl
    type: File
    format: edam:format_1929
    outputBinding: {glob: '*.embl'}
  - doc: Genome Sequences as Fasta
    id: sequences_fna
    type: File
    format: edam:format_1929
    outputBinding: {glob: '*.fna'}
  - doc: CDS as Fasta
    id: sequences_cds
    type: File
    format: edam:format_1929
    outputBinding: {glob: '*.faa'}

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0003-4216-2721
    s:email: mailto:oliver.schwengers@computational.bio.uni-giessen.de
    s:name: Oliver Schwengers

# s:citation:
s:codeRepository: https://github.com/oschwengers/bakta
s:license: https://spdx.org/licenses/GNU GPL3
s:programmingLanguage: Python

$namespaces:
  s: https://schema.org/
  edam: http://edamontology.org/

$schemas:
 - https://schema.org/version/latest/schema.rdf
 - http://edamontology.org/EDAM_1.18.owl
