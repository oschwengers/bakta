[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/oschwengers/bacanno/blob/master/LICENSE)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/bakta.svg)
![GitHub release](https://img.shields.io/github/release/oschwengers/bakta.svg)
![PyPI](https://img.shields.io/pypi/v/cb-bakta.svg)
![PyPI - Status](https://img.shields.io/pypi/status/bakta.svg)
![Conda](https://img.shields.io/conda/v/bioconda/bakta.svg)
![Conda](https://img.shields.io/conda/pn/bioconda/bakta.svg)

# Bakta: comprehensive annotation of bacterial genomes

## Contents

- [Description](#description)
- [Input/Output](#inputoutput)
- [Examples](#examples)
- [Installation](#installation)
  - [Bioconda](#bioconda)
  - [GitHub](#github)
  - [Pip](#pip)
  - [Dependencies](#dependencies)
- [Annotation workflow](#annoation_workflow)
- [Database](#database)
- [Usage](#usage)
- [Citation](#citation)
- [Issues & Feature Requests](#issues)

## Description

**TL;DR**
Bakta provides a comprehensive, standardized and *Dbxref*-rich annotation tailored to bacterial genomes
without the need & option to fine-tune parameters fitting the space between
large but computationally-demanding annotation pipelines like PGAP and fast & highly-customizable tools like Prokka.

Lorem ipsum ...

## Input/Output

### Input

Bakta accepts bacterial (complete / draft) assemblies in fasta format.

### Output

Bakta provides detailed information on each annotated feature in a standardized machine-readable JSON format.
In addition, Bakta supports the following standard file formats:

- `tsv`: dense human-readable information as tab separated values
- `GFF3`: standard GFF3
- `GenBank`: GenBank file format produced via BioPython

## Examples

Simple:

```bash
$ bakta -db ~/db genome.fasta
```

Expert: writing results to `results` directory with verbose output using 8 threads:

```bash
$ bakta -db ~/db --output results/ --verbose --threads 8 genome.fasta
```

## Installation

Bakta can be installed/used in 3 different ways.

In all cases, the custom database must be downloaded which we provide for download:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXX.svg)](https://doi.org/10.5281/zenodo.XXX)

To automatically install all required dependencies we highly encourage to use BioConda!

### BioConda

1. install Bakta via Conda
2. download & extract the database

Example:

```bash
$ conda install -c bioconda bakta
$ wget <XYZ>/db.tar.gz
$ tar -xzf db.tar.gz
$ rm db.tar.gz
$ bakta --db ./db genome.fasta
```

### GitHub

1. clone the the repository
2. download & extract the database
3. install 3rd party binaries (-> Dependencies)

Example:

```bash
$ git clone git@github.com:oschwengers/bakta.git
$ wget <XYZ>/db.tar.gz
$ tar -xzf db.tar.gz
$ rm db.tar.gz
$ bakta/bin/bakta --db ./db genome.fasta
```

### Pip

1. install Bakta per pip
2. download and extract the database
3. install 3rd party binaries (-> Dependencies)

Bakta/database (1./2.):

```bash
$ python3 -m pip install bakta
$ wget <XYZ>/db.tar.gz
$ tar -xzf db.tar.gz
$ rm db.tar.gz
$ bakta --db ./db genome.fasta
```

### Dependencies

Bakta depends on the following 3rd party executables which must be installed and available:

- tRNAscan-SE (2.0.6) <https://doi.org/10.1101/614032> <http://lowelab.ucsc.edu/tRNAscan-SE>
- Aragorn (1.2.38) <http://dx.doi.org/10.1093/nar/gkh152> <http://130.235.244.92/ARAGORN>
- INFERNAL (1.1.3) <https://dx.doi.org/10.1093%2Fbioinformatics%2Fbtt509> <http://eddylab.org/infernal>
- Prodigal (2.6.3) <https://dx.doi.org/10.1186%2F1471-2105-11-119> <https://github.com/hyattpd/Prodigal>
- Diamond (2.0.2) <https://doi.org/10.1038/nmeth.3176> <https://github.com/bbuchfink/diamond>

Ubuntu:

```bash
$ sudo apt install aragorn infernal prodigal diamond-aligner
```

tRNAscan-se must be installed manually as v2.0 is currently not yet available via standard Ubuntu packages.

## Annotation workflow

### RNAs

1. tRNA genes: tRNAscan-SE 2.0
2. tmRNA genes: Aragorn
3. rRNA genes: Infernal vs. Rfam rRNA covariance models
4. ncRNA genes: Infernal vs. Rfam ncRNA covariance models
5. ncRNA regions: Infernal vs. Rfam ncRNA covariance models

Bakta distinguishes ncRNA genes and regions in order to enable the distinct handling thereof during the annotation process.
ncRNA gene types:

- sRNA
- antisense
- ribozyme
- antitoxin

ncRNA region types:

- riboswitch
- thermoregulator
- leader
- frameshift element

### Coding sequences

The structural prediction is conducted via Prodigal and complemented by a
custom detection of short open reading freames (**sORF**) smaller than 30 aa.

**CDS**:

1. Prediction via Prodigal
2. Detection of unique protein sequences (**UPS**)s via **MD5** hashes and lookup of related (**IPS**)s.
3. Alignment of remainder via Diamond vs. UniProt's UniRef90 based protein sequence clusters (**PSC**)s.
4. Combine available **IPS** & **PSC** annotations erasing redundancy and obeying specificity hierarchies.

**sORFs**:

1. Custom detection & extraction of **sORF**s with amino acid lengths < 30
2. Filter via strict feature type-dependent overlap filters with annotated features
3. Detection of **UPS**s via **MD5** hashes and lookup of related **IPS**s.

Only **sORF** which are detected by their identity (100% coverage & 100% sequence identity)
will be included in the annotation. A more sensitive but also more false-positive prone approach is currently tested.

## Database

The Bakta database is built on **IPS**s and **PSC**s from:

- **IPS**: UniProt UniRef100
- **PSC**: UniProt UniRef90

which have been comprehensively pre-annotated integrating
annotations & database cross references (*Dbxrefs*) from:

- NCBI nonredundant proteins ('WP_*', exact matches)
- NCBI COG db (80% coverage & 90% identity)
- GO terms (via SwissProt)
- EC (via SwissProt)
- NCBI AMRFinderPlus (**IPS** exact matches, **PSC** HMM hits achieving the trusted cutoff)
- ISFinder db (90% coverage & 99% identity)

For the sake of performance and data size, all pre-annotation information is stored in a compact
SQLite database which can be downloaded here:
[![DOI](https://zenodo.org/badge/DOI/<DOI>.svg)](https://doi.org/<DOI>)

- [<DB_URL>](<DB_URL>)

DB size:

- zipped: 23 Gb
- unzipped 43 Gb

## Usage

Usage:

```bash
bakta --help
usage: bakta [--db DB] [--min-contig-length MIN_CONTIG_LENGTH] [--prefix PREFIX] [--output OUTPUT] [--tsv] [--gff3] [--genbank] [--embl]
             [--genus GENUS] [--species SPECIES] [--strain STRAIN] [--plasmid PLASMID]
             [--prodigal-tf PRODIGAL_TF] [--keep-contig-names] [--locus LOCUS] [--locus-tag LOCUS_TAG] [--gram {+,-,?}] [--complete]
             [--skip-trna] [--skip-tmrna] [--skip-rrna] [--skip-ncrna] [--skip-ncrna-region] [--skip-cds] [--skip-sorf] [--skip-gap]
             [--help] [--verbose] [--threads THREADS] [--tmp-dir TMP_DIR] [--version] [--citation]
             <genome>

Comprehensive and rapid annotation of bacterial genomes.

positional arguments:
  <genome>              (Draft) genome in fasta format

Input / Output:
  --db DB, -d DB        Database path (default = <bakta_path>/db)
  --min-contig-length MIN_CONTIG_LENGTH, -m MIN_CONTIG_LENGTH
                        Minimum contig size (default = 1)
  --prefix PREFIX, -p PREFIX
                        Prefix for output files
  --output OUTPUT, -o OUTPUT
                        Output directory (default = current working directory)
  --tsv                 Write TSV annotation file
  --gff3                Write GFF3 annotation file
  --genbank             Write GenBank annotation file
  --embl                Write EMBL annotation file

Organism:
  --genus GENUS         Genus name
  --species SPECIES     Species name
  --strain STRAIN       Strain name
  --plasmid PLASMID     Plasmid name

Annotation:
  --prodigal-tf PRODIGAL_TF
                        Path to existing Prodigal training file to use for CDS prediction
  --translation-table TRANSLATION_TABLE
                        Translation table to use (default = 11)
  --keep-contig-names   Keep original contig names
  --locus LOCUS         Locus prefix
  --locus-tag LOCUS_TAG
                        Locus tag prefix
  --gram {+,-,?}        Gram type: +/-/? (default = '?')
  --complete            Replicons (chromosome/plasmid[s]) are complete

Workflow:
  --skip-trna           Skip tRNA detection & annotation
  --skip-tmrna          Skip tmRNA detection & annotation
  --skip-rrna           Skip rRNA detection & annotation
  --skip-ncrna          Skip ncRNA detection & annotation
  --skip-ncrna-region   Skip ncRNA region detection & annotation
  --skip-cds            Skip CDS detection & annotation
  --skip-sorf           Skip sORF detection & annotation
  --skip-gap            Skip gap detection & annotation

General:
  --help, -h            Show this help message and exit
  --verbose, -v         Print verbose information
  --threads THREADS, -t THREADS
                        Number of threads to use (default = number of available CPUs)
  --tmp-dir TMP_DIR     Location for temporary files (default = system dependent auto detection)
  --version             Show version number and exit
  --citation            Print citation
```

## Citation

A manuscript is in preparation... stay tuned!
To temporarily cite our work, please transitionally refer to:
> Schwengers O., Goesmann A. (2020) Bakta: comprehensive annotation of bacterial genomes. GitHub https://github.com/oschwengers/bakta

As Bakta takes advantage of ISFinder database, please also cite:
> <ISFinder_citation>

## Issues and Feature Requests

If you run into any issues with Bakta, we'd be happy to hear about it!
Please, execute bakta in verbose mode (`-v`) and do not hesitate
to file an issue including as much of the following as possible:

- a detailed description of the issue
- command line output
- log file (`<prefix>.log`)
- result file (`<prefix>.json`) _if possible_
- a reproducible example of the issue with an input file that you can share _if possible_
