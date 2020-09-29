[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/oschwengers/bacanno/blob/master/LICENSE)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/bakta.svg)
![GitHub release](https://img.shields.io/github/release/oschwengers/bakta.svg)
![PyPI](https://img.shields.io/pypi/v/cb-bakta.svg)
![PyPI - Status](https://img.shields.io/pypi/status/bakta.svg)
![Conda](https://img.shields.io/conda/v/bioconda/bakta.svg)
![Conda](https://img.shields.io/conda/pn/bioconda/bakta.svg)

# Bakta: comprehensive annotation of bacterial genomes.

## Contents
- [Description](#description)
- [Input/Output](#inputoutput)
- [Annotation workflow](#annoation_workflow)
- [Installation](#installation)
  - [Bioconda](#bioconda)
  - [GitHub](#github)
  - [Pip](#pip)
- [Usage](#usage)
- [Examples](#examples)
- [Database](#database)
- [Dependencies](#dependencies)
- [Citation](#citation)

## Description
**TL;DR**
Bakta lorem ipsum ...

Bakta lorem ipsum ...

## Input/Output

### Input
Bakta accepts bacterial (draft) assemblies in fasta format.

### Output
Lorem ipsum

In addition, Bakta writes the following files into the output directory:
-   `<prefix>`.plasmid.fasta: contigs classified as plasmids or plasmodal origin
-   `<prefix>`.chromosome.fasta: contigs classified as chromosomal origin
-   `<prefix>`.tsv: dense information as printed to STDOUT (see above)
-   `<prefix>`.json: comprehensive results and information on each single plasmid contig.
All files are prefixed (`<prefix>`) as the input genome fasta file.

## Installation
Bakta can be installed/used in 3 different ways.

In all cases, the custom database must be downloaded which we provide for download:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3349651.svg)](https://doi.org/10.5281/zenodo.3349651)

### BioConda
1.  install Bakta via Conda
2.  download & extract the database

Example:
```
$ conda install -c conda-forge -c bioconda -c defaults bakta
$ wget https://zenodo.org/record/3358926/files/db.tar.gz
$ tar -xzf db.tar.gz
$ rm db.tar.gz
$ bakta --db ./db genome.fasta
```

### GitHub
1.  clone the the repository
2.  download & extract the database

Example:
```
$ git clone git@github.com:oschwengers/bakta.git
$ wget <XYZ>/db.tar.gz
$ tar -xzf db.tar.gz
$ rm db.tar.gz
$ bakta/bin/bakta --db ./db genome.fasta
```

### Pip
1.  install Bakta per pip
2.  download and extract the database
3.  install 3rd party binaries

Bakta/database (1./2.):
```
$ pip3 install bakta
$ wget https://zenodo.org/record/3358926/files/db.tar.gz
$ tar -xzf db.tar.gz
$ rm db.tar.gz
$ bakta --db ./db genome.fasta
```

3rd party dependencies on Ubuntu (3.):
```
$ sudo apt install ncbi-blast+ prodigal infernal hmmer diamond-aligner
```
If there are any issues compiling ghostz, please make sure you have everything
correctly setup, e.g. `$ sudo apt install build-essential`.

## Annotation workflow
### RNA
1.  tRNA: tRNAscan-SE 2.0
2.  tmRNA: Aragorn
3.  rRNA: Infernal vs. Rfam rRNA covariance models
4.  ncRNA: Infernal vs. Rfam ncRNA covariance models
5.  CDS: Prodigal
        a) IPS vs. UniRef100
        b) PSC vs. UniRef90
6.  sORFs: IPS vs. UniRef100

### Coding sequences
The structural annotation is conducted via Prodigal and is complemented by a
custom detection of short open reading freames (**sORF**) smaller than 30 aa.

CDS:
1. Prediction via Prodigal
2. Detection of unique protein sequences (**UPS**)s via **MD5** hashes and lookup of related (**IPS**)s.
3. Alignment of remainder via Diamond vs. UniProt's UniRef90 based protein sequence clusters (**PSC**)s.

sORFs:
1. Strict filtering by overlaps with other detected features.
2. Detection of **UPS**s via **MD5** hashes and lookup of related **IPS**s.

Only **sORF** which are detected by their identity (100% coverage & 100% sequence identity)
will be included in the annotation.


## Database
The Bakta database is built on **IPS**s and **PSC**s from:
- **IPS**: UniProt UniRef100
- **PSC**: UniProt UniRef90
- 
which have been comprehensively pre-annotated integrating
annotations & database cross references (db xrefs) from:
- NCBI nonredundant proteins ('WP_*', exact matches)
- NCBI COG db (80% coverage & 90% identity)
- GO terms (via SwissProt)
- EC (via SwissProt)
- NCBI AMRFinderPlus (**IPS** exact matches, **PSC** HMM hits achieving the trusted cutoff)
- ISFinder db (90% coverage & 99% identity)

For the sake of performance and data size, all pre-annotation information is stored in a compact 
SQLite database which can be downloaded here:
[![DOI](https://zenodo.org/badge/DOI/<DOI>.svg)](https://doi.org/<DOI>)
-   [<DB_URL>](<DB_URL>)

DB size:
- zipped: x.y Gb
- unzipped x.y Gb


## Usage
Usage:
```

```

## Examples
Simple:
```
$ bakta -db ~/db genome.fasta
```

Expert: writing results to `results` directory with verbose output using 8 threads:
```
$ bakta -db ~/db --output results/ --verbose --threads 8 genome.fasta
```

## Dependencies
Bakta was developed and tested in Python 3.5 and depends on BioPython (>=1.71).

Additionally, it depends on the following 3rd party executables:
-   Prodigal (2.6.3) <https://dx.doi.org/10.1186%2F1471-2105-11-119> <https://github.com/hyattpd/Prodigal>
-   Diamond (2.0.2) <https://doi.org/10.1038/nmeth.3176> <https://github.com/bbuchfink/diamond>
-   Blast+ (2.7.1) <https://doi.org/10.1016/s0022-2836(05)80360-2> <https://blast.ncbi.nlm.nih.gov>
-   HMMER (3.2.1) <https://dx.doi.org/10.1093%2Fnar%2Fgkt263> <http://hmmer.org/>
-   INFERNAL (1.1.2) <https://dx.doi.org/10.1093%2Fbioinformatics%2Fbtt509> <http://eddylab.org/infernal>

## Citation
A manuscript is in preparation... stay tuned!
To temporarily cite our work, please transitionally refer to:
> Schwengers O., Goesmann A. (2019) Bakta: comprehensive annotation of bacterial genomes. GitHub https://github.com/oschwengers/bakta

As Bakta takes advantage of ISFinder database, please also cite:
> <ISFinder_citation>

## Issues
If you run into any issues with Bakta, we'd be happy to hear about it!
Please, execute bakta in verbose mode (`-v`) and do not hesitate
to file an issue including as much of the following as possible:
-  a detailed description of the issue
-  command line output
-  log file (`<prefix>.log`)
-  result file (`<prefix>.json`) _if possible_
-  a reproducible example of the issue with an input file that you can share _if possible_
