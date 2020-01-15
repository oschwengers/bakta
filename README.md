[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/oschwengers/bacanno/blob/master/LICENSE)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/bakta.svg)
![GitHub release](https://img.shields.io/github/release/oschwengers/bakta.svg)
![PyPI](https://img.shields.io/pypi/v/cb-bakta.svg)
![PyPI - Status](https://img.shields.io/pypi/status/bakta.svg)
![Conda](https://img.shields.io/conda/v/bioconda/bakta.svg)
![Conda](https://img.shields.io/conda/pn/bioconda/bakta.svg)

# Bakta: comprehensive annotation of bacterial genomes.

## Contents
-   [Description](#description)
-   [Input/Output](#inputoutput)
-   [Installation](#installation)
-   [Usage](#usage)
-   [Examples](#examples)
-   [Database](#database)
-   [Dependencies](#dependencies)
-   [Citation](#citation)

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


### GitHub
1.  clone the the repository
2.  download & extract the database

Example:
```
$ git clone git@github.com:oschwengers/bacanno.git
$ wget <XYZ>/db.tar.gz
$ tar -xzf db.tar.gz
$ rm db.tar.gz
$ bacanno/bin/bacanno --db ./db genome.fasta
```

### Conda
1.  install Bakta via Conda
2.  download & extract the database

Example:
```
$ conda install -c conda-forge -c bioconda -c defaults bacanno
$ wget https://zenodo.org/record/3358926/files/db.tar.gz
$ tar -xzf db.tar.gz
$ rm db.tar.gz
$ bacanno --db ./db genome.fasta
```

### Pip
1.  install Bakta per pip
2.  download and extract the database
3.  install 3rd party binaries

Bakta/database (1./2.):
```
$ pip3 install cb-bacanno
$ wget https://zenodo.org/record/3358926/files/db.tar.gz
$ tar -xzf db.tar.gz
$ rm db.tar.gz
$ bacanno --db ./db genome.fasta
```

3rd party dependencies on Ubuntu (3.):
```
$ sudo apt install ncbi-blast+ prodigal infernal hmmer
$ wget http://www.bi.cs.titech.ac.jp/ghostz/releases/ghostz-1.0.2.tar.gz
$ tar -xzf ghostz-1.0.2.tar.gz
$ cd ghostz-1.0.2/
$ make
$ sudo cp ghostz /usr/bin/
```
If there are any issues compiling ghostz, please make sure you have everything
correctly setup, e.g. `$ sudo apt install build-essential`.

## Usage
Usage:
```

```

## Examples
Simple:
```
$ bacanno -db ~/db genome.fasta
```

Expert: writing results to `results` directory with verbose output using 8 threads:
```
$ bacanno -db ~/db --output results/ --verbose --threads 8 genome.fasta
```

## Database
Bakta depends on a custom database based on:
-   UniProt UniRef100
-   UniProt UniRef90

comprising annotations from:
-   NCBI non-redundant proteins
-   NCBI RefSeq
-   NCBI AMRFinderPlus
-   COG
-   GO
-   E.C.
-   ISFinder db

This database can be downloaded here:
(zipped x.y Gb, unzipped x.y Gb)
[![DOI](https://zenodo.org/badge/DOI/<DOI>.svg)](https://doi.org/<DOI>)
-   [<DB_URL>](<DB_URL>)

## Dependencies
Bakta was developed and tested in Python 3.5 and depends on BioPython (>=1.71).

Additionally, it depends on the following 3rd party executables:
-   Prodigal (2.6.3) <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2848648> <https://github.com/hyattpd/Prodigal>
-   Ghostz (1.0.2) <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4393512> <http://www.bi.cs.titech.ac.jp/ghostz>
-   Blast+ (2.7.1) <https://www.ncbi.nlm.nih.gov/pubmed/2231712> <https://blast.ncbi.nlm.nih.gov>
-   HMMER (3.2.1) <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3695513/> <http://hmmer.org/>
-   INFERNAL (1.1.2) <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3810854> <http://eddylab.org/infernal>

## Citation
A manuscript is in preparation... stay tuned!
To temporarily cite our work, please transitionally refer to:
> Schwengers O., Goesmann A. (2019) Bakta: comprehensive annotation of bacterial genomes. GitHub https://github.com/oschwengers/bacanno

As Bakta takes advantage of ISFinder database, please also cite:
> <ISFinder_citation>

## Issues
If you run into any issues with Bakta, we'd be happy to hear about it!
Please, start the pipeline with `-v` (verbose) and do not hesitate
to file an issue including as much of the following as possible:
-   a detailed description of the issue
-   the bacanno cmd line output
-   the `<prefix>.json` file if possible
-   if possible a reproducible example of the issue with an input file that you can share
(helps us identify whether the issue is specific to a particular computer, operating system, and/or dataset).
