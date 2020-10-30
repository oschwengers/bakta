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
Bakta is a local tool dedicated to the rapid & comprehensive annotation of bacteria & plasmids. It provides **dbxref**-rich and **sORF**-including annotations as well as machine-readble (`JSON`) & standard output files for automatic downstream analysis.

The annotation of microbial genomes is a diverse task comprising the structural & functional annotation of different feature types with distinct overlapping characteristics. Existing local annotation pipelines cover a broad range of microbial taxa, *e.g.* bacteria, aerchaea, viruses. To expand the range of supported feature types and to improve annotation quality and database cross-references, Bakta is dedicated to the annotation of bacteria and plasmids.

Knwon exact protein coding sequences (**CDS**) are identified via `MD5` aa digests and annotated with unique database cross-references (**dbxref**) to:

- RefSeq (`WP_*`)
- UniRef90/UniRef100 (`UniRef100_*`)
- UniParc (`UPI*`)

By this, identified exact sequences allow the surveillance of certain gene alleles and streamlining comparative analysis. Also, posterior (external) annotations reagarding `putative` & `hypothetical` protein sequences can easily be mapped back to existing `cds` via these exact identifiers (*E. coli* gene [ymiA](https://www.uniprot.org/uniprot/P0CB62) [...more](https://www.uniprot.org/help/dubious_sequences)).
Additionally, **CDS** are annotated via UniRef90 protein clusters. These as well as exact sequences are further annotated (`GO`, `COG`, `EC`).

Next to standard feature types (tRNA, tmRNA, rRNA, ncRNA, CRISPR, CDS, gaps) Bakta also detects and annotates:

- short ORFs (**sORF**) which are not predicted by tools like `Prodigal`
- ncRNA regulatory regions distinct from ncRNA genes
- origins of replication/transfer (oriC, oriV, oriT)

Bakta can annotate a typical bacterial genome within minutes and hence fits the niche between large & computationally-demanding (online) pipelines and rapid, highly-customizable offline tools like Prokka. If Bakta doesn't fit your needs, please consider using [Prokka](https://github.com/tseemann/prokka). The development of Bakta was highly inspired by Prokka and many command line options are mutually compatible for the sake of interoperability and user convenience.

## Input/Output

### Input

Bakta accepts bacterial and plasmid assemblies (complete / draft) in fasta format.

### Output

Bakta provides detailed information on each annotated feature in a standardized machine-readable JSON format.
In addition, Bakta supports the following standard file formats:

- `tsv`: dense human readable information as simple tab separated values
- `GFF3`: standard GFF3 format
- `GenBank`: standard GenBank format
- `fna`: genome sequences as FASTA
- `faa`: protein sequences as FASTA

## Examples

Simple:

```bash
$ bakta -db ~/db genome.fasta
```

Expert: verbose output writing results to *results* directory (`TSV`, `GFF3` and `GenBank`) with *ecoli123* file `prefix` and *eco634* `locus tag` using an existing prodigal training file and 8 threads:

```bash
$ bakta -db ~/db --verbose --output results/ --tsv --gff --genbank --prefix ecoli123 --locus-tag eco634 --prodigal-tf eco.tf --threads 8 genome.fasta
```

## Installation

Bakta can be installed via BioConda, Pip and GitHub.
To automatically install all required 3rd party dependencies, we encourage to use BioConda.

### BioConda

```bash
$ conda install -c bioconda bakta
```

### Pip

1. install Bakta per pip
2. install 3rd party binaries (-> Dependencies)

```bash
$ python3 -m pip install bakta
```

### GitHub

1. clone the the repository
2. install Python dependencies
3. install 3rd party binaries (-> Dependencies)

```bash
$ git clone git@github.com:oschwengers/bakta.git
$ python3 setup.py install --user
```

### Mandatory database

In all cases, Bakta requires a mandatory database which is publicly hosted at Zenodo:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXX.svg)](https://doi.org/10.5281/zenodo.XXX)
Further information is provided [below](#database).

```bash
$ wget <XYZ>/db.tar.gz
$ tar -xzf db.tar.gz
$ rm db.tar.gz
```

The, the database path can be provided via the `--db` parameter:

```bash
$ bakta --db <db-path>
```

It's also possible to set a `BAKTA_DIR` environment variable:

```bash
$ export BAKTA_DIR=<db-path>
```

Additionally, for a system-wide setup, the database can be copied to the Bakta base directory:

```bash
$ cp -r db/ <bakta-installation-dir>
```

### Dependencies

Bacta uses the following 3rd party executables which must be installed & executable:

- tRNAscan-SE (2.0.6) <https://doi.org/10.1101/614032> <http://lowelab.ucsc.edu/tRNAscan-SE>
- Aragorn (1.2.38) <http://dx.doi.org/10.1093/nar/gkh152> <http://130.235.244.92/ARAGORN>
- INFERNAL (1.1.3) <https://dx.doi.org/10.1093%2Fbioinformatics%2Fbtt509> <http://eddylab.org/infernal>
- PILER-CR (1.06) <https://doi.org/10.1186/1471-2105-8-18> <http://www.drive5.com/pilercr>
- Prodigal (2.6.3) <https://dx.doi.org/10.1186%2F1471-2105-11-119> <https://github.com/hyattpd/Prodigal>
- Diamond (2.0.2) <https://doi.org/10.1038/nmeth.3176> <https://github.com/bbuchfink/diamond>
- Blast+ (2.7.1) <https://www.ncbi.nlm.nih.gov/pubmed/2231712> <https://blast.ncbi.nlm.nih.gov>

On Ubuntu you can install these via:

```bash
$ sudo apt install aragorn infernal prodigal diamond-aligner ncbi-blast+
```

tRNAscan-se must be installed manually as v2.0 is currently not yet available via standard Ubuntu packages.

## Annotation workflow

### RNAs

1. tRNA genes: tRNAscan-SE 2.0
2. tmRNA genes: Aragorn
3. rRNA genes: Infernal vs. Rfam rRNA covariance models
4. ncRNA genes: Infernal vs. Rfam ncRNA covariance models
5. ncRNA regulatory regions: Infernal vs. Rfam ncRNA covariance models
6. CRISPR arrays: PILER-CR

Bakta distinguishes ncRNA genes and (regulatory) regions in order to enable the distinct handling thereof during the annotation process, *i.e.* feature overlap detection.
ncRNA gene types:

- sRNA
- antisense
- ribozyme
- antitoxin

ncRNA (regulatory) region types:

- riboswitch
- thermoregulator
- leader
- frameshift element

### Coding sequences

The structural prediction is conducted via Prodigal and complemented by a custom detection of short open reading freames (**sORF**) < 30 aa.

To rapidly conduct a comprehensive annotation while also identifing known protein sequences with exact sequence matches,
Bakta uses a comprehensive SQLite database comprising protein sequence digests and pre-annotations for millions of
known protein sequences and clusters.

Conceptual terms:

- **UPS**: unique protein sequences identified via length and **MD5** sequence digests (100% coverage & 100% sequence identity)
- **IPS**: identical protein sequences comprising representatives of UniProt's UniRef100 protein sequence clusters
- **PSC**: protein sequences clusters comprising representatives of UniProt's UniRef90 protein sequence clusters

**CDS**:

1. Prediction via Prodigal
2. Detection of **UPS**s via **MD5** digests and lookup of related **IPS** and **PCS**
3. Homology search of remainder via Diamond vs. **PSC**
4. Combination of available **IPS** & **PSC** information favouring more specific annotations and avoiding redundancy

**CDS** without **IPS** or **PSC** hits will be marked as `hypothetical`.
Additionally, all **CDS** without gene symbols or with product descriptions equal to `hypothetical` will be marked as `hypothetical`.

However, `hypothetical` **CDS** are included in the final annotation.

**sORFs**:

1. Custom detection & extraction of **sORF** with amino acid lengths < 30 aa
2. Filter via strict feature type-dependent overlap filters with annotated features
3. Detection of **UPS** via **MD5** hashes and lookup of related **IPS**
4. Homology search of remainder via Diamond vs. seed sequences of an sORF subset of UniProt's UniRef90 **PSC**
5. Exclude **sORF** without sufficient annotation information

**sORF** not identified via **IPS** or **PSC** will be discarded.
Additionally, all **sORF** without gene symbols or with product descriptions equal to `hypothetical` will be discarded.

Due due to uncertain nature of **sORF** prediction, only those identified via **IPS** / **PSC** hits exhibiting proper gene symbols or product descriptions different from `hypothetical` will be included in the final annotation.

### Miscellaneous

1. Gaps: in-mem detection & annotation of sequence gaps
2. oriC/oriV/oriT: Blast+ (blastn) vs. MOB-suite oriT & DoriC oriC/oriV sequences. Annotations of ori regions take into account overlapping Blast+ hits and are conducted based on a majority vote heuristic.

## Database

The Bakta database comprises a set of DNA & AA sequence databases as well as HMM & covariance models.
In addition, at its core Bakta uses a compact SQLite db storing protein sequence digests, lengths, pre-annotations and *dbxrefs* of **UPS**, **IPS** and **PSC** from:

- **UPS**: UniParc / UniProtKB
- **IPS**: UniProt UniRef100
- **PSC**: UniProt UniRef90

This allows the exact protein sequences identification via **MD5** digests & sequence lengths as well as the rapid subsequent lookup of related information.
**IPS** & **PSC** have been comprehensively pre-annotated integrating annotations & database *dbxrefs* from:

- NCBI nonredundant proteins ('WP_*', exact matches)
- NCBI COG db (80% coverage & 90% identity)
- GO terms (via **IPS**/**PSC** SwissProt entries)
- EC (via **IPS**/**PSC** SwissProt entries)
- NCBI AMRFinderPlus (**IPS** exact matches, **PSC** HMM hits reaching trusted cutoffs)
- ISFinder db (90% coverage & 99% identity)

Database (23 Gb zipped, 43 Gb unzipped) hosted at Zenodo:
[![DOI](https://zenodo.org/badge/DOI/<DOI>.svg)](https://doi.org/<DOI>)

## Usage

Usage:

```bash
bakta --help
usage: bakta [--db DB] [--min-contig-length MIN_CONTIG_LENGTH]
             [--prefix PREFIX] [--output OUTPUT] [--tsv] [--gff3] [--genbank]
             [--embl] [--fna] [--faa] [--genus GENUS] [--species SPECIES]
             [--strain STRAIN] [--plasmid PLASMID] [--prodigal-tf PRODIGAL_TF]
             [--translation-table {11,4}] [--complete] [--gram {+,-,?}]
             [--locus LOCUS] [--locus-tag LOCUS_TAG] [--keep-contig-headers]
             [--replicons REPLICONS] [--skip-trna] [--skip-tmrna]
             [--skip-rrna] [--skip-ncrna] [--skip-ncrna-region]
             [--skip-crispr] [--skip-cds] [--skip-sorf] [--skip-gap]
             [--skip-ori] [--help] [--verbose] [--threads THREADS]
             [--tmp-dir TMP_DIR] [--version] [--citation]
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
  --fna                 Write genome sequences as fasta file
  --faa                 Write translated CDS sequences as fasta file

Organism:
  --genus GENUS         Genus name
  --species SPECIES     Species name
  --strain STRAIN       Strain name
  --plasmid PLASMID     Plasmid name

Annotation:
  --prodigal-tf PRODIGAL_TF
                        Path to existing Prodigal training file to use for CDS
                        prediction
  --translation-table {11,4}
                        Translation table to use: 11/4 (default = 11)
  --complete            Replicons (chromosome/plasmid[s]) are complete
  --gram {+,-,?}        Gram type: +/-/? (default = '?')
  --locus LOCUS         Locus prefix (instead of 'contig')
  --locus-tag LOCUS_TAG
                        Locus tag prefix
  --keep-contig-headers
                        Keep original contig headers
  --replicons REPLICONS, -r REPLICONS
                        Replicon information table (TSV)

Workflow:
  --skip-trna           Skip tRNA detection & annotation
  --skip-tmrna          Skip tmRNA detection & annotation
  --skip-rrna           Skip rRNA detection & annotation
  --skip-ncrna          Skip ncRNA detection & annotation
  --skip-ncrna-region   Skip ncRNA region detection & annotation
  --skip-crispr         Skip CRISPR array detection & annotation
  --skip-cds            Skip CDS detection & annotation
  --skip-sorf           Skip sORF detection & annotation
  --skip-gap            Skip gap detection & annotation
  --skip-ori            Skip oriC/oriT detection & annotation

General:
  --help, -h            Show this help message and exit
  --verbose, -v         Print verbose information
  --threads THREADS, -t THREADS
                        Number of threads to use (default = number of
                        available CPUs)
  --tmp-dir TMP_DIR     Location for temporary files (default = system
                        dependent auto detection)
  --version             show program's version number and exit
  --citation            Print citation
```

## Citation

A manuscript is in preparation. To temporarily cite our work, please transitionally refer to:
> Schwengers O., Goesmann A. (2020) Bakta: comprehensive annotation of bacterial genomes. GitHub https://github.com/oschwengers/bakta

Bakta takes advantage of many publicly available databases. If you find any of the data used within Bakta useful, please also be sure to credit the primary source also:
- UniProt: <https://doi.org/10.1093/nar/gky1049>
- RefSeq: <DOI>
- Rfam: <https://doi.org/10.1002/cpbi.51>
- AMRFinder: <https://doi.org/10.1128/AAC.00483-19>
- ISFinder: <https://doi.org/10.1093/nar/gkj014>
- AntiFam: <DOI>
- Mob-suite: <https://doi.org/10.1099/mgen.0.000206>
- COG: <DOI>

## Issues and Feature Requests

If you run into any issues with Bakta, we'd be happy to hear about it!
Please, execute bakta in verbose mode (`-v`) and do not hesitate
to file an issue including as much information as possible:

- a detailed description of the issue
- command line output
- log file (`<prefix>.log`)
- result file (`<prefix>.json`) _if possible_
- a reproducible example of the issue with an input file that you can share _if possible_
