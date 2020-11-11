[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/oschwengers/bakta/blob/master/LICENSE)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/bakta.svg)
![GitHub release](https://img.shields.io/github/release/oschwengers/bakta.svg)
![PyPI](https://img.shields.io/pypi/v/bakta.svg)
![PyPI - Status](https://img.shields.io/pypi/status/bakta.svg)
![Conda](https://img.shields.io/conda/v/bioconda/bakta.svg)
![Conda](https://img.shields.io/conda/pn/bioconda/bakta.svg)

# Bakta: rapid & comprehensive annotation of bacterial genomes & plasmids

## Contents

- [Description](#description)
- [Input/Output](#inputoutput)
- [Examples](#examples)
- [Installation](#installation)
  - [Bioconda](#bioconda)
  - [Docker](#docker)
  - [Pip](#pip)
  - [Dependencies](#dependencies)
  - [Database](#mandatory_database)
- [Annotation workflow](#annoation_workflow)
- [Database](#database)
- [Usage](#usage)
- [Citation](#citation)
- [FAQ](#faq)
- [Issues & Feature Requests](#issues)

## Description

**TL;DR**

Bakta is an offline tool dedicated to the rapid & comprehensive annotation of bacteria & plasmids. It provides **dbxref**-rich and **sORF**-including annotations in machine-readble (`JSON`) & bioinformatics standard file formats for automatic downstream analysis.

The annotation of microbial genomes is a diverse task comprising the structural & functional annotation of different feature types with distinct overlapping characteristics. Existing local annotation pipelines cover a broad range of microbial taxa, *e.g.* bacteria, aerchaea, viruses. To streamline and foster the expansion of supported feature types, Bakta is strictly dedicated to the annotation of bacteria and plasmids. To standardize the annotation of bacterial sequences, Bakta uses a comprehensive annotation database based on UniProt's UniRef protein clusters enriched by cross-references and specialized niche databases.

Exact matches to known protein coding sequences (**CDS**), subsequently referred to as identical protein sequences (**IPS**) are identified via `MD5` digests and annotated with database cross-references (**dbxref**) to:

- RefSeq (`WP_*`)
- UniRef100/UniRef90 (`UniRef100_*`/`UniRef90_*`)
- UniParc (`UPI*`)

By doing so, **IPS** allow the surveillance of distinct gene alleles and streamline comparative analysis. Also, posterior (external) annotations of `putative` & `hypothetical` protein sequences can be mapped back to existing `cds` via these exact & stable identifiers (*E. coli* gene [ymiA](https://www.uniprot.org/uniprot/P0CB62) [...more](https://www.uniprot.org/help/dubious_sequences)).
Unidentified remaining **CDS** are annotated via UniRef90 protein sequence clusters (**PSC**).
**PSC** & **IPS** are enriched by pre-annotated and stored information (`GO`, `COG`, `EC`).

Next to standard feature types (tRNA, tmRNA, rRNA, ncRNA, CRISPR, CDS, gaps) Bakta also detects and annotates:

- short ORFs (**sORF**) which are not predicted by tools like `Prodigal`
- ncRNA regulatory regions distinct from ncRNA genes
- origins of replication/transfer (oriC, oriV, oriT)

Bakta can annotate a typical bacterial genome within minutes and hence fits the niche between large & computationally-demanding (online) pipelines and rapid, highly-customizable offline tools like Prokka. If Bakta does not fit your needs, please consider using [Prokka](https://github.com/tseemann/prokka). The development of Bakta was highly inspired by Prokka and many command line options are mutually compatible for the sake of interoperability and user convenience.

## Input/Output

### Input

Bakta accepts bacterial and plasmid assemblies (complete / draft) in (zipped) fasta format.

Further genome information and workflow customizations can be provided and set via a number of input parameters.
For a full description, please have a look at the [Usage](#usage) section.

Most important parameters:

- use a custom database location, *e.g.* a local instance for runtime improvements: `--db`
- genome parameters: `--min-contig-length`, `--complete`
- number of threads: `--threads`
- locus information `--locus`, `--locus-tag`

Replicon meta data table:

To fine-tune the very details of each sequence in the input fasta file, Bakta accepts a replicon meta data table provided in `tsv` file format: `--replicons <tsv-replicon-file>`.
Thus, for example, complete replicons within partially completed draft assemblies can be marked as such.

Table format:

original locus id  |  new locus id  |  type  |  topology  |  name
----|----------------|----------------|----------------|----------------
`old id` | [`new id` / `<empty>`] | [`chromosome` / `plasmid` / `contig` / `<empty>`] | [`circular` / `linear` / `<empty>`] | `name`

Thus, for each input sequence recognized via the `original locus id` a `new locus id`, the replicon `type` and the `topology` as well a `name` can be explicitly set.

Available short cuts:

- `chromosome`: `c`
- `plasmid`: `p`
- `circular`: `c`
- `linear`: `l`

`<empty>` values (`-` / ``) will be replaced by defaults. If **new locus id** is `empty`, a new contig name will be autogenerated.

Defaults:
replicon type: `contig`
topology: `linear`

Example:

original locus id  |  new locus id  |  type  |  topology  |  name
----|----------------|----------------|----------------|----------------
NODE_1 | chrom | `chromosome` | `circular` | `-`
NODE_2 | p1 | `plasmid` | `c` | `pXYZ1`
NODE_3 | p2 | `p`  |  `c` | `pXYZ2`
NODE_4 | special-contig-name-xyz |  `-` | -
NODE_5 | `` |  `-` | -

### Output

Bakta provides detailed information on each annotated feature in a standardized machine-readable JSON file.
In addition, the following standard file formats are supported:

- `tsv`: annotations as simple human readble tab separated values
- `GFF3`: annotations in GFF3 format
- `GenBank`: annotations in GenBank format
- `fna`: replicons/contigs as FASTA
- `faa`: CDS as FASTA

## Examples

Simple:

```bash
$ bakta --db ~/db genome.fasta
```

Expert: verbose output writing results to *results* directory with *ecoli123* file `prefix` and *eco634* `locus tag` using an existing prodigal training file, using additional replicon information and 8 threads:

```bash
$ bakta --db ~/db --verbose --output results/ --prefix ecoli123 --locus-tag eco634 --prodigal-tf eco.tf --replicons replicon.tsv --threads 8 genome.fasta
```

## Installation

Bakta can be installed via BioConda, Docker or Pip.
To automatically install all required 3rd party dependencies, we highly encourage to use [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).
In all cases a mandatory db must be downloaded (-> Mandatory database)

### BioConda

```bash
$ conda install -c conda-forge -c bioconda -c defaults bakta
```

### Docker
```bash
$ sudo docker pull oschwengers/bakta

$ sudo docker run oschwengers/bakta --help
$ bakta-docker.sh --help
```

### Pip

1. install Bakta per pip
2. install 3rd party binaries (-> Dependencies)

```bash
$ python3 -m pip install --user bakta
```

### Dependencies

Bacta requires Biopython (>=1.72), Xopen (0.9) and the following 3rd party executables which must be installed & executable:

- tRNAscan-SE (2.0.6) <https://doi.org/10.1101/614032> <http://lowelab.ucsc.edu/tRNAscan-SE>
- Aragorn (1.2.38) <http://dx.doi.org/10.1093/nar/gkh152> <http://130.235.244.92/ARAGORN>
- INFERNAL (1.1.2) <https://dx.doi.org/10.1093%2Fbioinformatics%2Fbtt509> <http://eddylab.org/infernal>
- PILER-CR (1.06) <https://doi.org/10.1186/1471-2105-8-18> <http://www.drive5.com/pilercr>
- Prodigal (2.6.3) <https://dx.doi.org/10.1186%2F1471-2105-11-119> <https://github.com/hyattpd/Prodigal>
- Hmmer (3.3.1) <https://doi.org/10.1093/nar/gkt263> <http://hmmer.org>
- Diamond (2.0.2) <https://doi.org/10.1038/nmeth.3176> <https://github.com/bbuchfink/diamond>
- Blast+ (2.7.1) <https://www.ncbi.nlm.nih.gov/pubmed/2231712> <https://blast.ncbi.nlm.nih.gov>

On Ubuntu/Debian/Mint you can install these via:

```bash
$ sudo apt install aragorn infernal prodigal diamond-aligner ncbi-blast+
```

tRNAscan-se must be installed manually as v2.0 is currently not yet available via standard Ubuntu packages.

### Mandatory database

In all cases, Bakta requires a mandatory database which is publicly hosted at Zenodo:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4247253.svg)](https://doi.org/10.5281/zenodo.4247253)
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

It's also possible to set a `BAKTA_DB` environment variable:

```bash
$ export BAKTA_DB=<db-path>
```

Additionally, for a system-wide setup, the database can be copied to the Bakta base directory:

```bash
$ cp -r db/ <bakta-installation-dir>
```

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
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4247253.svg)](https://doi.org/10.5281/zenodo.4247253)

## Usage

Usage:

```bash
bakta --help
usage: bakta [--db DB] [--min-contig-length MIN_CONTIG_LENGTH]
             [--prefix PREFIX] [--output OUTPUT] [--genus GENUS]
             [--species SPECIES] [--strain STRAIN] [--plasmid PLASMID]
             [--prodigal-tf PRODIGAL_TF] [--translation-table {11,4}]
             [--complete] [--gram {+,-,?}] [--locus LOCUS]
             [--locus-tag LOCUS_TAG] [--keep-contig-headers]
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
- RefSeq: <https://doi.org/10.1093/nar/gkx1068>
- Rfam: <https://doi.org/10.1002/cpbi.51>
- AMRFinder: <https://doi.org/10.1128/AAC.00483-19>
- ISFinder: <https://doi.org/10.1093/nar/gkj014>
- AntiFam: <https://doi.org/10.1093/database/bas003>
- Mob-suite: <https://doi.org/10.1099/mgen.0.000206>
- COG: <https://doi.org/10.1093/bib/bbx117>

## FAQ

* __Bakta is running too long without CPU load... why?__
Bakta takes advantage of an SQLite DB which results in high storage IO loads. If this DB is stored on a remote / network volume, the lookup of IPS/PSC annotations might take a long time. In these cases, please, consider moving the DB to a local volume/hard drive.

## Issues and Feature Requests

If you run into any issues with Bakta, we'd be happy to hear about it!
Please, execute bakta in verbose mode (`-v`) and do not hesitate
to file an issue including as much information as possible:

- a detailed description of the issue
- command line output
- log file (`<prefix>.log`)
- result file (`<prefix>.json`) _if possible_
- a reproducible example of the issue with an input file that you can share _if possible_
