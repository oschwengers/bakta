[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/oschwengers/bakta/blob/master/LICENSE)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/bakta.svg)
![GitHub release](https://img.shields.io/github/release/oschwengers/bakta.svg)
![PyPI](https://img.shields.io/pypi/v/bakta.svg)
![PyPI - Status](https://img.shields.io/pypi/status/bakta.svg)
![Conda](https://img.shields.io/conda/v/bioconda/bakta.svg)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4247252.svg)](https://doi.org/10.5281/zenodo.4247252)

# Bakta: Rapid & standardized annotation of bacterial genomes & plasmids

Bakta is a tool for the rapid & standardized local annotation of bacterial genomes & plasmids. It provides **dbxref**-rich and **sORF**-including annotations in machine-readble `JSON` & bioinformatics standard file formats for automatic downstream analysis.

> Bakta is young and still in a beta status! So please, test it, bend, sequeeze & smash it, try to crash it... but don't forget to file bug reports ;-)
Of course, feedback of any kind and feature requests are highly welcome & very much appreciated!

## Contents

- [Description](#description)
- [Installation](#installation)
- [Examples](#examples)
- [Input/Output](#inputoutput)
- [Usage](#usage)
- [Annotation workflow](#annoation_workflow)
- [Database](#database)
- [WIP](#wip)
- [Citation](#citation)
- [FAQ](#faq)
- [Issues & Feature Requests](#issues)

## Description

- **Bacteria & plasmids only**
Bakta is desined to annotate bacteria and plasmids, only. This decision by design has been made in order to tweak the annotation process regarding tools, preferences & databases and to streamline further development & maintenance of the software.

- **FAIR annotations**
To provide standardized annotations adhearing to [FAIR](https://www.go-fair.org/fair-principles) principles, Bakta utilizes a comprehensive & versioned custom annotation database based on UniProt's [UniRef100 & UniRef90](https://www.uniprot.org/uniref/) protein clusters (`FAIR` -> [DOI](http://dx.doi.org/10.1038/s41597-019-0180-9)/[DOI](https://doi.org/10.1093/nar/gkaa1100)) enriched with dbxrefs (`GO`, `COG`, `EC`) and annotated by specialized niche databases. For each db version we provide a comprehensive log file of imported sequences and conducted annotations.

- **Protein sequence identification**
Fostering the FAIR aspect, Bakta identifies identical protein sequences (**IPS**) via `MD5` digests which are annotated with database cross-references (**dbxref**) to RefSeq (`WP_*`), UniRef100 (`UniRef100_*`) and UniParc (`UPI*`).
By doing so, IPS allow the surveillance of distinct gene alleles and streamlining comparative analysis as well as posterior (external) annotations of `putative` & `hypothetical` protein sequences which can be mapped back to existing CDS via these exact & stable identifiers (*E. coli* gene [ymiA](https://www.uniprot.org/uniprot/P0CB62) [...more](https://www.uniprot.org/help/dubious_sequences)). Currently, Bakta identifies ~169 mio distinct UniRef100 sequences and for certain genomes, up to 99 % of all CDS can be identified this way, sparing expensive homology searches.

- **Short open reading frames**
Next to standard feature types (tRNA, tmRNA, rRNA, ncRNA, ncRNA cis-regulatory regions, CRISPR, CDS, oriC/V/T and gaps) Bakta also detects and annotates short open reading frames (**sORF**) which are not predicted by tools like `Prodigal`.

- **Fast**
Bakta can annotate a typical bacterial genome in 10 &plusmn;5 min on a laptop, plasmids in a couple of seconds/minutes.

- **Reasoning**
Annotating bacterial genomes in a standardized, taxon-independent, high-throughput and local manner, Bakta targets the niche between standardized, fully-featured but computationally-demanding pipelines like [PGAP](https://github.com/ncbi/pgap) and rapid highly-customizable offline tools like [Prokka](https://github.com/tseemann/prokka). Indeed, Bakta is heavily inspired by Prokka (kudos to [Torsten Seemann](https://github.com/tseemann)) and many command line options are compatible for the sake of interoperability and user convenience. Hence, if Bakta does not fit your needs, please try Prokka.

## Installation

Bakta can be installed via BioConda, Docker, Singularity and Pip.
However, to automatically install all required 3rd party dependencies, we encourage to use [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) or [Docker](https://www.docker.com/get-started)/[Singularity](https://sylabs.io/singularity).
In all cases the mandatory [database](#database_download) must be downloaded.

### BioConda

```bash
$ conda install -c conda-forge -c bioconda -c defaults bakta
```

### Docker

```bash
$ sudo docker pull oschwengers/bakta
$ sudo docker run oschwengers/bakta --help
```

Installation instructions, get-started and guides: Docker [docs](https://docs.docker.com)
For further convenience, we provide a shell script (`bakta-docker.sh`) handling stuff like volume mounting, etc:

```bash
$ bakta-docker.sh --help
```

### Singularity

```bash
$ singularity pull docker://oschwengers/bakta:latest
$ singularity run bakta-latest.simg --help
```

Installation instructions, get-started and guides: Singularity [docs](https://sylabs.io/docs)

### Pip

```bash
$ python3 -m pip install --user bakta
```

Bacta requires the following 3rd party executables which must be installed & executable:

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
$ sudo apt install trnascan-se aragorn infernal pilercr prodigal hmmer diamond-aligner ncbi-blast+
```

Tested with Ubuntu 20.04 - some older distributions might provide outdated versions, *e.g.* trnascan-se in Ubuntu 18.04. In these cases dependencies must be installed manually.

### Database download

Bakta requires a mandatory database which is publicly hosted at Zenodo:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4247252.svg)](https://doi.org/10.5281/zenodo.4247252)

Further information is provided in the [database](#database) section below.

```bash
$ wget https://zenodo.org/record/4247253/files/db.tar.gz
$ tar -xzf db.tar.gz
$ rm db.tar.gz
```

The db path can either be provided via parameter (`--db`) or environment variable (`BAKTA_DB`):

```bash
$ bakta --db <db-path> genome.fasta

$ export BAKTA_DB=<db-path>
$ bakta genome.fasta
```

Additionally, for a system-wide setup, the database can be copied to the Bakta base directory:

```bash
$ cp -r db/ <bakta-installation-dir>
```

## Examples

Simple:

```bash
$ bakta --db ~/db genome.fasta
```

Expert: verbose output writing results to *results* directory with *ecoli123* file `prefix` and *eco634* `locus tag` using an existing prodigal training file, using additional replicon information and 8 threads:

```bash
$ bakta --db ~/db --verbose --output results/ --prefix ecoli123 --locus-tag eco634 --prodigal-tf eco.tf --replicons replicon.tsv --threads 8 genome.fasta
```

## Input/Output

### Input

Bakta accepts bacterial and plasmid assemblies (complete / draft) in (zipped) fasta format.

For a full description of how further genome information can be provided and workflow customizations can be set, please have a look at the [Usage](#usage) section.

Replicon meta data table:

To fine-tune the very details of each sequence in the input fasta file, Bakta accepts a replicon meta data table provided in `tsv` file format: `--replicons <file.tsv>`.
Thus, complete replicons within partially completed draft assemblies can be marked & handled as such, *e.g.* detection & annotation of features spanning sequence edges.

Table format:

original sequence id  |  new sequence id  |  type  |  topology  |  name
----|----------------|----------------|----------------|----------------
`old id` | [`new id` / `<empty>`] | [`chromosome` / `plasmid` / `contig` / `<empty>`] | [`circular` / `linear` / `<empty>`] | [`name` / `<empty>`]

For each input sequence recognized via the `original locus id` a `new locus id`, the replicon `type` and the `topology` as well a `name` can be explicitly set.

Shortcuts:

- `chromosome`: `c`
- `plasmid`: `p`
- `circular`: `c`
- `linear`: `l`

`<empty>` values (`-` / ``) will be replaced by defaults. If **new locus id** is `empty`, a new contig name will be autogenerated.

Defaults:

- type: `contig`
- topology: `linear`

Example:

original locus id  |  new locus id  |  type  |  topology  |  name
----|----------------|----------------|----------------|----------------
NODE_1 | chrom | `chromosome` | `circular` | `-`
NODE_2 | p1 | `plasmid` | `c` | `pXYZ1`
NODE_3 | p2 | `p`  |  `c` | `pXYZ2`
NODE_4 | special-contig-name-xyz |  `-` | -
NODE_5 | `` |  `-` | -

### Output

Bakta provides detailed information on each annotated feature in a standardized machine-readable JSON file:

```json
{
    "genome": {
        "genus": "Escherichia",
        "species": "coli",
        ...
    },
    "stats": {
        "size": 5594605,
        "gc": 0.497,
        ...
    },
    "features": [
        {
            "type": "cds",
            ...
        },
        ...
    ],
    "sequences": [
        {
            "id": "c1",
            "description": "[organism=Escherichia coli] [completeness=complete] [topology=circular]",
            "sequence": "AGCTTT...",
            "length": 5498578,
            "complete": true,
            "type": "chromosome",
            "topology": "circular"
            ...
        },
        ...
    ]
}
```

Additionally, the following output files are written:

- `.tsv`: annotations as simple human readble tab separated values
- `.gff3`: annotations & sequences in GFF3 format
- `.gbff`: annotations & sequences in (multi) GenBank format
- `.fna`: replicon/contig DNA sequences as FASTA
- `.faa`: CDS/sORF amino acid sequences as FASTA

## Usage

Usage:

```bash
bakta --help
usage: bakta [--db DB] [--min-contig-length MIN_CONTIG_LENGTH] [--prefix PREFIX] [--output OUTPUT] [--genus GENUS] [--species SPECIES] [--strain STRAIN] [--plasmid PLASMID] [--complete] [--prodigal-tf PRODIGAL_TF] [--translation-table {11,4}] [--gram {+,-,?}] [--locus LOCUS]
             [--locus-tag LOCUS_TAG] [--keep-contig-headers] [--replicons REPLICONS] [--skip-trna] [--skip-tmrna] [--skip-rrna] [--skip-ncrna] [--skip-ncrna-region] [--skip-crispr] [--skip-cds] [--skip-sorf] [--skip-gap] [--skip-ori] [--help] [--verbose] [--threads THREADS]
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
  --complete            All sequences are complete replicons (chromosome/plasmid[s])
  --prodigal-tf PRODIGAL_TF
                        Path to existing Prodigal training file to use for CDS prediction
  --translation-table {11,4}
                        Translation table to use: 11/4 (default = 11)
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
                        Number of threads to use (default = number of available CPUs)
  --tmp-dir TMP_DIR     Location for temporary files (default = system dependent auto detection)
  --version             show program's version number and exit
  --citation            Print citation
```

## Annotation workflow

### RNAs

1. tRNA genes: tRNAscan-SE 2.0
2. tmRNA genes: Aragorn
3. rRNA genes: Infernal vs. Rfam rRNA covariance models
4. ncRNA genes: Infernal vs. Rfam ncRNA covariance models
5. ncRNA cis-regulatory regions: Infernal vs. Rfam ncRNA covariance models
6. CRISPR arrays: PILER-CR

Bakta distinguishes ncRNA genes and (regulatory) regions in order to enable the distinct handling thereof during the annotation process, *i.e.* feature overlap detection.

ncRNA gene types:

- sRNA
- antisense
- ribozyme
- antitoxin

ncRNA (cis-regulatory) region types:

- riboswitch
- thermoregulator
- leader
- frameshift element

### Coding sequences

The structural prediction is conducted via Prodigal and complemented by a custom detection of sORF < 30 aa.

To rapidly identify known protein sequences with exact sequence matches and to conduct a comprehensive annotations, Bakta utilizes a compact SQLite database comprising protein sequence digests and pre-annotations for millions of known protein sequences and clusters.

Conceptual terms:

- **UPS**: unique protein sequences identified via length and MD5 digests (100% coverage & 100% sequence identity)
- **IPS**: identical protein sequences comprising seeds of UniProt's UniRef100 protein sequence clusters
- **PSC**: protein sequences clusters comprising seeds of UniProt's UniRef90 protein sequence clusters

**CDS**:

1. Prediction via Prodigal respecting sequences' completeness (distinct prediction for complete replicons and uncompleted contigs)
2. discard spurious CDS via AntiFam
3. Detection of UPSs via MD5 digests and lookup of related IPS and PCS
4. Homology search of remainder via Diamond vs. PSC (coverage=0.8, identity=0.9)
5. Combination of available IPS & PSC information favouring more specific annotations and avoiding redundancy

CDS without IPS or PSC hits as well as those without gene symbols or product descriptions different from `hypothetical` will be marked as `hypothetical`.

**sORFs**:

1. Custom sORF detection & extraction with amino acid lengths < 30 aa
2. Apply strict feature type-dependent overlap filters
3. discard spurious sORF via AntiFam
4. Detection of UPS via MD5 hashes and lookup of related IPS
5. Homology search of remainder via Diamond vs. an sORF subset of PSCs (coverage=0.9, identity=0.9)
6. Exclude sORF without sufficient annotation information

sORF not identified via IPS or PSC will be discarded. Additionally, all sORF without gene symbols or product descriptions different from `hypothetical` will be discarded.
Due due to uncertain nature of sORF prediction, only those identified via IPS / PSC hits exhibiting proper gene symbols or product descriptions different from `hypothetical` will be included in the final annotation.

### Miscellaneous

1. Gaps: in-mem detection & annotation of sequence gaps
2. oriC/oriV/oriT: Blast+ (cov=0.8, id=0.8) vs. [MOB-suite](https://github.com/phac-nml/mob-suite) oriT & [DoriC](http://tubic.org/doric/public/index.php) oriC/oriV sequences. Annotations of ori regions take into account overlapping Blast+ hits and are conducted based on a majority vote heuristic. Region edges are fuzzy - use with caution!

## Database

The Bakta database comprises a set of AA & DNA sequence databases as well as HMM & covariance models.
At its core Bakta utilizes a compact SQLite db storing protein sequence digests, lengths, pre-annotations and dbxrefs of UPS, IPS and PSC from:

- **UPS**: UniParc / UniProtKB (192,795,177)
- **IPS**: UniProt UniRef100 (169,958,214)
- **PSC**: UniProt UniRef90 (77,128,011)

This allows the exact protein sequences identification via MD5 digests & sequence lengths as well as the rapid subsequent lookup of related information. Protein sequence digests are checked for hash collisions while the db creation process.
IPS & PSC have been comprehensively pre-annotated integrating annotations & database *dbxrefs* from:

- NCBI nonredundant proteins ('WP_*' -> 139,330,543)
- NCBI COG db (coverage=80%, identity=90% -> 1,893,080)
- GO terms (via SwissProt IPS/PSC exact matches)
- EC (via SwissProt IPS/PSC exact matches)
- NCBI AMRFinderPlus (IPS exact matches, PSC HMM hits reaching trusted cutoffs)
- ISFinder db (coverage=90%, identify=99% -> 2,981)

Rfam covariance models:

- ncRNA: 750
- ncRNA cis-regulatory regions: 107

To provide FAIR annotations, the database releases are SemVer versioned (w/o patch level), *i.e.* `<major>.<minor>`. For each version we provide a comprehensive log file tracking all imported sequences as well as annotations thereof.
The db schema is represented by the `<major>` digit and automatically checked at runtime by Bakta in order to ensure compatibility. Content updates are tracked by the `<minor>` digit.

All database releases (latest 1.0, 23 Gb zipped, 43 Gb unzipped) are hosted at Zenodo:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4247252.svg)](https://doi.org/10.5281/zenodo.4247252)

## WIP

We're keen to constantly improve and further expand Bakta. If you miss any features, please do not hesitate to ask for it!
There are several features & improvements which we're currently working on or have plans to do so as well as barely ideas:

Workflow:

- analysis & annotation of `hypothetical` CDS [#28](https://github.com/oschwengers/bakta/issues/28)
- detections of pseudo genes [#4](https://github.com/oschwengers/bakta/issues/4)

Annotation:

- (idea) KEGG KofamKOALA annotation [#9](https://github.com/oschwengers/bakta/issues/9)
- (idea) protein classification of IPS/PSC in terms of (`amr`, `virulence`, `mobilization`, `conjugation`, `replication`, *etc*) for subsequent visualization/analysis

Technical:

- 3rd party dependency version checks at runtime [#21](https://github.com/oschwengers/bakta/issues/21)
- CWL description file [#27](https://github.com/oschwengers/bakta/issues/27)
- Expand tests [#29](https://github.com/oschwengers/bakta/issues/29)
- (idea) download/update the database within Bakta

## Citation

A manuscript is not yet in preparation, but might be soon... To temporarily cite our work, please transitionally refer to:
> Schwengers O., Goesmann A. (2020) Bakta: Rapid & standardized annotation of bacterial genomes & plasmids. GitHub https://github.com/oschwengers/bakta

Bakta takes advantage of many publicly available databases. If you find any of the data used within Bakta useful, please also be sure to credit the primary source also:

- UniProt: <https://doi.org/10.1093/nar/gky1049>
- RefSeq: <https://doi.org/10.1093/nar/gkx1068>
- Rfam: <https://doi.org/10.1002/cpbi.51>
- AMRFinder: <https://doi.org/10.1128/AAC.00483-19>
- ISFinder: <https://doi.org/10.1093/nar/gkj014>
- AntiFam: <https://doi.org/10.1093/database/bas003>
- Mob-suite: <https://doi.org/10.1099/mgen.0.000206>
- DoriC: <https://doi.org/10.1093/nar/gky1014>
- COG: <https://doi.org/10.1093/bib/bbx117>

## FAQ

* __Bakta is running too long without CPU load... why?__
Bakta takes advantage of an SQLite DB which results in high storage IO loads. If this DB is stored on a remote / network volume, the lookup of IPS/PSC annotations might take a long time. In these cases, please, consider moving the DB to a local volume/hard drive. Setting POSIX permissions of the db directory to read/access only (`555`) and files to read only (`444`) might also help:

```bash
chmod 444 <db-path>/*
chmod 555 <db-path>
```

## Issues and Feature Requests

Bakta is brand new and like in every software, expect some bugs lurking around. So, if you run into any issues with Bakta, we'd be happy to hear about it.
Therefore, please, execute bakta in verbose mode (`-v`) and do not hesitate to file an issue including as much information as possible:

- a detailed description of the issue
- command line output
- log file (`<prefix>.log`)
- result file (`<prefix>.json`) _if possible_
- a reproducible example of the issue with an input file that you can share _if possible_
