[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/oschwengers/bakta/blob/master/LICENSE)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/bakta.svg)
![GitHub release](https://img.shields.io/github/release/oschwengers/bakta.svg)
[![PyPI](https://img.shields.io/pypi/v/bakta.svg)](https://pypi.org/project/bakta)
![PyPI - Status](https://img.shields.io/pypi/status/bakta.svg)
[![Conda](https://img.shields.io/conda/v/bioconda/bakta.svg)](https://bioconda.github.io/recipes/bakta/README.html)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4247252.svg)](https://doi.org/10.5281/zenodo.4247252)

# Bakta: Rapid & standardized annotation of bacterial genomes & plasmids

Bakta is a tool for the rapid & standardized annotation of bacterial genomes & plasmids. It provides **dbxref**-rich and **sORF**-including annotations in machine-readable `JSON` & bioinformatics standard file formats for automatic downstream analysis.

## Contents

- [Description](#description)
- [Installation](#installation)
- [Examples](#examples)
- [Input/Output](#inputoutput)
- [Usage](#usage)
- [Annotation workflow](#annoation_workflow)
- [Database](#database)
- [Citation](#citation)
- [FAQ](#faq)
- [Issues & Feature Requests](#issues)

## Description

- **Bacteria & plasmids only**
Bakta was designed to annotate bacteria and plasmids, only. This decision by design has been made in order to tweak the annotation process regarding tools, preferences & databases and to streamline further development & maintenance of the software.

- **FAIR annotations**
To provide standardized annotations adhearing to [FAIR](https://www.go-fair.org/fair-principles) principles, Bakta utilizes a comprehensive & versioned custom annotation database based on UniProt's [UniRef100 & UniRef90](https://www.uniprot.org/uniref/) protein clusters (`FAIR` -> [DOI](http://dx.doi.org/10.1038/s41597-019-0180-9)/[DOI](https://doi.org/10.1093/nar/gkaa1100)) enriched with dbxrefs (`GO`, `COG`, `EC`) and annotated by specialized niche databases. For each db version we provide a comprehensive log file of all imported sequences and annotations.

- **Protein sequence identification**
Fostering the FAIR aspect, Bakta identifies identical protein sequences (**IPS**) via `MD5` digests which are annotated with database cross-references (**dbxref**) to RefSeq (`WP_*`), UniRef100 (`UniRef100_*`) and UniParc (`UPI*`).
By doing so, IPS allow the surveillance of distinct gene alleles and streamlining comparative analysis as well as posterior (external) annotations of `putative` & `hypothetical` protein sequences which can be mapped back to existing CDS via these exact & stable identifiers (*E. coli* gene [ymiA](https://www.uniprot.org/uniprot/P0CB62) [...more](https://www.uniprot.org/help/dubious_sequences)). Currently, Bakta identifies ~198 mio, ~185 mio and ~150 mio distinct protein sequences from UniParc, UniRef100 and RefSeq, respectively. Hence, for certain genomes, up to 99 % of all CDS can be identified this way, skipping expensive homology searches.

- **Small proteins / short open reading frames**
Bakta detects and annotates small proteins/short open reading frames (**sORF**) which are not predicted by tools like `Prodigal`.

- **Fast**
Bakta can annotate a typical bacterial genome in 10 &plusmn;5 min on a laptop, plasmids in a couple of seconds/minutes.

- **Expert annotation systems**
To provide high quality annotations for certain proteins of higher interest, *e.g.* AMR & VF genes, Bakta includes & merges different expert annotation systems. Currently, Bakta uses NCBI's AMRFinderPlus for AMR gene annotations as well as an generalized protein sequence expert system with distinct coverage, identity and priority values for each sequence, currenlty comprising the [VFDB](http://www.mgc.ac.cn/VFs/main.htm) as well as NCBI's [BlastRules](https://ftp.ncbi.nih.gov/pub/blastrules/).

- **Comprehensive workflow**
Bakta annotates all standard/mandatory feature types (tRNA, tmRNA, rRNA, ncRNA genes, CRISPR, CDS) as well as ncRNA cis-regulatory regions, oriC/oriV/oriT and assembly gaps.

- **Reasoning**
By annotating bacterial genomes in a standardized, taxon-independent, high-throughput and local manner, Bakta targets the niche between fully-featured but computationally-demanding pipelines like [PGAP](https://github.com/ncbi/pgap) and rapid highly-customizable offline tools like [Prokka](https://github.com/tseemann/prokka). Indeed, Bakta is heavily inspired by Prokka (kudos to [Torsten Seemann](https://github.com/tseemann)) and many command line options are compatible for the sake of interoperability and user convenience. Hence, if Bakta does not fit your needs, please try Prokka.

## Installation

Bakta can be installed via BioConda, Docker, Singularity and Pip.
However, we encourage to use [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) or [Docker](https://www.docker.com/get-started)/[Singularity](https://sylabs.io/singularity) to automatically install all required 3rd party dependencies. In all cases a mandatory [database](#database_download) must be downloaded.

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
$ singularity build bakta.sif docker://oschwengers/bakta:latest
$ singularity run bakta.sif --help
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
- Diamond (2.0.6) <https://doi.org/10.1038/nmeth.3176> <https://github.com/bbuchfink/diamond>
- Blast+ (2.7.1) <https://www.ncbi.nlm.nih.gov/pubmed/2231712> <https://blast.ncbi.nlm.nih.gov>
- AMRFinderPlus (3.10.1) <https://github.com/ncbi/amr>

On Ubuntu/Debian/Mint you can install these via:

```bash
$ sudo apt install trnascan-se aragorn infernal pilercr prodigal hmmer diamond-aligner ncbi-blast+
```

Tested with Ubuntu 20.04 - some older distributions might provide outdated versions, *e.g.* trnascan-se in Ubuntu 18.04. In these cases dependencies must be installed manually.

### Database download

Bakta requires a mandatory database which is publicly hosted at Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4247252.svg)](https://doi.org/10.5281/zenodo.4247252)
Further information is provided in the [database](#database) section below.

List available DB versions:

```bash
$ bakta_db list
...
```

Download the most recent compatible database version:

```bash
$ bakta_db download --output <output-path>
```

Update an existing database:

```bash
$ bakta_db update --db <existing-db-path> [--tmp-dir <tmp-directory>]
```

or download it manually:

```bash
$ wget https://zenodo.org/record/4337350/files/db.tar.gz
$ tar -xzf db.tar.gz
$ rm db.tar.gz
```

The db path can be provided either via parameter (`--db`) or environment variable (`BAKTA_DB`):

```bash
$ bakta --db <db-path> genome.fasta

$ export BAKTA_DB=<db-path>
$ bakta genome.fasta
```

For system-wide setups, the database can also be copied to the Bakta base directory:

```bash
$ cp -r db/ <bakta-installation-dir>
```

## Examples

Simple:

```bash
$ bakta --db <db-path> genome.fasta
```

Expert: verbose output writing results to *results* directory with *ecoli123* file `prefix` and *eco634* `locus tag` using an existing prodigal training file, using additional replicon information and 8 threads:

```bash
$ bakta --db <db-path> --verbose --output results/ --prefix ecoli123 --locus-tag eco634 --prodigal-tf eco.tf --replicons replicon.tsv --threads 8 genome.fasta
```

## Input/Output

### Input

Bakta accepts bacterial genomes and plasmids (complete / draft assemblies) in (zipped) fasta format. For a full description of how further genome information can be provided and workflow customizations can be set, please have a look at the [Usage](#usage) section.

Replicon meta data table:

To fine-tune the very details of each sequence in the input fasta file, Bakta accepts a replicon meta data table provided in `csv`/`tsv` file format: `--replicons <file.tsv>`. Thus, complete replicons within partially completed draft assemblies can be marked & handled as such, *e.g.* detection & annotation of features spanning sequence edges.

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

Bakta provides detailed information on each annotated feature in a standardized machine-readable JSON file `<prefix>.json`:

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

- `<prefix>.tsv`: annotations as simple human readble tab separated values
- `<prefix>.gff3`: annotations & sequences in GFF3 format
- `<prefix>.gbff`: annotations & sequences in (multi) GenBank format
- `<prefix>.embl`: annotations & sequences in (multi) EMBL format
- `<prefix>.fna`: replicon/contig DNA sequences as FASTA
- `<prefix>.faa`: CDS/sORF amino acid sequences as FASTA
- `<prefix>.hypotheticals.tsv`: further information on hypothetical protein CDS as simple human readble tab separated values
- `<prefix>.hypotheticals.faa`: hypothetical protein CDS amino acid sequences as FASTA

The `<prefix>` can be set via `--prefix <prefix>`. If no prefix is set, Bakta uses the input file prefix.

## Usage

Usage:

```bash
usage: bakta [--db DB] [--min-contig-length MIN_CONTIG_LENGTH] [--prefix PREFIX] [--output OUTPUT] [--genus GENUS] [--species SPECIES] [--strain STRAIN] [--plasmid PLASMID] [--complete] [--prodigal-tf PRODIGAL_TF] [--translation-table {11,4}] [--gram {+,-,?}] [--locus LOCUS]
             [--locus-tag LOCUS_TAG] [--keep-contig-headers] [--replicons REPLICONS] [--skip-trna] [--skip-tmrna] [--skip-rrna] [--skip-ncrna] [--skip-ncrna-region] [--skip-crispr] [--skip-cds] [--skip-sorf] [--skip-gap] [--skip-ori] [--help] [--verbose] [--threads THREADS]
             [--tmp-dir TMP_DIR] [--version]
             <genome>

Rapid & standardized annotation of bacterial genomes & plasmids.

positional arguments:
  <genome>              Genome sequences in (zipped) fasta format

Input / Output:
  --db DB, -d DB        Database path (default = <bakta_path>/db). Can also be provided as BAKTA_DB environment variable.
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
                        Translation table: 11/4 (default = 11)
  --gram {+,-,?}        Gram type: +/-/? (default = '?')
  --locus LOCUS         Locus prefix (default = 'contig')
  --locus-tag LOCUS_TAG
                        Locus tag prefix (default = autogenerated)
  --keep-contig-headers
                        Keep original contig headers
  --replicons REPLICONS, -r REPLICONS
                        Replicon information table (tsv)

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
4. Homology search of remainder via Diamond vs. PSC (query/subject coverage=0.8, identity=0.9)
5. Execution of expert systems:
  - AMR: AMRFinderPlus
  - Alignments: NCBI BlastRules, VFDB
6. Combination of available IPS, PSC & expert system information favouring more specific annotations and avoiding redundancy

CDS without IPS or PSC hits as well as those without gene symbols or product descriptions different from `hypothetical` will be marked as `hypothetical`.

Such hypothetical CDS are further analyzed:

1. Detection of Pfam domains, repeats & motifs
2. Calculation of protein sequence statistics, *i.e.* molecular weight, isoelectric point

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

- **UPS**: UniParc / UniProtKB (198,764,035)
- **IPS**: UniProt UniRef100 (185,077,759)
- **PSC**: UniProt UniRef90 (83,486,930)

This allows the exact protein sequences identification via MD5 digests & sequence lengths as well as the rapid subsequent lookup of related information. Protein sequence digests are checked for hash collisions while the db creation process.
IPS & PSC have been comprehensively pre-annotated integrating annotations & database *dbxrefs* from:

- NCBI nonredundant proteins (IPS: 150,463,165)
- NCBI COG db (PSC: 3,353,131)
- SwissProt EC/GO terms (PSC: 335,256)
- NCBI AMRFinderPlus (IPS: 6,255, PSC: 43,413)
- ISFinder db (IPS: 27,841, PSC: 10,161)

To provide high quality annotations for distinct protein sequences of high importance (AMR, VF, *etc*) which cannot sufficiently be covered by the IPS/PSC approach, Bakta provides additional expert systems. For instance, AMR genes, are annotated via NCBI's AMRFinderPlus.
An expandable alignment-based expert system supports the incorporation of high quality annotations from multiple sources. This currenlty comprises NCBI's BlastRules as well as VFDB and will be complemented with more expert annotation sources over time. Internally, this expert system is based on a Diamond DB comprising the following information in a standardized format:

- source: *e.g.* BlastRules
- rank: a precedence rank
- min identity
- min query coverage
- min model coverage
- gene lable
- product description
- dbxrefs

Rfam covariance models:

- ncRNA: 750
- ncRNA cis-regulatory regions: 107

To provide FAIR annotations, the database releases are SemVer versioned (w/o patch level), *i.e.* `<major>.<minor>`. For each version we provide a comprehensive log file tracking all imported sequences as well as annotations thereof. The db schema is represented by the `<major>` digit and automatically checked at runtime by Bakta in order to ensure compatibility. Content updates are tracked by the `<minor>` digit.

All database releases (latest 1.0, 25 Gb zipped, 48 Gb unzipped) are hosted at Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo..svg)](https://doi.org/10.5281/zenodo.)

## Citation

A manuscript is in preparation. To temporarily cite our work, please transitionally refer to:
> Schwengers O., Goesmann A. (2021) Bakta: Rapid & standardized annotation of bacterial genomes & plasmids. GitHub https://github.com/oschwengers/bakta

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
- VFDB: <https://doi.org/10.1093/nar/gky1080>

## FAQ

* __Nice, but I'm mising XYZ...__
Bakta is quite new and we're keen to constantly improve it and further expand its feature set. In case there's anything missing, please do not hesitate to open an issue and ask for it!

* __Bakta is running too long without CPU load... why?__
Bakta takes advantage of an SQLite DB which results in high storage IO loads. If this DB is stored on a remote / network volume, the lookup of IPS/PSC annotations might take a long time. In these cases, please, consider moving the DB to a local volume or hard drive. Setting POSIX permissions of the db directory to read/access only (`555`) and files to read only (`444`) might also help:

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
