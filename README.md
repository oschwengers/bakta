[![DOI:10.1099/mgen.0.000685](https://zenodo.org/badge/DOI/10.1099/mgen.0.000685.svg)](https://doi.org/10.1099/mgen.0.000685)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4247252.svg)](https://doi.org/10.5281/zenodo.4247252)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/oschwengers/bakta/blob/master/LICENSE)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/bakta.svg)
![PyPI - Status](https://img.shields.io/pypi/status/bakta.svg)
![GitHub release](https://img.shields.io/github/release/oschwengers/bakta.svg)

[![PyPI](https://img.shields.io/pypi/v/bakta.svg)](https://pypi.org/project/bakta)
[![Conda](https://img.shields.io/conda/vn/bioconda/bakta.svg)](https://bioconda.github.io/recipes/bakta/README.html)
[![Docker Image Version](https://img.shields.io/docker/v/oschwengers/bakta?sort=semver&label=docker)](https://hub.docker.com/r/oschwengers/bakta)
[![Spack](https://img.shields.io/spack/v/py-bakta)](https://packages.spack.io/package.html?name=py-bakta)
[![Galaxy Toolshed - Tool Version](https://img.shields.io/galaxytoolshed/v/bakta/iuc/bakta?label=usegalaxy.eu)](https://usegalaxy.eu/root?tool_id=bakta)
[![Static Badge](https://img.shields.io/badge/bio.tools-v1.9.4-blue?link=https%3A%2F%2Fbio.tools%2Fbakta)](https://bio.tools/bakta)

# Bakta: rapid & standardized annotation of bacterial genomes, MAGs & plasmids

Bakta is a tool for the rapid & standardized annotation of bacterial genomes and plasmids from both isolates and MAGs. It provides **dbxref**-rich, **sORF**-including and taxon-independent annotations in machine-readable `JSON` & bioinformatics standard file formats for automated downstream analysis.

## Contents

- [Description](#description)
- [Installation](#installation)
- [Examples](#examples)
- [Input & Output](#input-and-output)
- [Usage](#usage)
- [Annotation Workflow](#annotation-workflow)
- [Database](#database)
- [Genome Submission](#genome-submission)
- [Protein bulk annotation](#protein-bulk-annotation)
- [Genome plots](#genome-plots)
- [Auxiliary scripts](#auxiliary-scripts)
- [Web version](#web-version)
- [Citation](#citation)
- [FAQ](#faq)
- [Issues & Feature Requests](#issues-and-feature-requests)

## Description

- **Comprehensive & taxonomy-independent database**
Bakta provides a large and taxonomy-independent database using UniProt's entire [UniRef](https://www.uniprot.org/uniref/) protein sequence cluster universe. Thus, it achieves favourable annotations in terms of sensitivity and specificity along the broad continuum ranging from well-studied species to unknown genomes from MAGs.

- **Protein sequence identification**
Bakta exactly identifies known identical protein sequences (**IPS**) from RefSeq and UniProt allowing the fine-grained annotation of gene alleles (`AMR`) or closely related but distinct protein families. This is achieved via an alignment-free sequence identification (**AFSI**) approach using full-length `MD5` protein sequence hash digests.

- **Fast**
This AFSI approach substantially accellerates the annotation process by avoiding computationally expensive homology searches for identified genes. Thus, Bakta can annotate a typical bacterial genome in 10 &plusmn;5 min on a laptop, plasmids in a couple of seconds/minutes.

- **Database cross-references**
Fostering the [FAIR](https://www.go-fair.org/fair-principles) principles, Bakta exploits its AFSI approach to annotate CDS with database cross-references (**dbxref**) to RefSeq (`WP_*`), UniRef100 (`UniRef100_*`) and UniParc (`UPI*`). By doing so, IPS allow the surveillance of distinct gene alleles and streamlining comparative analysis as well as posterior (external) annotations of `putative` & `hypothetical` protein sequences which can be mapped back to existing CDS via these exact & stable identifiers (*E. coli* gene [ymiA](https://www.uniprot.org/uniprot/P0CB62) [...more](https://www.uniprot.org/help/dubious_sequences)). Currently, Bakta identifies ~350 mio, ~330 mio and ~290 mio distinct protein sequences from UniParc, UniRef100 and RefSeq, respectively. Hence, for certain genomes, up to 99 % of all CDS can be identified this way, skipping computationally expensive sequence alignments.

- **FAIR annotations**
To provide standardized annotations adhearing to FAIR principles, Bakta utilizes a versioned custom annotation database comprising UniProt's [UniRef100 & UniRef90](https://www.uniprot.org/uniref/) protein clusters (FAIR -> [DOI](http://dx.doi.org/10.1038/s41597-019-0180-9)/[DOI](https://doi.org/10.1093/nar/gkaa1100)) enriched with dbxrefs (`GO`, `COG`, `EC`) and annotated by specialized niche databases. For each DB version we provide a comprehensive log file of all imported sequences and annotations.

- **Small proteins / short open reading frames**
Bakta detects and annotates small proteins/short open reading frames (**sORF**) which are not predicted by tools like `Prodigal`.

- **Expert annotation systems**
To provide high quality annotations for certain proteins of higher interest, *e.g.* AMR & VF genes, Bakta includes & merges different expert annotation systems. Currently, Bakta uses NCBI's AMRFinderPlus for AMR gene annotations as well as an generalized protein sequence expert system with distinct coverage, identity and priority values for each sequence, currenlty comprising the [VFDB](http://www.mgc.ac.cn/VFs/main.htm) as well as NCBI's [BlastRules](https://ftp.ncbi.nih.gov/pub/blastrules/).

- **Comprehensive workflow**
Bakta annotates ncRNA cis-regulatory regions, oriC/oriV/oriT and assembly gaps as well as standard feature types: tRNA, tmRNA, rRNA, ncRNA genes, CRISPR, CDS and pseudogenes.

- **GFF3 & INSDC conform annotations**
Bakta writes GFF3 and INSDC-compliant (Genbank & EMBL) annotation files ready for submission (checked via [GenomeTools GFF3Validator](http://genometools.org/cgi-bin/gff3validator.cgi), [table2asn_GFF](https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/#run) and [ENA Webin-CLI](https://github.com/enasequence/webin-cli) for GFF3 and EMBL file formats, respectively for representative genomes of all ESKAPE species).

- **Bacteria & plasmids only**
Bakta was designed to annotate bacteria (isolates & MAGs) and plasmids, only. This decision by design has been made in order to tweak the annotation process regarding tools, preferences & databases and to streamline further development & maintenance of the software.

- **Reasoning**
By annotating bacterial genomes in a standardized, taxonomy-independent, high-throughput and local manner, Bakta aims at a well-balanced tradeoff between fully featured but computationally demanding pipelines like [PGAP](https://github.com/ncbi/pgap) and rapid highly customizable offline tools like [Prokka](https://github.com/tseemann/prokka). Indeed, Bakta is heavily inspired by Prokka (kudos to [Torsten Seemann](https://github.com/tseemann)) and many command line options are compatible for the sake of interoperability and user convenience. Hence, if Bakta does not fit your needs, please consider trying Prokka.

## Installation

Bakta can be installed via BioConda, Docker, Singularity and Pip. However, we encourage to use [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) or [Docker](https://www.docker.com/get-started)/[Singularity](https://sylabs.io/singularity) to automatically install all required 3rd party dependencies.

In all cases a mandatory [database](#database-download) must be downloaded.

### BioConda

```bash
conda install -c conda-forge -c bioconda bakta
```

### Podman (Docker)

We maintain a Docker image `oschwengers/bakta` providing an entrypoint, so that containers can be used like an executable:

```bash
podman pull oschwengers/bakta
podman run oschwengers/bakta --help
```

Installation instructions and get-started guides: Podman [docs](https://podman.io/docs). For further convenience, we provide a shell script (`bakta-podman.sh`) handling Podman related parameters (volume mounting, user IDs, etc):

```bash
bakta-podman.sh --db <db-path> --output <output-path> <input>
```

For experienced users and full functionality (`bakta_db` & `bakta_proteins`), an image without entrypoint might be a better option. For these cases, please use one of the [Biocontainer](https://quay.io/repository/biocontainers/bakta?tab=tags) images:

```bash
export CONTAINER="quay.io/biocontainers/bakta:1.8.2--pyhdfd78af_0"
podman run -it --rm $CONTAINER bakta --help
podman run -it --rm $CONTAINER bakta_db --help
```

### Pip

```bash
python3 -m pip install --user bakta
```

Bakta requires the following 3rd party software tools which must be installed and executable to use the full set of features:

- tRNAscan-SE (2.0.11) <https://doi.org/10.1101/614032> <http://lowelab.ucsc.edu/tRNAscan-SE>
- Aragorn (1.2.41) <http://dx.doi.org/10.1093/nar/gkh152> <http://130.235.244.92/ARAGORN>
- INFERNAL (1.1.4) <https://dx.doi.org/10.1093%2Fbioinformatics%2Fbtt509> <http://eddylab.org/infernal>
- PILER-CR (1.06) <https://doi.org/10.1186/1471-2105-8-18> <http://www.drive5.com/pilercr>
- Pyrodigal (3.5.0) <https://doi.org/10.21105/joss.04296> <https://github.com/althonos/pyrodigal>
- PyHMMER (0.10.15) <https://doi.org/10.21105/joss.04296> <https://github.com/althonos/pyhmmer>
- Diamond (2.1.10) <https://doi.org/10.1038/nmeth.3176> <https://github.com/bbuchfink/diamond>
- Blast+ (2.14.0) <https://www.ncbi.nlm.nih.gov/pubmed/2231712> <https://blast.ncbi.nlm.nih.gov>
- AMRFinderPlus (4.0.3) <https://github.com/ncbi/amr>
- pyCirclize (1.7.0) https://github.com/moshi4/pyCirclize

### Database download

Bakta requires a mandatory database which is publicly hosted at Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4247252.svg)](https://doi.org/10.5281/zenodo.4247252)
We provide 2 types: `full` and `light`. To get best annotation results and to use all features, we recommend using the `full` (default). If you seek for maximum runtime performance or if download time/storage requirements are an issue, please try the `light` version. Further information is provided in the [database](#database) section below.

List available DB versions (available as either `full` or `light`):

```bash
bakta_db list
...
```

To download the most recent compatible database version we recommend to use the internal database download & setup tool:

```bash
bakta_db download --output <output-path> --type [light|full]
```

Of course, the database can also be downloaded and installed manually:

```bash
wget https://zenodo.org/record/14916843/files/db-light.tar.xz
bakta_db install -i db-light.tar.xz
```

If required, or desired, the AMRFinderPlus DB can also be updated manually:

```bash
amrfinder_update --force_update --database db-light/amrfinderplus-db/
```

If you're using bakta on Docker:

```bash
docker run -v /path/to/desired-db-path:/db --entrypoint /bin/bash oschwengers/bakta:latest -c "bakta_db download --output /db --type [light|full]"
```

As an additional data repository backup, we provide the most recent database version via our institute servers: [full](https://s3.computational.bio.uni-giessen.de/bakta-db/db-v6.0.tar.xz), [light](https://s3.computational.bio.uni-giessen.de/bakta-db/db-light-v6.0.tar.xz). However, the bandwith is limited. Hence, please use it with caution and only if Zenodo might be temporarily unreachable or slow.

Update an existing database:

```bash
bakta_db update --db <existing-db-path> [--tmp-dir <tmp-directory>]
```

Update using Docker:

```bash
docker run -v /path/to/desired-db-path:/db --entrypoint /bin/bash oschwengers/bakta:latest -c "bakta_db update --db /db/db-[light|full]"
```

The database path can be provided either via parameter (`--db`) or environment variable (`BAKTA_DB`):

```bash
bakta --db <db-path> genome.fasta

export BAKTA_DB=<db-path>
bakta genome.fasta
```

For system-wide setups, the database can also be copied to the Bakta base directory:

```bash
cp -r db/ <bakta-installation-dir>
```

As Bakta takes advantage of AMRFinderPlus for the annotation of AMR genes, AMRFinder is required to setup its own internal databases in a `<amrfinderplus-db>` subfolder within the Bakta database `<db-path>`, once via `amrfinder_update --force_update --database <db-path>/amrfinderplus-db/`. To ease this process we recommend to use Bakta's internal download procedure.

## Examples

Simple:

```bash
bakta --db <db-path> genome.fasta
```

Expert: verbose output writing results to *results* directory with *ecoli123* file `prefix` and *eco634* `locus tag` using an existing prodigal training file, using additional replicon information and 8 threads:

```bash
bakta --db <db-path> --verbose --output results/ --prefix ecoli123 --locus-tag eco634 --prodigal-tf eco.tf --replicons replicon.tsv --threads 8 genome.fasta
```

## Input and Output

### Input

Bakta accepts bacterial genomes and plasmids (complete / draft assemblies) in (zipped) fasta format. For a full description of how further genome information can be provided and workflow customizations can be set, please have a look at the [Usage](#usage) section or this [manual](https://bakta.readthedocs.io/).

#### Replicon meta data table

To fine-tune the very details of each sequence in the input fasta file, Bakta accepts a replicon meta data table provided in `csv` or `tsv` file format: `--replicons <file.tsv>`. Thus, complete replicons within partially completed draft assemblies can be marked & handled as such, *e.g.* detection & annotation of features spanning sequence edges.

Table format:

original sequence id  |  new sequence id  |  type  |  topology  |  name
----|----------------|----------------|----------------|----------------
`old id` | `new id`, `<empty>` | `chromosome`, `plasmid`, `contig`, `<empty>` | `circular`, `linear`, `<empty>` | `name`, `<empty>`

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
NODE_4 | special-contig-name-xyz |  `-` | `-` | `-`
NODE_5 | `` |  `-` | `-` | `-`

#### User-provided regions

Bakta accepts pre-annotated (*a priori*), user-provided feature regions via `--regions` in either GFF3 or GenBank format. These regions supersede all de novo-predicted regions, but are equally subject to the internal functional annotation process. Currently, only `CDS` are supported. A maximum overlap with *de novo*-predicted CDS of 30 bp is allowed. If you would like to provide custom functional annotations, you can provide these via `--proteins` which is described in the following section.

#### User-provided protein sequences

Bakta accepts user-provided trusted protein sequences via `--proteins` in either GenBank (CDS features) or Fasta format which are used in the functional annotation process. Using the Fasta format, each reference sequence can be provided in a short or long format:

```bash
# short:
>id gene~~~product~~~dbxrefs
MAQ...

# long:
>id min_identity~~~min_query_cov~~~min_subject_cov~~~gene~~~product~~~dbxrefs
MAQ...
```

Allowed values:

field  |  value(s)  |  example
----|----------------|----------------
min_identity | `int`, `float` | 80, 90.3
min_query_cov | `int`, `float` | 80, 90.3
min_subject_cov | `int`, `float` | 80, 90.3
gene | `<empty>`, `string` | msp
product | `string` | my special protein
dbxrefs | `<empty>`, `db:id`, `,` separated list  | `VFDB:VF0511`

Protein sequences provided in short Fasta or GenBank format are searched with default thresholds of 90%, 80% and 80% for minimal identity, query and subject coverage, respectively.

#### User-provided HMMs

Bakta accepts user-provided trusted HMMs via `--hmms` in HMMER's text format. If set, Bakta will adhere to the *trusted cutoff* specified in the HMM header. In addition, a max. evalue threshold of 1e-6 is applied. By default, Bakta uses the HMM description line as a product description. Further information can be provided via the HMM description line using the *short* format as explained above in the [User-provided protein sequences](####user-provided-protein-sequences) section.

```bash
# default
HMMER3/f [3.1b2 | February 2015]
NAME  id
ACC   id
DESC  product
LENG  435
TC    600 600

# short
NAME  id
ACC   id
DESC  gene~~~product~~~dbxrefs
LENG  435
TC    600 600
```

### Output

Annotation results are provided in standard bioinformatics file formats:

- `<prefix>.tsv`: annotations as simple human readble TSV
- `<prefix>.gff3`: annotations & sequences in GFF3 format
- `<prefix>.gbff`: annotations & sequences in (multi) GenBank format
- `<prefix>.embl`: annotations & sequences in (multi) EMBL format
- `<prefix>.fna`: replicon/contig DNA sequences as FASTA
- `<prefix>.ffn`: feature nucleotide sequences as FASTA
- `<prefix>.faa`: CDS/sORF amino acid sequences as FASTA
- `<prefix>.inference.tsv`: inference metrics (score, evalue, coverage, identity) for annotated accessions as TSV
- `<prefix>.hypotheticals.tsv`: further information on hypothetical protein CDS as simple human readble tab separated values
- `<prefix>.hypotheticals.faa`: hypothetical protein CDS amino acid sequences as FASTA
- `<prefix>.txt`: summary as TXT
- `<prefix>.png`: circular genome annotation plot as PNG
- `<prefix>.svg`: circular genome annotation plot as SVG
- `<prefix>.json`: all (internal) annotation & sequence information as JSON

The `<prefix>` can be set via `--prefix <prefix>`. If no prefix is set, Bakta uses the input file prefix.

Of note, Bakta provides all detailed (internal) information on each annotated feature in a standardized machine-readable JSON file `<prefix>.json`:

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
            "contig": "contig_1",
            "start": 971,
            "stop": 1351,
            "strand": "-",
            "gene": "lsoB",
            "product": "type II toxin-antitoxin system antitoxin LsoB",
            ...
        },
        ...
    ],
    "sequences": [
        {
            "id": "c1",
            "description": "[organism=Escherichia coli] [completeness=complete] [topology=circular]",
            "nt": "AGCTTT...",
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

Bakta provides a helper function to create above mentioned output files from the (GNU-zipped) *JSON* result file, thus helping potential long-term or large-scale annotation projects to reduce overall storage requirements.

```bash
bakta_io --output <output-path> --prefix <prefix> result.json.gz

bakta_io --help
```

Exemplary annotation result files for several genomes (mostly ESKAPE species) are hosted at Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4770026.svg)](https://doi.org/10.5281/zenodo.4770026)

## Usage

```bash
usage: bakta [--db DB] [--min-contig-length MIN_CONTIG_LENGTH] [--prefix PREFIX] [--output OUTPUT] [--force]
             [--genus GENUS] [--species SPECIES] [--strain STRAIN] [--plasmid PLASMID]
             [--complete] [--prodigal-tf PRODIGAL_TF] [--translation-table {11,4,25}] [--gram {+,-,?}]
             [--locus LOCUS] [--locus-tag LOCUS_TAG] [--locus-tag-increment {1,5,10}] [--keep-contig-headers] [--compliant]
             [--replicons REPLICONS] [--regions REGIONS] [--proteins PROTEINS] [--hmms HMMS] [--meta]
             [--skip-trna] [--skip-tmrna] [--skip-rrna] [--skip-ncrna] [--skip-ncrna-region]
             [--skip-crispr] [--skip-cds] [--skip-pseudo] [--skip-sorf] [--skip-gap] [--skip-ori] [--skip-filter] [--skip-plot]
             [--help] [--verbose] [--debug] [--threads THREADS] [--tmp-dir TMP_DIR] [--version]
             <genome>

Rapid & standardized annotation of bacterial genomes, MAGs & plasmids

positional arguments:
  <genome>              Genome sequences in (zipped) fasta format

Input / Output:
  --db DB, -d DB        Database path (default = <bakta_path>/db). Can also be provided as BAKTA_DB environment variable.
  --min-contig-length MIN_CONTIG_LENGTH, -m MIN_CONTIG_LENGTH
                        Minimum contig/sequence size (default = 1; 200 in compliant mode)
  --prefix PREFIX, -p PREFIX
                        Prefix for output files
  --output OUTPUT, -o OUTPUT
                        Output directory (default = current working directory)
  --force, -f           Force overwriting existing output folder (except for current working directory)

Organism:
  --genus GENUS         Genus name
  --species SPECIES     Species name
  --strain STRAIN       Strain name
  --plasmid PLASMID     Plasmid name

Annotation:
  --complete            All sequences are complete replicons (chromosome/plasmid[s])
  --prodigal-tf PRODIGAL_TF
                        Path to existing Prodigal training file to use for CDS prediction
  --translation-table {11,4,25}
                        Translation table: 11/4/25 (default = 11)
  --gram {+,-,?}        Gram type for signal peptide predictions: +/-/? (default = ?)
  --locus LOCUS         Locus prefix (default = 'contig')
  --locus-tag LOCUS_TAG
                        Locus tag prefix (default = autogenerated)
  --locus-tag-increment {1,5,10}
                        Locus tag increment: 1/5/10 (default = 1)

  --keep-contig-headers
                        Keep original contig/sequence headers
  --compliant           Force Genbank/ENA/DDJB compliance
  --replicons REPLICONS, -r REPLICONS
                        Replicon information table (tsv/csv)
  --regions REGIONS     Path to pre-annotated regions in GFF3 or Genbank format (regions only, no functional annotations).
  --proteins PROTEINS   Fasta file of trusted protein sequences for CDS annotation
  --hmms HMMS           HMM file of trusted hidden markov models in HMMER format for CDS annotation
  --meta                Run in metagenome mode. This only affects CDS prediction.

Workflow:
  --skip-trna           Skip tRNA detection & annotation
  --skip-tmrna          Skip tmRNA detection & annotation
  --skip-rrna           Skip rRNA detection & annotation
  --skip-ncrna          Skip ncRNA detection & annotation
  --skip-ncrna-region   Skip ncRNA region detection & annotation
  --skip-crispr         Skip CRISPR array detection & annotation
  --skip-cds            Skip CDS detection & annotation
  --skip-pseudo         Skip pseudogene detection & annotation
  --skip-sorf           Skip sORF detection & annotation
  --skip-gap            Skip gap detection & annotation
  --skip-ori            Skip oriC/oriT detection & annotation
  --skip-filter         Skip feature overlap filters
  --skip-plot           Skip generation of circular genome plots

General:
  --help, -h            Show this help message and exit
  --verbose, -v         Print verbose information
  --debug               Run Bakta in debug mode. Temp data will not be removed.
  --threads THREADS, -t THREADS
                        Number of threads to use (default = number of available CPUs)
  --tmp-dir TMP_DIR     Location for temporary files (default = system dependent auto detection)
  --version             show program's version number and exit
```

## Annotation Workflow

### RNAs

1. tRNA genes: tRNAscan-SE 2.0
2. tmRNA genes: Aragorn
3. rRNA genes: Infernal vs. Rfam rRNA covariance models
4. ncRNA genes: Infernal vs. Rfam ncRNA covariance models
5. ncRNA cis-regulatory regions: Infernal vs. Rfam ncRNA covariance models
6. CRISPR arrays: PILER-CR

Bakta distinguishes ncRNA genes and (cis-regulatory) regions in order to enable the distinct handling thereof during the annotation process, *i.e.* feature overlap detection.

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

The structural prediction is conducted via Pyrodigal and complemented by a custom detection of sORF < 30 aa. In addition, superseding regions of pre-predicted CDS can be provided via `--regions`.

To rapidly identify known protein sequences with exact sequence matches and to conduct a comprehensive annotations, Bakta utilizes a compact read-only SQLite database comprising protein sequence digests and pre-assigned annotations for millions of known protein sequences and clusters.

Conceptual terms:

- **UPS**: unique protein sequences identified via length and MD5 hash digests (100% coverage & 100% sequence identity)
- **IPS**: identical protein sequences comprising seeds of UniProt's UniRef100 protein sequence clusters
- **PSC**: protein sequences clusters comprising seeds of UniProt's UniRef90 protein sequence clusters
- **PSCC**: protein sequences clusters of clusters comprising annotations of UniProt's UniRef50 protein sequence clusters

**CDS**:

1. De novo-prediction via Pyrodigal respecting sequences' completeness (distinct prediction for complete replicons and uncompleted contigs)
2. Discard spurious CDS via AntiFam
3. Detect translational exceptions (selenocysteines)
4. Import of superseding user-provided CDS regions (optional)
5. Detection of UPSs via MD5 digests and lookup of related IPS and PCS
6. Sequence alignments of remainder via Diamond vs. PSC (query/subject coverage=0.8, identity=0.5)
7. Assignment to UniRef90 or UniRef50 clusters if alignment hits achieve identities larger than 0.9 or 0.5, respectively
8. Execution of expert systems:
   - AMR: AMRFinderPlus
   - Expert proteins: NCBI BlastRules, VFDB
   - User proteins (optionally via `--proteins <Fasta/GenBank>`)
9. Prediction of signal peptides (optionally via `--gram <+/->`)
10. Detection of pseudogenes:
   1. Search for reference PCSs using `hypothetical` CDS as seed sequences
   2. Translated alignment (blastx) of reference PCSs against up-/downstream-elongated CDS regions
   3. Analysis of translated alignments and detection of pseudogenization causes & effects
11. Combination of IPS, PSC, PSCC and expert system information favouring more specific annotations and avoiding redundancy

CDS without IPS or PSC hits as well as those without gene symbols or product descriptions different from `hypothetical` will be marked as `hypothetical`.

Such hypothetical CDS are further analyzed:

1. Detection of Pfam domains, repeats & motifs
2. Calculation of protein sequence statistics, *i.e.* molecular weight, isoelectric point

**sORFs**:

1. Custom sORF detection & extraction with amino acid lengths < 30 aa
2. Apply strict feature type-dependent overlap filters
3. discard spurious sORF via AntiFam
4. Detection of UPS via MD5 hashes and lookup of related IPS
5. Sequence alignments of remainder via Diamond vs. an sORF subset of PSCs (coverage=0.9, identity=0.9)
6. Exclude sORF without sufficient annotation information
7. Prediction of signal peptides (optionally via `--gram <+/->`)

sORF not identified via IPS or PSC will be discarded. Additionally, all sORF without gene symbols or product descriptions different from `hypothetical` will be discarded.
Due due to uncertain nature of sORF prediction, only those identified via IPS / PSC hits exhibiting proper gene symbols or product descriptions different from `hypothetical` will be included in the final annotation.

### Miscellaneous

1. Gaps: in-mem detection & annotation of sequence gaps
2. oriC/oriV/oriT: Blast+ (cov=0.8, id=0.8) vs. [MOB-suite](https://github.com/phac-nml/mob-suite) oriT & [DoriC](http://tubic.org/doric/public/index.php) oriC/oriV sequences. Annotations of ori regions take into account overlapping Blast+ hits and are conducted based on a majority vote heuristic. Region edges are fuzzy - use with caution!

## Database

The Bakta database comprises a set of AA & DNA sequence databases as well as HMM & covariance models.
At its core Bakta utilizes a compact read-only SQLite DB storing protein sequence digests, lengths, pre-assigned annotations and dbxrefs of UPS, IPS and PSC from:

- **UPS**: UniParc / UniProtKB (350,631,327)
- **IPS**: UniProt UniRef100 (330,865,009)
- **PSC**: UniProt UniRef90 (135,274,518)
- **PSCC**: UniProt UniRef50 (37,008,138)

This allows the exact protein sequences identification via MD5 digests & sequence lengths as well as the rapid subsequent lookup of related information. Protein sequence digests are checked for hash collisions while the DB creation process. IPS & PSC have been comprehensively pre-annotated integrating annotations & database *dbxrefs* from:

- NCBI nonredundant proteins (UPS: 290,693,966)
- NCBI COG DB (PSC: 3,513,643)
- KEGG Kofams (PSC: 24,267,514)
- SwissProt EC/GO terms (PSC: 337,264)
- NCBI NCBIfams (PSC: 21,758,901)
- PHROG (PSC: 11,717)
- NCBI AMRFinderPlus (IPS: 8,382)
- ISFinder DB (IPS: 155,449, PSC: 14,481)
- Pfam families (PSC: 659,781)

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

- ncRNA: 779
- ncRNA cis-regulatory regions: 288

ori sequences:

- oriC/V: 6,690
- oriT: 502

To provide FAIR annotations, the database releases are SemVer versioned (w/o patch level), *i.e.* `<major>.<minor>`. For each version we provide a comprehensive log file tracking all imported sequences as well as annotations thereof. The DB schema is represented by the `<major>` digit and automatically checked at runtime by Bakta in order to ensure compatibility. Content updates are tracked by the `<minor>` digit.

As this taxonomic-untargeted database is fairly demanding in terms of storage consumption, we also provide a lightweight DB type providing all non-coding feature information but only PSCC information from UniRef50 clusters for CDS. If download bandwiths or storage requirements become an issue or if shorter runtimes are favored over more-specific annotation, the `light` DB will do the job.

Latest database version: 6.0
DB types:

- `light`: 1.3 Gb zipped, 3.9 Gb unzipped, MD5: 4a6e059ded39e9c5537ef4137d2f5648
- `full`: 30 Gb zipped, 84 Gb unzipped, MD5: 4c1115e40abfa2b464ae5dd988bdd88e

All database releases are hosted at Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14916843.svg)](https://doi.org/10.5281/zenodo.14916843)

## Genome Submission

Most genomes annotated with Bakta should be ready-to-submid to INSDC member databases GenBank and ENA. As a first step, please register your BioProject (e.g. PRJNA123456) and your locus_tag prefix (*e.g.* ESAKAI).

```bash
# annotate your genome in `--compliant` mode:
$ bakta --db <db-path> -v --genus Escherichia --species "coli O157:H7" --strain Sakai --complete --compliant --locus-tag ESAKAI test/data/GCF_000008865.2.fna.gz
```

### GenBank

Genomes are submitted to GenBank via Fasta (`.fna`) and SQN files. Therefore, `.sqn` files can be created with NCBI's new [table2asn](https://ftp.ncbi.nlm.nih.gov/asn1-converters/by_program/table2asn/) tool via Bakta's `.gff3` files.
Please, have a look at the [documentation]((https://www.ncbi.nlm.nih.gov/genbank/genomes_gff)) and have all additional files (template.txt) prepared:

```bash
# download table2asn for Linux
$ wget https://ftp.ncbi.nlm.nih.gov/asn1-converters/by_program/table2asn/linux64.table2asn.gz
$ gunzip linux64.table2asn.gz

# or MacOS
$ wget https://ftp.ncbi.nlm.nih.gov/asn1-converters/by_program/table2asn/mac.table2asn.gz
$ gunzip mac.table2asn.gz

$ chmod 755 linux64.table2asn.gz mac.table2asn.gz

# create the SQN file:
$ linux64.table2asn -Z -W -M n -J -c w -t template.txt -V vbt -l paired-ends -i GCF_000008865.2.fna -f GCF_000008865.2.gff3 -o GCF_000008865.2.sqn
```

### ENA

Genomes are submitted to ENA as EMBL (`.embl`) files via EBI's [Webin-CLI](https://ena-docs.readthedocs.io/en/latest/submit/general-guide/webin-cli.html) tool.
Please have all additional files (manifest.tsv, chrom-list.tsv) prepared as described [here](https://ena-docs.readthedocs.io/en/latest/submit/fileprep/assembly.html#flat-file).

```bash
# download ENA Webin-CLI
$ wget https://github.com/enasequence/webin-cli/releases/download/8.1.0/webin-cli-8.1.0.jar

$ gzip -k GCF_000008865.2.embl
$ gzip -k chrom-list.tsv
$ java -jar webin-cli-8.1.0.jar -submit -userName=<LOGIN> -password <PWD> -context genome -manifest manifest.tsv
```

Exemplarey manifest.tsv and chrom-list.tsv files might look like:

```bash
$ cat manifest.tsv
STUDY    PRJEB44484
SAMPLE    ERS6291240
ASSEMBLYNAME    GCF
ASSEMBLY_TYPE    isolate
COVERAGE    100
PROGRAM    SPAdes
PLATFORM    Illumina
MOLECULETYPE    genomic DNA
FLATFILE    GCF_000008865.2.embl.gz
CHROMOSOME_LIST    chrom-list.tsv.gz

$ cat chrom-list.tsv
contig_1    contig_1    circular-chromosome
contig_2    contig_2    circular-plasmid
contig_3    contig_3    circular-plasmid
```

## Protein bulk annotation

For the direct bulk annotation of protein sequences aside from the genome, Bakta provides a dedicated CLI entry point `bakta_proteins`:

Examples:

```bash
bakta_proteins --db <db-path> input.fasta

bakta_proteins --db <db-path> --prefix test --output test --proteins special.faa --threads 8 input.fasta
```

### Output

Annotation results are provided in standard bioinformatics file formats:

- `<prefix>.tsv`: annotations as simple human readble TSV
- `<prefix>.faa`: protein sequences as FASTA
- `<prefix>.hypotheticals.tsv`: further information on hypothetical proteins as simple human readble tab separated values
- `<prefix>.json`: all (internal) annotation & sequence information as JSON

The `<prefix>` can be set via `--prefix <prefix>`. If no prefix is set, Bakta uses the input file prefix.

### Usage

```bash
usage: bakta_proteins [--db DB] [--output OUTPUT] [--prefix PREFIX] [--force]
                      [--proteins PROTEINS]
                      [--help] [--verbose] [--debug] [--threads THREADS] [--tmp-dir TMP_DIR] [--version]
                      <input>

Rapid & standardized annotation of bacterial genomes, MAGs & plasmids

positional arguments:
  <input>               Protein sequences in (zipped) fasta format

Input / Output:
  --db DB, -d DB        Database path (default = <bakta_path>/db). Can also be provided as BAKTA_DB environment variable.
  --output OUTPUT, -o OUTPUT
                        Output directory (default = current working directory)
  --prefix PREFIX, -p PREFIX
                        Prefix for output files
  --force, -f           Force overwriting existing output folder

Annotation:
  --proteins PROTEINS   Fasta file of trusted protein sequences for annotation

General:
  --help, -h            Show this help message and exit
  --verbose, -v         Print verbose information
  --debug               Run Bakta in debug mode. Temp data will not be removed.
  --threads THREADS, -t THREADS
                        Number of threads to use (default = number of available CPUs)
  --tmp-dir TMP_DIR     Location for temporary files (default = system dependent auto detection)
  --version, -V         show program's version number and exit
```

## Genome plots

Bakta allows the creation of circular genome plots via [pyCirclize](https://github.com/moshi4/pyCirclize). Plots are generated as part of the default workflow and saved as `PNG` and `SVG` files. In addition to the default workflow, Bakta provides a dedicated CLI entry point `bakta_plot`:

Examples:

```bash
bakta_plot input.json

bakta_plot --output test --prefix test --config config.yaml --sequences 1,2 input.json
```

It accepts the results of a former annotation process in JSON format and allows the selection of distinct sequences, either denoted by their `FASTA` identifiers or sequential number starting by 1. Colors for each feature type can be adopted via a simple configuration file in `YAML` format, *e.g.* [config.yaml](config.yaml). Currently, two default plot types are supported, *i.e.* `features` and `cog`. Examples for chromosomes and plasmids are provided in [here](examples/)

### Usage

```bash
usage: bakta_plot [--config CONFIG] [--output OUTPUT] [--prefix PREFIX]
                  [--sequences SEQUENCES] [--type {features,cog}] [--label LABEL] [--size {4,8,16}] [--dpi {150,300,600}]
                  [--help] [--verbose] [--debug] [--tmp-dir TMP_DIR] [--version]
                  <input>

Rapid & standardized annotation of bacterial genomes, MAGs & plasmids

positional arguments:
  <input>               Bakta annotations in (zipped) JSON format

Input / Output:
  --config CONFIG, -c CONFIG
                        Plotting configuration in YAML format
  --output OUTPUT, -o OUTPUT
                        Output directory (default = current working directory)
  --prefix PREFIX, -p PREFIX
                        Prefix for output files

Plotting:
  --sequences SEQUENCES
                        Sequences to plot: comma separated number or name (default = all, numbers one-based)
  --type {features,cog}
                        Plot type: feature/cog (default = features)
  --label LABEL         Plot center label (for line breaks use '|')
  --size {4,8,16}       Plot size in inches: 4/8/16 (default = 8)
  --dpi {150,300,600}   Plot resolution as dots per inch: 150/300/600 (default = 300)

General:
  --help, -h            Show this help message and exit
  --verbose, -v         Print verbose information
  --debug               Run Bakta in debug mode. Temp data will not be removed.
  --tmp-dir TMP_DIR     Location for temporary files (default = system dependent auto detection)
  --version             show program's version number and exit
```

### Description

Currently, there are two types of plots: `features` (the default) and `cog`. In default mode (`features`), all features are plotted on two rings representing the forward and reverse strand from outer to inner, respectively using the following feature colors:

- CDS: `#cccccc`
- tRNA/tmRNA: `#b2df8a`
- rRNA: `#fb8072`
- ncRNA: `#fdb462`
- ncRNA-region: `#80b1d3`
- CRISPR: `#bebada`
- Gap: `#000000`
- Misc: `#666666`

In the `cog` mode, all protein-coding genes (CDS) are colored due to assigned COG functional categories. To better distinguish non-coding genes, these are plotted on an additional 3rd ring.

In addition, both plot types share two innermost GC content and GC skew rings. The first ring represents the GC content per sliding window over the entire sequence(s) in green (`#33a02c`) and red `#e31a1c` representing GC above and below average, respectively. The 2nd ring represents the GC skew in orange (`#fdbf6f`) and blue (`#1f78b4`). The GC skew gives hints on a replicon's replication bubble and hence, on the completeness of the assembly. On a complete & circular bacterial chromosome, you normally see two inflection points at the origin of replication and at its opposite region -> [Wikipedia](https://en.wikipedia.org/wiki/GC_skew)

Custom plot labels (text in the center) can be provided via `--label`:

```bash
bakta_plot --sequences 2 --dpi 300 --size 8 --prefix plot-cog-p2 --type cog --label="pO157|plasmid, 92.7 kbp"
```

![Plot example of Bakta test genome.](/examples/plot-cog-p2.png)

## Auxiliary scripts

Often, the usage of Bakta is a necessary upfront task followed by deeper analyses implemented in custom scripts. In [scripts](scripts) we'd like to collect & offer a pool of scripts addressing common tasks:

- `collect-annotation-stats.py`: Collect annotation stats for a cohort of genomes and print a condensed `TSV`.
- `extract-region.py`: Extract genome features within a given genomic range and export them as `GFF3`, `Embl`, `Genbank`, `FAA` and `FFN`

Of course, pull requests are welcome ;-)

## Web version

For further convenience, we developed an accompanying web application available at https://bakta.computational.bio.

This interactive web application provides an interactive genome browsers, aggregated feature counts and a searchable data table with detailed information on each predicted feature as well as dbxref-linked records to public databases.

Of note, this web application can also be used to visualize offline annotation results conducted by using the command line version. Therefore, the web application provides an offline viewer accepting JSON result files which are parsed and visualized locally within the browser without sending any data to the server.

## Citation

If you use Bakta in your research, please cite this paper:
> Schwengers O., Jelonek L., Dieckmann M. A., Beyvers S., Blom J., Goesmann A. (2021). Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification. Microbial Genomics, 7(11). https://doi.org/10.1099/mgen.0.000685

Bakta is *standing on the shoulder of giants* taking advantage of many great software tools and databases. If you find any of these useful for your research, please cite these primary sources, as well.

### Tools

- tRNAscan-SE 2.0 <https://doi.org/10.1093/nar/gkab688>
- Aragorn <https://doi.org/10.1093/nar/gkh152>
- Infernal <https://doi.org/10.1093/bioinformatics/btt509>
- PilerCR <https://doi.org/10.1186/1471-2105-8-18>
- Pyrodigal <https://doi.org/10.21105/joss.04296> Prodigal <https://doi.org/10.1186/1471-2105-11-119>
- Diamond <https://doi.org/10.1038/s41592-021-01101-x>
- BLAST+ <https://doi.org/10.1186/1471-2105-10-421>
- PyHMMER <https://doi.org/10.21105/joss.04296> HMMER <https://doi.org/10.1371/journal.pcbi.1002195>
- AMRFinderPlus <https://doi.org/10.1038/s41598-021-91456-0>
- pyCirclize https://github.com/moshi4/pyCirclize

### Databases

- Rfam: <https://doi.org/10.1002/cpbi.51>
- Mob-suite: <https://doi.org/10.1099/mgen.0.000206>
- DoriC: <https://doi.org/10.1093/nar/gky1014>
- AntiFam: <https://doi.org/10.1093/database/bas003>
- UniProt: <https://doi.org/10.1093/nar/gky1049>
- RefSeq: <https://doi.org/10.1093/nar/gkx1068>
- COG: <https://doi.org/10.1093/bib/bbx117>
- KEGG: <https://doi.org/10.1093/bioinformatics/btz859>
- PHROG: <https://doi.org/10.1093/nargab/lqab067>
- AMRFinder: <https://doi.org/10.1128/AAC.00483-19>
- ISFinder: <https://doi.org/10.1093/nar/gkj014>
- Pfam: <https://doi.org/10.1093/nar/gky995>
- VFDB: <https://doi.org/10.1093/nar/gky1080>

## FAQ

- **AMRFinder fails**
If AMRFinder constantly crashes even on fresh setups and Bakta's database was downloaded manually, then AMRFinder needs to setup its own internal database. This is required only once: `amrfinder_update --force_update --database <bakta-db>/amrfinderplus-db`. You could also try Bakta's internal database download logic automatically taking care of this: `bakta_db download --output <bakta-db>`

- **DeepSig not found in Conda environment**
For the prediction of signal predictions, Bakta uses DeepSig that is currently not available for MacOS and only up to Bakta v1.9.4. Therefore, we decided to exclude DeepSig from Bakta's default Conda dependencies because otherwise it would not be installable on MacOS systems. On Linux systems it can be installed via `conda install -c conda-forge -c bioconda python=3.8 deepsig`.

- **Nice, but I'm mising XYZ...**
Bakta is quite new and we're keen to constantly improve it and further expand its feature set. In case there's anything missing, please do not hesitate to open an issue and ask for it!

- **Bakta is running too long without CPU load... why?**
Bakta takes advantage of an SQLite DB which results in high storage IO loads. If this DB is stored on a remote / network volume, the lookup of IPS/PSC annotations might take a long time. In these cases, please, consider moving the DB to a local volume or hard drive.

## Issues and Feature Requests

Bakta is new and like in every software, expect some bugs lurking around. So, if you run into any issues with Bakta, we'd be happy to hear about it.
Therefore, please, execute bakta in debug mode (`--debug`) and do not hesitate to file an issue including as much information as possible:

- a detailed description of the issue
- command line output
- log file (`<prefix>.log`)
- result file (`<prefix>.json`) *if possible*
- a reproducible example of the issue with an input file that you can share *if possible*
