# Archaeal annotation support plan

## Scope and compatibility

Archaeal annotation is an explicit `--domain archaea` profile. `bacteria`
remains the default and must retain its current results and database layout.
One database build always creates both profiles; the profile is selected by the
user at runtime and is not inferred from sequence content.

## Implemented

- Domain-aware CLI configuration, JSON metadata, database checks, and protein
  bulk-annotation configuration.
- Archaeal tRNA mode, rRNA models, and domain-selected Rfam ncRNA/regulatory
  models.
- Domain-specific tmRNA labels, sORF initiation policy, and Sec/Pyl recoding
  policy.
- Archaeal oriC prediction through Ori-Finder-Arch; bacterial oriC/oriT stays
  on the existing BLAST/DoriC/MOB-suite path.
- UniRef/Swiss-Prot cluster builders accept multiple root taxa and include both
  NCBI taxon `2` (Bacteria) and `2157` (Archaea).
- An online Ori-Finder-Arch installer plus MEME and HMMER environment support.
- Bacterial DoriC/MOB-suite, AMRFinderPlus, and ISfinder assets are retained
  in the shared database, but are not used by the archaeal runtime profile.
- Ori-Finder-Arch provenance, RIP/origin attributes, scores, and tRNA intron
  coordinates are retained in GFF and INSDC exports.

## Remaining implementation work

### Database and functional annotation

1. Build, version, and publish dual-profile full and light database archives.
   Required bacterial and archaeal assets are checked by `bakta.db`; no
   production archive exists yet.
2. Curate a versioned archaeal expert-protein FASTA and supply it through
   `BAKTA_ARCHAEAL_EXPERT_PROTEINS` during database construction.
3. Replace the bacterial AntiFam-only spurious-ORF filter with a validated
   archaeal or domain-neutral model.
4. Review selected archaeal Rfam families and thresholds for precision; do not
   assume that bacterial regulatory-RNA selection rules transfer unchanged.

### Feature fidelity

1. Validate tRNA intron coordinates and exports against curated archaeal
   records.
2. Validate tmRNA calls and sORF start-codon policies across archaeal clades.
3. Validate Sec/Pyl recoding calls against known loci and add suitable Rfam
   model coverage where required.
4. Validate Ori-Finder-Arch provenance, RIP type, and origin class in every
   output format against downstream consumers.

### Validation and release

1. Add fixture databases and curated archaeal genomes for end-to-end tests.
2. Benchmark CDS boundaries, tRNAs, rRNAs, ncRNAs, pseudogenes, and oriCs;
   report recall/precision by archaeal lineage and assembly quality.
3. Pin and checksum the upstream Ori-Finder-Arch binary release.
4. Update CLI/package wording, manuals, examples, and release notes.

## Acceptance criteria

The archaeal profile is release-ready only when a versioned database is
available, all profile dependencies are installable, end-to-end tests cover
complete and fragmented archaeal assemblies, and benchmark results meet the
project's agreed accuracy thresholds.
