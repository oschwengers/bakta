#!/usr/bin/env bash

set -euo pipefail

# Build one database containing both domain profiles. The runtime
# --domain option selects its profile-specific assets and predictors.
BAKTA_TAXON_ARGS="--taxon 2 --taxon 2157"

# Database helper scripts reside beside this entry point. This used to rely on
# an undocumented BAKTA_DB_SCRIPTS environment variable, making direct runs
# fail with paths such as /init-db.py. Resolve it from the script location so
# the build can be started from any working directory.
BAKTA_DB_SCRIPTS="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"

mkdir -p db
cd db

# Guard against two invocations building into the same db/ directory at once:
# concurrent runs share unnamespaced intermediate files (e.g. rfam-genes.txt)
# and can silently interleave writes into them, corrupting the build.
exec 9>.build.lock
if ! flock -n 9; then
    printf "Another buid-db.sh build is already running against %s -- aborting.\n" "$(pwd)" >&2
    # Distinct from a genuine failure (exit 1): callers that auto-resubmit on
    # failure (e.g. a PBS orchestrator job) should treat this as "another
    # instance already has it covered", not as something to retry.
    exit 99
fi

# Every numbered step below is guarded by a .stepNN.done marker, touched only
# as the very last action of a successful step. A killed/crashed run can
# therefore be resumed by simply re-invoking this script from the same db/
# directory: completed steps are skipped, and steps 1-4 additionally skip
# already-finished bacteria/archaea profiles within an interrupted step.
have_outputs() {
    local f
    for f in "$@"; do
        [ -s "${f}" ] || return 1
    done
    return 0
}

if have_outputs .build.done; then
    printf "Bakta database already fully built in %s -- nothing to do.\n" "$(pwd)"
    exit 0
fi

printf "Create Bakta database\n"

# download rRNA covariance models from Rfam
printf "\n1/19: download rRNA covariance models from Rfam ...\n"
if have_outputs .step01.done; then
    printf "    already built, skipping.\n"
else
    if [ ! -s Rfam.cm ]; then
        wget -c https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
        pigz -df Rfam.cm.gz
    fi
    for profile in bacteria archaea; do
        if [ "${profile}" = "bacteria" ]; then
            rRNA_name=rRNA
            rRNA_ssu=RF00177
            rRNA_lsu=RF02541
        else
            rRNA_name=rRNA-archaea
            rRNA_ssu=RF01959
            rRNA_lsu=RF02540
        fi
        if have_outputs "${rRNA_name}.i1m" "${rRNA_name}.i1i" "${rRNA_name}.i1f" "${rRNA_name}.i1p"; then
            continue
        fi
        cmfetch Rfam.cm RF00001 > "${rRNA_name}"
        cmfetch Rfam.cm "${rRNA_ssu}" >> "${rRNA_name}"
        cmfetch Rfam.cm "${rRNA_lsu}" >> "${rRNA_name}"
        cmpress -F "${rRNA_name}"
        rm "${rRNA_name}"
    done
    date > .step01.done
fi


# download and extract ncRNA gene covariance models from Rfam
printf "\n2/19: download ncRNA gene covariance models from Rfam ...\n"
if have_outputs .step02.done; then
    printf "    already built, skipping.\n"
else
    if [ ! -s Rfam.cm ]; then
        wget -c https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
        pigz -df Rfam.cm.gz
    fi
    for profile in bacteria archaea; do
        if [ "${profile}" = "bacteria" ]; then
            RFAM_TAXON=Bacteria
            NCRNA_GENES_NAME=ncRNA-genes
        else
            RFAM_TAXON=Archaea
            NCRNA_GENES_NAME=ncRNA-genes-archaea
        fi
        if have_outputs "${NCRNA_GENES_NAME}.i1m" "${NCRNA_GENES_NAME}.i1i" "${NCRNA_GENES_NAME}.i1f" "${NCRNA_GENES_NAME}.i1p"; then
            continue
        fi
        sed "s/LIKE 'Bacteria%'/LIKE '${RFAM_TAXON}%'/" ${BAKTA_DB_SCRIPTS}/ncRNA-genes.sql | mysql --user rfamro --host mysql-rfam-public.ebi.ac.uk --port 4497 --database Rfam | tail -n +2 > rfam-genes.raw.txt
        rm -f rfam-genes.txt
        # Not every domain has hits in every category (e.g. archaea has no
        # antitoxin entries): grep exits 1 on zero matches, which is fatal
        # under `set -e` unless guarded with `|| true`.
        grep "antitoxin;" rfam-genes.raw.txt >> rfam-genes.txt || true
        grep "antisense;" rfam-genes.raw.txt >> rfam-genes.txt || true
        grep "ribozyme;" rfam-genes.raw.txt >> rfam-genes.txt || true
        grep "sRNA;" rfam-genes.raw.txt >> rfam-genes.txt || true
        cut -f1 ${BAKTA_DB_SCRIPTS}/ncRNA-genes.blocklist.txt > ncRNA-genes.blocklist
        grep -e "Gene;[^ ]" rfam-genes.raw.txt | grep -v -f ncRNA-genes.blocklist >> rfam-genes.txt || true
        sort -u rfam-genes.txt > rfam-genes.uniq.txt
        cmfetch -o "${NCRNA_GENES_NAME}" -f Rfam.cm rfam-genes.uniq.txt
        cmpress -F "${NCRNA_GENES_NAME}"
        rm rfam-genes.raw.txt rfam-genes.txt rfam-genes.uniq.txt ncRNA-genes.blocklist "${NCRNA_GENES_NAME}"
    done
    if [ ! -s rfam-go.tsv ]; then
        wget -c http://current.geneontology.org/ontology/external2go/rfam2go
        awk -F ' ' '{print $1 "\t" $NF}' rfam2go > rfam-go.tsv
        rm rfam2go
    fi
    date > .step02.done
fi


# download and extract ncRNA regions (cis reg elements) covariance models from Rfam
printf "\n3/19: download ncRNA region covariance models from Rfam ...\n"
if have_outputs .step03.done; then
    printf "    already built, skipping.\n"
else
    if [ ! -s Rfam.cm ]; then
        wget -c https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
        pigz -df Rfam.cm.gz
    fi
    for profile in bacteria archaea; do
        if [ "${profile}" = "bacteria" ]; then
            RFAM_TAXON=Bacteria
            NCRNA_REGIONS_NAME=ncRNA-regions
        else
            RFAM_TAXON=Archaea
            NCRNA_REGIONS_NAME=ncRNA-regions-archaea
        fi
        if have_outputs "${NCRNA_REGIONS_NAME}.i1m" "${NCRNA_REGIONS_NAME}.i1i" "${NCRNA_REGIONS_NAME}.i1f" "${NCRNA_REGIONS_NAME}.i1p"; then
            continue
        fi
        sed "s/LIKE 'Bacteria%'/LIKE '${RFAM_TAXON}%'/" ${BAKTA_DB_SCRIPTS}/ncRNA-regions.sql | mysql --user rfamro --host mysql-rfam-public.ebi.ac.uk --port 4497 --database Rfam | tail -n +2 > rfam-regions.raw.txt
        rm -f rfam-regions.txt
        grep "riboswitch;" rfam-regions.raw.txt >> rfam-regions.txt || true
        grep "thermoregulator;" rfam-regions.raw.txt >> rfam-regions.txt || true
        grep "leader;" rfam-regions.raw.txt >> rfam-regions.txt || true
        grep "frameshift_element;" rfam-regions.raw.txt >> rfam-regions.txt || true
        cut -f1 ${BAKTA_DB_SCRIPTS}/ncRNA-regions.blocklist.txt > ncRNA-regions.blocklist
        grep -e "Cis-reg;[^ ]" rfam-regions.raw.txt | grep -v -f ncRNA-regions.blocklist >> rfam-regions.txt || true
        sort -u rfam-regions.txt > rfam-regions.uniq.txt
        cmfetch -o "${NCRNA_REGIONS_NAME}" -f Rfam.cm rfam-regions.uniq.txt
        cmpress -F "${NCRNA_REGIONS_NAME}"
        rm rfam-regions.raw.txt rfam-regions.txt rfam-regions.uniq.txt ncRNA-regions.blocklist "${NCRNA_REGIONS_NAME}"
    done
    date > .step03.done
fi
rm -f Rfam.cm


# download and extract spurious ORF HMMs from AntiFam
printf "\n4/19: download and extract spurious ORF HMMs from AntiFam ...\n"
if have_outputs .step04.done; then
    printf "    already built, skipping.\n"
else
    if [ ! -s antifam-dir/AntiFam_Bacteria.hmm ] || [ ! -s antifam-dir/AntiFam_Archaea.hmm ]; then
        rm -rf antifam-dir
        mkdir antifam-dir
        (cd antifam-dir && wget -c https://ftp.ebi.ac.uk/pub/databases/Pfam/AntiFam/current/Antifam.tar.gz && tar -xzf Antifam.tar.gz)
    fi
    if ! have_outputs antifam.h3m antifam.h3i antifam.h3f antifam.h3p; then
        cp antifam-dir/AntiFam_Bacteria.hmm antifam
        hmmpress -f antifam
    fi
    if ! have_outputs antifam-archaea.h3m antifam-archaea.h3i antifam-archaea.h3f antifam-archaea.h3p; then
        cp antifam-dir/AntiFam_Archaea.hmm antifam-archaea
        hmmpress -f antifam-archaea
    fi
    rm -f antifam antifam-archaea
    rm -rf antifam-dir/
    date > .step04.done
fi


# DoriC and MOB-suite records are retained for the bacterial profile only;
# archaeal annotation uses Ori-Finder-Arch at runtime.
printf "\n5/19: download and extract oriT sequences from Mob-suite ...\n"
if have_outputs orit.fna; then
    printf "    already built, skipping.\n"
else
    wget -c https://zenodo.org/records/10304948/files/data.tar.gz
    tar -xvzf data.tar.gz
    mv data/orit.fas ./orit.fna
    rm -r data/ data.tar.gz
fi
printf "\n5/19: download oriC/V sequences from DoriC ...\n"
if have_outputs oric.fna; then
    printf "    already built, skipping.\n"
else
    curl 'https://tubic.org/doric/search/bacteria' \
      -H 'content-type: multipart/form-data; boundary=----WebKitFormBoundaryBDBZTWpS3orCjS0m' \
      --data-raw $'------WebKitFormBoundaryBDBZTWpS3orCjS0m\r\nContent-Disposition: form-data; name="assembly_level"\r\n\r\nComplete\r\n------WebKitFormBoundaryBDBZTWpS3orCjS0m\r\nContent-Disposition: form-data; name="topology"\r\n\r\nAll\r\n------WebKitFormBoundaryBDBZTWpS3orCjS0m\r\nContent-Disposition: form-data; name="chromosome_type"\r\n\r\nAll\r\n------WebKitFormBoundaryBDBZTWpS3orCjS0m\r\nContent-Disposition: form-data; name="oric_type"\r\n\r\nSingle\r\n------WebKitFormBoundaryBDBZTWpS3orCjS0m\r\nContent-Disposition: form-data; name="organism"\r\n\r\n\r\n------WebKitFormBoundaryBDBZTWpS3orCjS0m\r\nContent-Disposition: form-data; name="lineage"\r\n\r\n\r\n------WebKitFormBoundaryBDBZTWpS3orCjS0m\r\nContent-Disposition: form-data; name="download1"\r\n\r\nDownload\r\n------WebKitFormBoundaryBDBZTWpS3orCjS0m--\r\n' \
      --compressed > oric.csv
    curl 'https://tubic.org/doric/search/plasmid' \
      -H 'content-type: multipart/form-data; boundary=----WebKitFormBoundaryMu32WgFUyqC7TO0d' \
      --data-raw $'------WebKitFormBoundaryMu32WgFUyqC7TO0d\r\nContent-Disposition: form-data; name="topology"\r\n\r\nAll\r\n------WebKitFormBoundaryMu32WgFUyqC7TO0d\r\nContent-Disposition: form-data; name="organism"\r\n\r\n\r\n------WebKitFormBoundaryMu32WgFUyqC7TO0d\r\nContent-Disposition: form-data; name="lineage"\r\n\r\n\r\n------WebKitFormBoundaryMu32WgFUyqC7TO0d\r\nContent-Disposition: form-data; name="download1"\r\n\r\nDownload\r\n------WebKitFormBoundaryMu32WgFUyqC7TO0d--\r\n' \
      --compressed > oriv.csv
    python3 ${BAKTA_DB_SCRIPTS}/extract-ori.py --doric oric.csv --fasta ori.chromosome.fna
    python3 ${BAKTA_DB_SCRIPTS}/extract-ori.py --doric oriv.csv --fasta ori.plasmid.fna
    cat ori.chromosome.fna > oric.raw.fna
    cat ori.plasmid.fna >> oric.raw.fna
    cd-hit-est -i oric.raw.fna -o oric.fna -c 0.99 -s 0.99 -aS 0.99 -g 1 -r 1
    rm *.csv ori.*.fna oric.raw.fna oric.fna.clstr
fi


# download NCBI Taxonomy DB
printf "\n6/19: download NCBI Taxonomy DB ...\n"
if have_outputs nodes.dmp; then
    printf "    already built, skipping.\n"
else
    rm -rf taxonomy
    mkdir taxonomy
    (cd taxonomy && wget -c https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz && tar -I pigz -xf taxdump.tar.gz)
    mv taxonomy/nodes.dmp .
    rm -rf taxonomy
fi


############################################################################
# Setup SQLite Bakta db
############################################################################
printf "\n7/19: setup SQLite Bakta db ...\n"
if have_outputs bakta.db; then
    printf "    already built, skipping.\n"
else
    python3 ${BAKTA_DB_SCRIPTS}/init-db.py --db bakta.db
fi


############################################################################
# Build protein sequence clusters (PSCCs) based on UniRef50 entries
# - download UniProt UniRef50
# - read and transform UniRef50 XML file to DB and Fasta file
# - build PSCC Diamond db
############################################################################
printf "\n8/19: download UniProt UniRef50 ...\n"
if have_outputs .step08.done; then
    printf "    already built, skipping.\n"
else
    wget -c https://ftp.expasy.org/databases/uniprot/current_release/uniref/uniref50/uniref50.xml.gz
    printf "\n8/19: download and merge UniParc active sequences (200 parts, resumable) ...\n"
    # Each part is merged and deleted immediately, one at a time, so at
    # most one ~1.2GB compressed part sits on disk alongside the growing
    # merged fasta -- downloading all 200 parts (~250GB) before merging
    # would nearly double peak disk usage for no benefit. This whole block
    # is gated only by the .uniparc_parts.merged manifest (not by whether
    # uniparc_active.fasta is merely non-empty, since it grows incrementally
    # and would look "done" long before all 200 parts are in it). wget -c
    # runs unconditionally so a truncated leftover from an earlier
    # interrupted download is resumed and completed rather than silently
    # treated as a finished part.
    touch .uniparc_parts.merged
    for i in {1..200}; do
        grep -qx "${i}" .uniparc_parts.merged && continue
        part="uniparc_active_p${i}.fasta.gz"
        wget -c "https://ftp.expasy.org/databases/uniprot/current_release/uniparc/fasta/active/${part}"
        pigz -dc "${part}" >> uniparc_active.fasta
        echo "${i}" >> .uniparc_parts.merged
        rm -f "${part}"
    done
    # .uniparc_parts.merged is deliberately NOT removed here: if the python
    # step below fails, a retry must still see the merge as complete rather
    # than starting the 200-part loop over and re-appending already-merged
    # parts on top of the existing uniparc_active.fasta.
    printf "\n8/19: read UniRef90 entries and build Protein Sequence Cluster sequence and information databases:\n"
    python3 ${BAKTA_DB_SCRIPTS}/init-pscc.py --taxonomy nodes.dmp --uniref50 uniref50.xml.gz --uniparc uniparc_active.fasta --db bakta.db --pscc pscc.faa --pscc_sorf pscc_sorf.faa ${BAKTA_TAXON_ARGS}
    printf "\n8/19: build PSCC Diamond db ...\n"
    diamond makedb --in pscc.faa --db pscc
    diamond makedb --in pscc_sorf.faa --db sorf
    mkdir -p db-light
    cp bakta.db db-light/
    mv pscc.dmnd sorf.dmnd db-light/
    (cd db-light && python3 ${BAKTA_DB_SCRIPTS}/optimize-db.py --db bakta.db)
    rm uniref50.xml.gz
    date > .step08.done
fi


############################################################################
# Build protein sequence clusters (PSCs) based on UniRef90 entries
# - download UniProt UniRef90
# - read and transform UniRef90 XML file to DB and Fasta file
# - build PSC Diamond db
############################################################################
printf "\n9/19: download UniProt UniRef90 ...\n"
if have_outputs .step09.done; then
    printf "    already built, skipping.\n"
else
    wget -c https://ftp.expasy.org/databases/uniprot/current_release/uniref/uniref90/uniref90.xml.gz
    printf "\n9/19: read UniRef90 entries and build Protein Sequence Cluster sequence and information databases:\n"
    python3 ${BAKTA_DB_SCRIPTS}/init-psc.py --taxonomy nodes.dmp --uniref90 uniref90.xml.gz --uniparc uniparc_active.fasta --db bakta.db --psc psc.faa --psc_sorf sorf.faa ${BAKTA_TAXON_ARGS}
    printf "\n9/19: build PSC Diamond db ...\n"
    diamond makedb --in psc.faa --db psc
    diamond makedb --in sorf.faa --db sorf
    rm uniref90.xml.gz
    date > .step09.done
fi


############################################################################
# Build unique protein sequences (IPSs) based on UniRef100 entries
# - download UniProt UniRef100
# - read, filter and transform UniRef100 entries and store to ips.db
############################################################################
printf "\n10/19: download UniProt UniRef100 ...\n"
if have_outputs .step10.done; then
    printf "    already built, skipping.\n"
else
    wget -c https://ftp.expasy.org/databases/uniprot/current_release/uniref/uniref100/uniref100.xml.gz
    printf "\n10/19: read, filter and store UniRef100 entries ...:\n"
    python3 ${BAKTA_DB_SCRIPTS}/init-ups-ips.py --taxonomy nodes.dmp --uniref100 uniref100.xml.gz --uniparc uniparc_active.fasta --db bakta.db --ips ips.faa ${BAKTA_TAXON_ARGS}
    rm uniref100.xml.gz uniparc_active.fasta
    date > .step10.done
fi


############################################################################
# Integrate NCBI nonredundant protein identifiers and COG db
# - download bacterial and archaeal RefSeq nonredundant proteins and COG files
# - annotate UPSs with NCBI nrp IDs (WP_*)
# - annotate IPSs/PSCs with COG IDs, gene symbols, product descriptions (seq -> hash -> UniParc/WP_* -> UniRef100 -> UniRef90 -> PSC)
############################################################################
printf "\n11/19: download NCBI COG clusters and RefSeq nonredundant proteins ...\n"
if have_outputs .step11.done; then
    printf "    already built, skipping.\n"
else
    wget -c https://ftp.ncbi.nih.gov/pub/COG/COG2024/data/cog-24.def.tab  # COG IDs and functional class
    wget -c https://ftp.ncbi.nih.gov/pub/COG/COG2024/data/cog-24.cog.csv  # Mapping GenBank IDs -> COG IDs
    mkdir -p refseq-nrp
    touch .refseq_files.merged
    for refseq_domain in bacteria archaea; do
        refseq_url=https://ftp.ncbi.nih.gov/refseq/release/${refseq_domain}
        wget -qO- "${refseq_url}/" | sed -n "s/.*href=\"\(${refseq_domain}\.wp_protein\.[0-9]*\.protein\.faa\.gz\)\".*/\1/p" > "refseq-nrp/${refseq_domain}.files.txt"
        # Each file is downloaded, merged and deleted immediately (like the
        # UniParc loop above) so only one file sits on disk at a time instead
        # of the whole RefSeq release alongside the merged fasta.
        while read -r refseq_file; do
            grep -qx "${refseq_file}" .refseq_files.merged && continue
            wget -c -O "refseq-nrp/${refseq_file}" "${refseq_url}/${refseq_file}"
            pigz -dc "refseq-nrp/${refseq_file}" | seqtk seq -CU >> refseq-nrp.trimmed.faa
            echo "${refseq_file}" >> .refseq_files.merged
            rm -f "refseq-nrp/${refseq_file}"
        done < "refseq-nrp/${refseq_domain}.files.txt"
    done
    rm -rf refseq-nrp
    # .refseq_files.merged is deliberately kept for the same reason as
    # .uniparc_parts.merged above: if the python step below fails, a retry
    # must not re-merge already-processed files on top of the existing
    # refseq-nrp.trimmed.faa.
    printf "\n11/19: annotate IPSs and PSCs ...\n"
    python3 ${BAKTA_DB_SCRIPTS}/annotate-ncbi-nrp-cog.py --db bakta.db --nrp refseq-nrp.trimmed.faa --cog-ids cog-24.def.tab --nrp-cog-mapping cog-24.cog.csv
    rm refseq-nrp.trimmed.faa .refseq_files.merged
    date > .step11.done
fi


############################################################################
# Integrate KEGG kofams
# - download KEGG kofams
# - select eligible HMMs (f measure>0.77)
# - annotate PSCs
############################################################################
printf "\n12/19: download KEGG kofams HMM models...\n"
if have_outputs .step12.done; then
    printf "    already built, skipping.\n"
else
    wget -c https://www.genome.jp/ftp/db/kofam/ko_list.gz
    wget -c https://www.genome.jp/ftp/db/kofam/profiles.tar.gz
    zcat ko_list.gz | grep full | awk '{ if($5>=0.77) print $0}' > hmms.kofam.selected.tsv
    cut -f1 hmms.kofam.selected.tsv > hmms.ids.txt
    rm -rf profiles
    tar -I pigz -xf profiles.tar.gz
    rm -f kofam-prok
    for kofam in `cat profiles/prokaryote.hal`; do cat profiles/$kofam >> kofam-prok; done
    hmmfetch -f -o kofams kofam-prok hmms.ids.txt
    hmmpress -f kofams
    printf "\n12/19: annotate PSCs...\n"
    mkdir -p work/tblout work/domtblout
    nextflow run ${BAKTA_DB_SCRIPTS}/hmmsearch.nf -resume --in psc.faa --db kofams --no_tc --out hmmsearch.kofam.tblout
    python3 ${BAKTA_DB_SCRIPTS}/annotate-kofams.py --db bakta.db --hmms hmms.kofam.selected.tsv --hmm-results hmmsearch.kofam.tblout
    rm -rf profiles ko_list.gz profiles.tar.gz kofam* hmmsearch.kofam.* hmms*
    date > .step12.done
fi


############################################################################
# Integrate UniProt Swissprot information
# - download SwissProt annotation xml file
# - annotate PSCs if IPS have PSC UniRef90 identifier (seq -> hash -> UPS -> IPS -> PSC)
# - annotate IPSs if IPS have no PSC UniRef90 identifier (seq -> hash -> UPS -> IPS)
############################################################################
printf "\n13/19: download UniProt/SwissProt ...\n"
if have_outputs .step13.done; then
    printf "    already built, skipping.\n"
else
    wget -c https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.xml.gz
    printf "\n13/19: annotate IPSs and PSCs ...\n"
    python3 ${BAKTA_DB_SCRIPTS}/annotate-swissprot.py --taxonomy nodes.dmp --xml uniprot_sprot.xml.gz --db bakta.db ${BAKTA_TAXON_ARGS}
    rm uniprot_sprot.xml.gz
    date > .step13.done
fi


############################################################################
# Integrate NCBIfams HMM models
# - download NCBIfams HMM models
# - annotate PSCs
############################################################################
printf "\n14/19: download NCBIfams HMM models...\n"
if have_outputs .step14.done; then
    printf "    already built, skipping.\n"
else
    wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.LIB
    wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.tsv
    grep -v "(Provisional)" hmm_PGAP.tsv > hmms.non-prov.tsv
    grep exception hmms.non-prov.tsv > hmms.ncbi.selected.tsv || true
    grep equivalog hmms.non-prov.tsv >> hmms.ncbi.selected.tsv || true
    sort hmms.ncbi.selected.tsv | uniq | cut -f1 > hmms.ids.txt
    hmmfetch -f -o ncbifams hmm_PGAP.LIB hmms.ids.txt
    hmmpress -f ncbifams
    printf "\n14/19: annotate PSCs...\n"
    mkdir -p work/tblout work/domtblout
    nextflow run ${BAKTA_DB_SCRIPTS}/hmmsearch.nf -resume --in psc.faa --db ncbifams --block 10000 --out hmmsearch.ncbifams.tblout
    python3 ${BAKTA_DB_SCRIPTS}/annotate-ncbi-fams.py --db bakta.db --hmms hmms.ncbi.selected.tsv --hmm-results hmmsearch.ncbifams.tblout
    rm ncbifams* hmms.* hmm_PGAP.* hmmsearch.ncbifams.tblout
    date > .step14.done
fi


############################################################################
# Integrate PHROG DB of phage orthologous
# - download PHROG protein sequences
# - filter unannotated PHROGs
# - annotate PSCs
############################################################################
printf "\n15/19: download PHROGs ...\n"
if have_outputs .step15.done; then
    printf "    already built, skipping.\n"
else
    wget -c https://phrogs.lmge.uca.fr/downloads_from_website/FAA_phrog.tar.gz
    wget -c https://phrogs.lmge.uca.fr/downloads_from_website/phrog_annot_v4.tsv
    rm -rf FAA_phrog
    tar -xzf FAA_phrog.tar.gz
    rm -f phrogs-raw.faa
    cat FAA_phrog/*.faa >> phrogs-raw.faa
    python3 ${BAKTA_DB_SCRIPTS}/extract-phrogs.py --annotation phrog_annot_v4.tsv --proteins phrogs-raw.faa --filtered-proteins phrogs.faa
    diamond makedb --in phrogs.faa --db phrog
    printf "\n15/19: annotate PSCs...\n"
    python3 ${BAKTA_DB_SCRIPTS}/extract-hypotheticals.py --psc psc.faa --db bakta.db --hypotheticals hypotheticals.faa
    nextflow run ${BAKTA_DB_SCRIPTS}/diamond.nf -resume --in hypotheticals.faa --db phrog.dmnd --block 100000 --id 90 --qcov 80 --scov 80 --out diamond.phrog.psc.tsv
    python3 ${BAKTA_DB_SCRIPTS}/annotate-phrogs.py --db bakta.db --annotation phrog_annot_v4.tsv --psc-alignments diamond.phrog.psc.tsv
    rm -r FAA_phrog.tar.gz phrog_annot_v4.tsv FAA_phrog phrogs-raw.faa phrogs.faa phrog.dmnd hypotheticals.faa diamond.phrog.psc.tsv
    date > .step15.done
fi


############################################################################
# Integrate NCBI Pathogen AMR db
# - download AMR gene WP_* annotations from NCBI Pathogen ReferenceGeneCatalog
# - annotate IPSs with AMR info
############################################################################
printf "\n16/19: download AMR gene WP_* annotations from NCBI Pathogen AMR db ...\n"
if have_outputs .step16.done; then
    printf "    already built, skipping.\n"
else
    wget -c https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/ReferenceGeneCatalog.txt
    printf "\n16/19: annotate PSCs...\n"
    python3 ${BAKTA_DB_SCRIPTS}/annotate-ncbi-amr.py --db bakta.db --genes ReferenceGeneCatalog.txt
    rm ReferenceGeneCatalog.txt
    date > .step16.done
fi


############################################################################
# Integrate ISfinder db
# - download IS protein sequences from GitHub (oschwengers/ISfinder-sequences)
# - extract IS transposase sequences and mark ORF A/B transposases
# - annotate IPSs/PCSs with IS info
############################################################################
printf "\n17/19: download & extract ISfinder protein sequences ...\n"
if have_outputs .step17.done; then
    printf "    already built, skipping.\n"
else
    wget -c https://github.com/oschwengers/ISfinder-sequences/raw/2e9162bd5e3448c86ec1549a55315e498bef72fc/IS.faa
    python3 ${BAKTA_DB_SCRIPTS}/extract-is.py --input IS.faa --output is.transposase.faa
    printf "\n17/19: annotate IPSs/PCSs ...\n"
    diamond makedb --in is.transposase.faa --db is
    nextflow run ${BAKTA_DB_SCRIPTS}/diamond.nf -resume --in ips.faa --db is.dmnd --block 100000 --id 95 --qcov 90 --scov 90 --out diamond.is.ips.tsv
    nextflow run ${BAKTA_DB_SCRIPTS}/diamond.nf -resume --in psc.faa --db is.dmnd --block 100000 --id 90 --qcov 80 --scov 80 --out diamond.is.psc.tsv
    python3 ${BAKTA_DB_SCRIPTS}/annotate-is.py --db bakta.db --ips-alignments diamond.is.ips.tsv --psc-alignments diamond.is.psc.tsv
    rm is.transposase.faa is.dmnd diamond.is.ips.tsv diamond.is.psc.tsv
    date > .step17.done
fi


############################################################################
# Integrate Pfam A
# - download all Pfam A HMM models
# - extract families & domains
# - compress HMM models
# - annotate hypothetical PSC via Pfam families
############################################################################
printf "\n18/19: download HMM models from Pfam ...\n"
if have_outputs .step18.done; then
    printf "    already built, skipping.\n"
else
    wget -c https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz
    python3 ${BAKTA_DB_SCRIPTS}/extract-pfam.py --pfam Pfam-A.hmm.dat.gz --family pfam.families.tsv --non-family pfam.non-families.tsv
    if [ ! -s Pfam-A.hmm ]; then
        if [ ! -s Pfam-A.hmm.gz ]; then
            wget -c https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
        fi
        pigz -df Pfam-A.hmm.gz
    fi
    hmmfetch -o pfam-families -f Pfam-A.hmm pfam.families.tsv
    hmmpress -f pfam-families
    hmmfetch -o pfam -f Pfam-A.hmm pfam.non-families.tsv
    hmmpress -f pfam
    python3 ${BAKTA_DB_SCRIPTS}/extract-hypotheticals.py --psc psc.faa --db bakta.db --hypotheticals hypotheticals.faa
    mkdir -p work/tblout work/domtblout
    nextflow run ${BAKTA_DB_SCRIPTS}/hmmsearch.nf -resume --in hypotheticals.faa --db pfam-families --block 10000 --out hmmsearch.pfam-families.tblout
    python3 ${BAKTA_DB_SCRIPTS}/annotate-pfam.py --db bakta.db --hmms pfam-families --hmm-results hmmsearch.pfam-families.tblout
    # `pfam.h3*` (the non-family Pfam HMMs) are a runtime bakta asset and are
    # kept; only the family HMMs used for this step's own hmmsearch, and the
    # intermediate tsv/source files, are discarded here.
    rm -f pfam-families* pfam pfam.families.tsv pfam.non-families.tsv Pfam-A.hmm.dat.gz Pfam-A.hmm hmmsearch.pfam-families.tblout hypotheticals.faa
    date > .step18.done
fi


############################################################################
# Setup expert protein sequences
# - import IS sequences
# - import NCBI BlastRules models
# - import VFDB sequences
############################################################################
printf "\n19/19: download AA sequences for expert annotation system ...\n"
if have_outputs .step19.done; then
    printf "    already built, skipping.\n"
else
    wget -c https://ftp.ncbi.nlm.nih.gov/pub/blastrules/4.2.2.tgz
    rm -rf 4.2.2
    tar -xzf 4.2.2.tgz
    if [ ! -s VFDB_setA_pro.fas ]; then
        if [ ! -s VFDB_setA_pro.fas.gz ]; then
            wget -c http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz
        fi
        gunzip -f VFDB_setA_pro.fas.gz
    fi
    rm -f expert-protein-sequences.faa
    python3 ${BAKTA_DB_SCRIPTS}/expert/setup-is.py --expert-sequence expert-protein-sequences.faa --proteins IS.faa
    python3 ${BAKTA_DB_SCRIPTS}/expert/setup-ncbiblastrules.py --expert-sequence expert-protein-sequences.faa --ncbi-blastrule-tsv 4.2.2/data/blast-rules_4.2.2.tsv --proteins 4.2.2/data/proteins.fasta
    python3 ${BAKTA_DB_SCRIPTS}/expert/setup-vfdb.py --expert-sequence expert-protein-sequences.faa --proteins VFDB_setA_pro.fas
    diamond makedb --in expert-protein-sequences.faa --db expert-protein-sequences
    rm -r 4.2.2/ 4.2.2.tgz IS.faa VFDB_setA_pro.fas expert-protein-sequences.faa
    : "${BAKTA_ARCHAEAL_EXPERT_PROTEINS:?Set BAKTA_ARCHAEAL_EXPERT_PROTEINS to a curated archaeal protein FASTA.}"
    diamond makedb --in "${BAKTA_ARCHAEAL_EXPERT_PROTEINS}" --db expert-protein-sequences-archaea
    date > .step19.done
fi

# The light database keeps its PSCC SQLite snapshot but requires the same
# profile-specific HMMs and expert databases as the full database.
cp antifam* ncRNA-genes* ncRNA-regions* rRNA* oric.fna orit.fna pfam* rfam-go.tsv expert-protein-sequences*.dmnd db-light/

# Cleanup
ls -l bakta.db
python3 ${BAKTA_DB_SCRIPTS}/optimize-db.py --db bakta.db --tmp /var/scratch/
ls -l bakta.db
rm -f psc.faa sorf.faa nodes.dmp
date > .build.done
