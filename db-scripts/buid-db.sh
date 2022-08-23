#!/bin/bash

set -e

mkdir db
cd db

printf "Create Bakta database\n"

# download rRNA covariance models from Rfam
printf "\n1/18: download rRNA covariance models from Rfam ...\n"
wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
pigz -d Rfam.cm.gz
cmfetch Rfam.cm RF00001 >  rRNA
cmfetch Rfam.cm RF00177 >> rRNA
cmfetch Rfam.cm RF02541 >> rRNA
cmpress rRNA
rm rRNA


# download and extract ncRNA gene covariance models from Rfam
printf "\n2/18: download ncRNA gene covariance models from Rfam ...\n"
mysql --user rfamro --host mysql-rfam-public.ebi.ac.uk --port 4497 --database Rfam < ${BAKTA_DB_SCRIPTS}/ncRNA-genes.sql | tail -n +2 > rfam-genes.raw.txt
grep "antitoxin;" rfam-genes.raw.txt >> rfam-genes.txt
grep "antisense;" rfam-genes.raw.txt >> rfam-genes.txt
grep "ribozyme;" rfam-genes.raw.txt >> rfam-genes.txt
grep "sRNA;" rfam-genes.raw.txt >> rfam-genes.txt
cut -f1 ${BAKTA_DB_SCRIPTS}/ncRNA-genes.blocklist.txt > ncRNA-genes.blocklist
grep -e "Gene;[^ ]" rfam-genes.raw.txt | grep -v -f ncRNA-genes.blocklist >> rfam-genes.txt
cmfetch -o ncRNA-genes -f Rfam.cm rfam-genes.txt
cmpress ncRNA-genes
wget http://current.geneontology.org/ontology/external2go/rfam2go
awk -F ' ' '{print $1 "\t" $NF}' rfam2go > rfam-go.tsv
rm rfam-genes.raw.txt rfam-genes.txt ncRNA-genes.blocklist ncRNA-genes rfam2go


# download and extract ncRNA regions (cis reg elements) covariance models from Rfam
printf "\n3/18: download ncRNA region covariance models from Rfam ...\n"
mysql --user rfamro --host mysql-rfam-public.ebi.ac.uk --port 4497 --database Rfam < ${BAKTA_DB_SCRIPTS}/ncRNA-regions.sql | tail -n +2 > rfam-regions.raw.txt
grep "riboswitch;" rfam-regions.raw.txt >> rfam-regions.txt
grep "thermoregulator;" rfam-regions.raw.txt >> rfam-regions.txt
grep "leader;" rfam-regions.raw.txt >> rfam-regions.txt
grep "frameshift_element;" rfam-regions.raw.txt >> rfam-regions.txt
cut -f1 ${BAKTA_DB_SCRIPTS}/ncRNA-regions.blocklist.txt > ncRNA-regions.blocklist
grep -e "Cis-reg;[^ ]" rfam-regions.raw.txt | grep -v -f ncRNA-regions.blocklist >> rfam-regions.txt
cmfetch -o ncRNA-regions -f Rfam.cm rfam-regions.txt
cmpress ncRNA-regions
rm rfam-regions.raw.txt rfam-regions.txt ncRNA-regions.blocklist ncRNA-regions Rfam.cm


# download and extract spurious ORF HMMs from AntiFam
printf "\n4/18: download and extract spurious ORF HMMs from AntiFam ...\n"
mkdir antifam-dir
cd antifam-dir
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/AntiFam/current/Antifam.tar.gz
tar -xzf Antifam.tar.gz
cd ..
mv antifam-dir/AntiFam_Bacteria.hmm antifam
hmmpress antifam
rm -r antifam antifam-dir/


# download & extract oriT sequences
printf "\n5/18: download and extract oriT sequences from Mob-suite ...\n"
wget https://zenodo.org/record/3786915/files/data.tar.gz
tar -xvzf data.tar.gz
mv data/orit.fas ./orit.fna
rm -r data/ data.tar.gz
printf "\n5/18: download and extract oriC sequences from DoriC ...\n"
wget http://tubic.org/doric/public/static/download/doric10.rar
unrar e doric10.rar
python3 ${BAKTA_DB_SCRIPTS}/extract-ori.py --doric tubic_bacteria.csv --fasta oric.chromosome.fna
python3 ${BAKTA_DB_SCRIPTS}/extract-ori.py --doric tubic_plasmid.csv --fasta oric.plasmid.fna
cat oric.plasmid.fna >> oric.chromosome.fna
mv oric.chromosome.fna oric.fna
rm doric10.rar tubic* oric.plasmid.fna


# download NCBI Taxonomy DB
printf "\n6/18: download NCBI Taxonomy DB ...\n"
mkdir taxonomy
cd taxonomy
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -I pigz -xf taxdump.tar.gz
cd ..
mv taxonomy/nodes.dmp .
rm -rf taxonomy


############################################################################
# Setup SQLite Bakta db
############################################################################
printf "\n7/18: setup SQLite Bakta db ...\n"
python3 ${BAKTA_DB_SCRIPTS}/init-db.py --db bakta.db


############################################################################
# Build protein sequence clusters (PSCs) based on UniRef90 entries
# - download UniProt UniRef90
# - read and transform UniRef90 XML file to DB and Fasta file
# - build PSC Diamond db
############################################################################
printf "\n8/18: download UniProt UniRef90 ...\n"
wget https://ftp.expasy.org/databases/uniprot/current_release/uniref/uniref90/uniref90.xml.gz
wget https://ftp.expasy.org/databases/uniprot/current_release/uniref/uniref50/uniref50.xml.gz
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/uniparc/uniparc_active.fasta.gz
printf "\n8/18: read UniRef90 entries and build Protein Sequence Cluster sequence and information databases:\n"
python3 ${BAKTA_DB_SCRIPTS}/init-psc.py --taxonomy nodes.dmp --uniref90 uniref90.xml.gz --uniref50 uniref50.xml.gz --uniparc uniparc_active.fasta.gz --db bakta.db --psc psc.faa --sorf sorf.faa
printf "\n8/18: build PSC Diamond db ...\n"
diamond makedb --in psc.faa --db psc
diamond makedb --in sorf.faa --db sorf
rm uniref90.xml.gz


############################################################################
# Build unique protein sequences (IPSs) based on UniRef100 entries
# - download UniProt UniRef100
# - read, filter and transform UniRef100 entries and store to ips.db
############################################################################
printf "\n9/18: download UniProt UniRef100 ...\n"
wget https://ftp.expasy.org/databases/uniprot/current_release/uniref/uniref100/uniref100.xml.gz
printf "\n9/18: read, filter and store UniRef100 entries ...:\n"
python3 ${BAKTA_DB_SCRIPTS}/init-ups-ips.py --taxonomy nodes.dmp --uniref100 uniref100.xml.gz --uniparc uniparc_active.fasta.gz --db bakta.db --ips ips.faa
rm uniref100.xml.gz uniparc_active.fasta.gz


############################################################################
# Integrate NCBI COG db
# - download NCBI COG db
# - align UniRef90 proteins to COG protein sequences
# - annotate PSCs with COG info
############################################################################
printf "\n10/18: download COG db ...\n"
wget https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.def.tab  # COG IDs and functional class
wget https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.cog.csv # Mapping GenBank IDs -> COG IDs
for i in $(seq -f "%04g" 1 5950)
do
    wget https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/fasta/COG${i}.fa.gz
    pigz -dc COG${i}.fa.gz | seqtk seq -CU >> cog.faa
    rm COG${i}.fa.gz
done
printf "\n10/18: annotate PSCs ...\n"
diamond makedb --in cog.faa --db cog.dmnd
nextflow run ${BAKTA_DB_SCRIPTS}/diamond.nf --in psc.faa --db cog.dmnd --block 1000000 --id 90 --qcov 80 --scov 80 --out diamond.cog.psc.tsv
python3 ${BAKTA_DB_SCRIPTS}/annotate-cog.py --db bakta.db --alignments diamond.cog.psc.tsv --cog-ids cog-20.def.tab --gi-cog-mapping cog-20.cog.csv
nextflow run ${BAKTA_DB_SCRIPTS}/diamond.nf --in sorf.faa --db cog.dmnd --block 1000000 --id 90 --qcov 90 --scov 90 --out diamond.cog.sorf.tsv
python3 ${BAKTA_DB_SCRIPTS}/annotate-cog.py --db bakta.db --alignments diamond.cog.sorf.tsv --cog-ids cog-20.def.tab --gi-cog-mapping cog-20.cog.csv
rm cognames2003-2015.tab cog2003-2015.csv prot2003-2015.fa.gz diamond.cog.tsv cog.*


############################################################################
# Integrate KEGG kofams
# - download KEGG kofams
# - select eligible HMMs (f measure>0.77)
# - annotate PSCs
############################################################################
printf "\n11/18: download KEGG kofams HMM models...\n"
wget https://www.genome.jp/ftp/db/kofam/ko_list.gz
wget https://www.genome.jp/ftp/db/kofam/profiles.tar.gz
zcat ko_list.gz | grep full | awk '{ if($5>=0.77) print $0}' > hmms.selected.tsv
cut -f1 hmms.selected.tsv > hmms.ids.txt
tar -I pigz -xf profiles.tar.gz
for kofam in `cat profiles/prokaryote.hal`; do cat profiles/$kofam >> kofam-prok; done
hmmfetch -f -o kofams kofam-prok hmms.ids.txt
hmmpress kofams
printf "\n11/18: annotate PSCs...\n"
mkdir -p work/collect
nextflow run ${BAKTA_DB_SCRIPTS}/hmmsearch.nf --in psc.faa --db kofams --no_tc --dom
python3 ${BAKTA_DB_SCRIPTS}/annotate-kofams.py --db bakta.db --hmms hmms.selected.tsv --hmm-results hmmsearch.tblout
rm -rf profiles ko_list.gz kofam* hmmsearch.tblout hmms*


############################################################################
# Integrate NCBI nonredundant protein identifiers and PCLA cluster information
# - download bacterial RefSeq nonredundant proteins and cluster files
# - annotate UPSs with NCBI nrp IDs (WP_*)
# - annotate IPSs/PSCs with NCBI gene names (WP_* -> hash -> UniRef100 -> UniRef90 -> PSC)
############################################################################
printf "\n12/18: download RefSeq nonredundant proteins and clusters ...\n"
wget https://ftp.ncbi.nlm.nih.gov/genomes/CLUSTERS/PCLA_proteins.txt
wget https://ftp.ncbi.nlm.nih.gov/genomes/CLUSTERS/PCLA_clusters.txt
for i in {1..1573}; do
    wget https://ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/bacteria.nonredundant_protein.${i}.protein.faa.gz
    pigz -dc bacteria.nonredundant_protein.${i}.protein.faa.gz | seqtk seq -CU >> refseq-bacteria-nrp.trimmed.faa
    rm bacteria.nonredundant_protein.${i}.protein.faa.gz
done
printf "\n12/18: annotate IPSs and PSCs ...\n"
python3 ${BAKTA_DB_SCRIPTS}/annotate-ncbi-nrp.py --db bakta.db --nrp refseq-bacteria-nrp.trimmed.faa --pcla-proteins PCLA_proteins.txt --pcla-clusters PCLA_clusters.txt
rm refseq-bacteria-nrp.trimmed.faa PCLA_proteins.txt PCLA_clusters.txt


############################################################################
# Integrate UniProt Swissprot information
# - download SwissProt annotation xml file
# - annotate PSCs if IPS have PSC UniRef90 identifier (seq -> hash -> UPS -> IPS -> PSC)
# - annotate IPSs if IPS have no PSC UniRef90 identifier (seq -> hash -> UPS -> IPS)
############################################################################
printf "\n13/18: download UniProt/SwissProt ...\n"
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.xml.gz
printf "\n13/18: annotate IPSs and PSCs ...\n"
python3 ${BAKTA_DB_SCRIPTS}/annotate-swissprot.py --taxonomy nodes.dmp --xml uniprot_sprot.xml.gz --db bakta.db
rm uniprot_sprot.xml.gz


############################################################################
# Integrate NCBIfams HMM models
# - download NCBIfams HMM models
# - annotate PSCs
############################################################################
printf "\n14/18: download NCBIfams HMM models...\n"
wget https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.LIB
wget https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.tsv
grep -v "(Provisional)" hmm_PGAP.tsv > hmms.non-prov.tsv
grep exception hmms.non-prov.tsv >> hmms.selected.tsv
grep equivalog hmms.non-prov.tsv >> hmms.selected.tsv
cut -f1 hmms.selected.tsv > hmms.ids.txt
hmmfetch -f -o ncbifams hmm_PGAP.LIB hmms.ids.txt
hmmpress ncbifams
printf "\n14/18: annotate PSCs...\n"
mkdir -p work/collect
nextflow run ${BAKTA_DB_SCRIPTS}/hmmsearch.nf --in psc.faa --db ncbifams
python3 ${BAKTA_DB_SCRIPTS}/annotate-ncbi-fams.py --db bakta.db --hmms hmms.selected.tsv --hmm-results hmmsearch.tblout
rm ncbifams* hmms.* hmm_PGAP.* hmmsearch.tblout


############################################################################
# Integrate NCBI Pathogen AMR db
# - download AMR gene WP_* annotations from NCBI Pathogen ReferenceGeneCatalog
# - annotate IPSs with AMR info
############################################################################
printf "\n15/18: download AMR gene WP_* annotations from NCBI Pathogen AMR db ...\n"
wget https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/ReferenceGeneCatalog.txt
printf "\n15/18: annotate PSCs...\n"
python3 ${BAKTA_DB_SCRIPTS}/annotate-ncbi-amr.py --db bakta.db --genes ReferenceGeneCatalog.txt
rm ReferenceGeneCatalog.txt


############################################################################
# Integrate ISfinder db
# - download IS protein sequences from GitHub (oschwengers/ISfinder-sequences)
# - annotate IPSs with IS info
############################################################################
printf "\n16/18: download ISfinder protein sequences ...\n"
wget https://github.com/oschwengers/ISfinder-sequences/raw/2e9162bd5e3448c86ec1549a55315e498bef72fc/IS.faa
printf "\n16/18: annotate IPSs ...\n"
grep -A 1 ~~~Transposase~~~ IS.faa | tr -d - | tr -s "\n" > is.transposase.faa
diamond makedb --in is.transposase.faa --db is
nextflow run ${BAKTA_DB_SCRIPTS}/diamond.nf --in ips.faa --db is.dmnd --block 1000000 --id 98 --qcov 99 --scov 99 --out diamond.is.ips.tsv
nextflow run ${BAKTA_DB_SCRIPTS}/diamond.nf --in psc.faa --db is.dmnd --block 1000000 --id 90 --qcov 80 --scov 80 --out diamond.is.psc.tsv
python3 ${BAKTA_DB_SCRIPTS}/annotate-is.py --db bakta.db --ips-alignments diamond.is.ips.tsv --psc-alignments diamond.is.psc.tsv
rm IS.faa is.transposase.faa is.dmnd diamond.is.ips.tsv diamond.is.psc.tsv


############################################################################
# Integrate Pfam A
# - download all Pfam A HMM models
# - extract families & domains
# - compress HMM models
# - annotate hypothetical PSC via Pfam families
############################################################################
printf "\n17/18: download HMM models from Pfam ...\n"
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz
python3 ${BAKTA_DB_SCRIPTS}/extract-pfam.py --pfam Pfam-A.hmm.dat.gz --family pfam.families.tsv --non-family pfam.non-families.tsv
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
pigz -d Pfam-A.hmm.gz
hmmfetch -o pfam-families -f Pfam-A.hmm pfam.families.tsv
hmmpress pfam-families
hmmfetch -o pfam -f Pfam-A.hmm pfam.non-families.tsv
hmmpress pfam
python3 ${BAKTA_DB_SCRIPTS}/extract-hypotheticals.py --psc psc.faa --db bakta.db --hypotheticals hypotheticals.faa
mkdir -p work/collect
nextflow run ${BAKTA_DB_SCRIPTS}/hmmsearch.nf --in hypotheticals.faa --db pfam-families
python3 ${BAKTA_DB_SCRIPTS}/annotate-pfam.py --db bakta.db --hmms pfam-families --hmm-results hmmsearch.tblout
rm pfam-families* pfam *.tsv Pfam* hmmsearch.tblout


############################################################################
# Setup expert protein sequences
# - import NCBI BlastRules models
# - import VFDB sequences
############################################################################
printf "\n18/18: download AA sequences for expert annotation system ...\n"
wget https://ftp.ncbi.nlm.nih.gov/pub/blastrules/4.2.2.tgz
tar -xzf 4.2.2.tgz
python3 ${BAKTA_DB_SCRIPTS}/expert/setup-ncbiblastrules.py --expert-sequence expert-protein-sequences.faa --ncbi-blastrule-tsv 4.2.2/data/blast-rules_4.2.2.tsv --proteins 4.2.2/data/proteins.fasta
wget http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz
gunzip VFDB_setA_pro.fas.gz
python3 ${BAKTA_DB_SCRIPTS}/expert/setup-vfdb.py --expert-sequence expert-protein-sequences.faa --proteins VFDB_setA_pro.fas
diamond makedb --in expert-protein-sequences.faa --db expert-protein-sequences
rm -r 4.2.2/ 4.2.2.tgz VFDB_setA_pro.fas expert-protein-sequences.faa

# Cleanup
ls -l bakta.db
python3 ${BAKTA_DB_SCRIPTS}/optimize-db.py --db bakta.db
ls -l bakta.db
rm psc.faa sorf.faa node.dmp
