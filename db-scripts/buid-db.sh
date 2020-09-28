#!/bin/bash

set -e

mkdir db
cd db

printf "Create Bakta database\n"

# download rRNA covariance models from Rfam
printf "\n1/11: download rRNA covariance models from Rfam ...\n"
mkdir Rfam
cd Rfam
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.tar.gz
tar -I pigz -xf Rfam.tar.gz
cd ..
cat Rfam/RF00001.cm >  rRNA
cat Rfam/RF00177.cm >> rRNA
cat Rfam/RF02541.cm >> rRNA
cmpress rRNA
rm -r Rfam rRNA


# download and extract ncRNA gene covariance models from Rfam
printf "\n2/11: download ncRNA gene covariance models from Rfam ...\n"
mysql --user rfamro --host mysql-rfam-public.ebi.ac.uk --port 4497 --database Rfam < ${BAKTA_DB_SCRIPTS}/ncRNA-genes.sql > rfam-genes.raw.txt
tail -n +2 rfam-genes.raw.txt > rfam-genes.txt
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
pigz -d Rfam.cm.gz
cmfetch -o ncRNA-genes -f Rfam.cm rfam-genes.txt
cmpress ncRNA-genes
wget http://current.geneontology.org/ontology/external2go/rfam2go
awk -F ' ' '{print $1 "\t" $NF}' rfam2go > rfam-go.tsv
rm rfam-genes.raw.txt rfam-genes.txt rfam2go ncRNA-genes


# download and extract ncRNA regions (cis reg elements) covariance models from Rfam
printf "\n2/11: download ncRNA region covariance models from Rfam ...\n"
mysql --user rfamro --host mysql-rfam-public.ebi.ac.uk --port 4497 --database Rfam < ${BAKTA_DB_SCRIPTS}/ncRNA-regions.sql > rfam-regions.raw.txt
tail -n +2 rfam-regions.raw.txt > rfam-regions.txt
cmfetch -o ncRNA-regions -f Rfam.cm rfam-regions.txt
cmpress ncRNA-regions
rm rfam-regions.raw.txt rfam-regions.txt Rfam.cm ncRNA-regions


# download NCBI Taxonomy DB
printf "\n3/11: download NCBI Taxonomy DB ...\n"
mkdir taxonomy
cd taxonomy
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -I pigz -xf taxdump.tar.gz
cd ..
mv taxonomy/nodes.dmp .
rm -rf taxonomy


############################################################################
# Setup SQLite Bakta db
############################################################################
printf "\n4/11: setup SQLite Bakta db ...\n"
python3 ${BAKTA_DB_SCRIPTS}/init-db.py --db bakta.db


############################################################################
# Build unique protein sequences (IPSs) based on UniRef100 entries
# - download UniProt UniRef100
# - read, filter and transform UniRef100 entries and store to ips.db
############################################################################
printf "\n5/11: download UniProt UniRef100 ...\n"
wget -nv ftp://ftp.expasy.org/databases/uniprot/current_release/uniref/uniref100/uniref100.xml.gz
wget -nv ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/uniparc/uniparc_active.fasta.gz
printf "\n5/11: read, filter and store UniRef100 entries ...:\n"
python3 ${BAKTA_DB_SCRIPTS}/init-ups-ips.py --taxonomy nodes.dmp --xml uniref100.xml.gz --uniparc uniparc_active.fasta.gz --db bakta.db
rm uniref100.xml.gz


############################################################################
# Build protein sequence clusters (PSCs) based on UniRef90 entries
# - download UniProt UniRef90
# - read and transform UniRef90 XML file to DB and Fasta file
# - build PSC Diamond db
############################################################################
printf "\n6/11: download UniProt UniRef90 ...\n"
wget -nv ftp://ftp.expasy.org/databases/uniprot/current_release/uniref/uniref90/uniref90.xml.gz
printf "\n6/11: read UniRef90 entries and build Protein Sequence Cluster sequence and information databases:\n"
python3 ${BAKTA_DB_SCRIPTS}/init-psc.py --taxonomy nodes.dmp --xml uniref90.xml.gz --uniparc uniparc_active.fasta.gz --db bakta.db --fasta psc.faa
printf "\n5/11: build PSC Diamond db ...\n"
diamond makedb --in psc.faa --db psc
rm uniref90.xml.gz uniparc_active.fasta.gz


############################################################################
# Integrate NCBI nonredundant protein identifiers and PCLA cluster information
# - download bacterial RefSeq nonredundant proteins and cluster files
# - annotate UPSs with NCBI nrp IDs (WP_*)
# - annotate IPSs/PSCs with NCBI gene names (WP_* -> hash -> UniRef100 -> UniRef90 -> PSC)
############################################################################
printf "\n7/11: download RefSeq nonredundant proteins and clusters ...\n"
wget -nv ftp://ftp.ncbi.nlm.nih.gov/genomes/CLUSTERS/PCLA_proteins.txt
wget -nv ftp://ftp.ncbi.nlm.nih.gov/genomes/CLUSTERS/PCLA_clusters.txt
for i in {1..1169}; do
    wget -nv ftp://ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/bacteria.nonredundant_protein.${i}.protein.faa.gz
    pigz -dc bacteria.nonredundant_protein.${i}.protein.faa.gz | seqtk seq -CU >> refseq-bacteria-nrp.trimmed.faa
    rm bacteria.nonredundant_protein.${i}.protein.faa.gz
done
printf "\n7/11: annotate IPSs and PSCs ...\n"
python3 ${BAKTA_DB_SCRIPTS}/annotate-ncbi-nrp.py --db bakta.db --nrp refseq-bacteria-nrp.trimmed.faa --pcla-proteins PCLA_proteins.txt --pcla-clusters PCLA_clusters.txt
rm refseq-bacteria-nrp.trimmed.faa PCLA_proteins.txt PCLA_clusters.txt


############################################################################
# Integrate UniProt Swissprot information
# - download SwissProt annotation xml file
# - annotate PSCs if IPS have PSC UniRef90 identifier (seq -> hash -> UPS -> IPS -> PSC)
# - annotate IPSs if IPS have no PSC UniRef90 identifier (seq -> hash -> UPS -> IPS)
############################################################################
# printf "\n8/11: download UniProt/SwissProt ...\n"
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.xml.gz
printf "\n8/11: annotate IPSs and PSCs ...\n"
python3 ${BAKTA_DB_SCRIPTS}/annotate-swissprot.py --taxonomy nodes.dmp --xml uniprot_sprot.xml.gz --db bakta.db
rm uniprot_sprot.xml.gz


############################################################################
# Integrate NCBI COG db
# - download NCBI COG db
# - align UniRef90 proteins to COG protein sequences
# - annotate PSCs with COG info
############################################################################
printf "\n9/11: download COG db ...\n"
wget -nv ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/cognames2003-2014.tab  # COG IDs and functional class
wget -nv ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/cog2003-2014.csv # Mapping GenBank IDs -> COG IDs
wget -nv ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/prot2003-2014.fa.gz  # annotated protein sequences
pigz -dc prot2003-2014.fa.gz | seqtk seq -CU > cog.faa
printf "\n9/11: annotate PSCs ...\n"
diamond makedb --in cog.faa --db cog.dmnd
diamond blastp --query psc.faa --db cog.dmnd --id 90 --query-cover 80 --subject-cover 80 --max-target-seqs 1 --out diamond.tsv --outfmt 6 qseqid sseqid length pident qlen slen evalue
python3 ${BAKTA_DB_SCRIPTS}/annotate-cog.py --db bakta.db --alignments diamond.tsv --cog-ids cognames2003-2014.tab --gi-cog-mapping cog2003-2014.csv
rm cognames2003-2014.tab cog2003-2014.csv prot2003-2014.fa.gz diamond.tsv cog.*


############################################################################
# Integrate NCBI Pathogen AMR db
# - download AMR gene WP_* annotations from NCBI Pathogen AMR db
# - annotate IPSs with AMR info
############################################################################
printf "\n10/11: download AMR gene WP_* annotations from NCBI Pathogen AMR db ...\n"
wget -nv ftp://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/Data/latest/ReferenceGeneCatalog.txt
wget -nv -nH ftp://ftp.ncbi.nlm.nih.gov/hmm/NCBIfam-AMRFinder/latest/NCBIfam-AMRFinder.LIB
wget -nv -nH ftp://ftp.ncbi.nlm.nih.gov/hmm/NCBIfam-AMRFinder/latest/NCBIfam-AMRFinder.tsv
printf "\n10/11: annotate IPSs and PSCs...\n"
mv NCBIfam-AMRFinder.LIB ncbifam-amr
hmmpress ncbifam-amr
nextflow run ${BAKTA_DB_SCRIPTS}/annotate-ncbi-amr.nf --psc psc.faa --db ncbifam-amr
python3 ${BAKTA_DB_SCRIPTS}/annotate-ncbi-amr.py --db bakta.db --genes ReferenceGeneCatalog.txt --hmms NCBIfam-AMRFinder.tsv --hmm-results hmm-amr.out
rm ncbifam-amr* NCBIfam-AMRFinder.tsv ReferenceGeneCatalog.txt hmm-amr.out


############################################################################
# Integrate ISfinder db
# - download IS protein sequences from GitHub (thanhleviet/ISfinder-sequences)
# - annotate IPSs with IS info
############################################################################
printf "\n11/11: download ISfinder protein sequences ...\n"
wget -nv https://raw.githubusercontent.com/thanhleviet/ISfinder-sequences/master/IS.faa
printf "\n11/11: annotate IPSs ...\n"
diamond makedb --in IS.faa --db is
diamond blastp --query psc.faa --db is.dmnd --id 98 --query-cover 99 --subject-cover 99 --max-target-seqs 1 --out diamond.tsv --outfmt 6 qseqid sseqid stitle length pident qlen slen evalue
python3 ${BAKTA_DB_SCRIPTS}/annotate-is.py --db bakta.db --alignments diamond.tsv
rm IS.faa is.dmnd


# Cleanup
python3 ${BAKTA_DB_SCRIPTS}/optimize-db.py --db bakta.db
rm psc.faa node.dmp
