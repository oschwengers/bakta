#!/bin/bash

mkdir db
cd db

printf "Create Bakta database\n"

# download rRNA covariance models from Rfam
printf "\n1/X: download rRNA covariance models from Rfam ...\n"
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


# download and extract ncRNA covariance models from Rfam
printf "\n2/X: download ncRNA covariance models from Rfam ...\n"
mysql --user rfamro --host mysql-rfam-public.ebi.ac.uk --port 4497 --database Rfam < ncRNA.sql > rfam.raw.txt
tail -n +2 rfam.raw.txt > rfam.txt
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
pigz -d Rfam.cm.gz
cmfetch -o ncRNA -f Rfam.cm rfam.txt
cmpress ncRNA
wget http://current.geneontology.org/ontology/external2go/rfam2go
awk -F ' ' '{print $1 "\t" $NF}' rfam2go > rfam-go.tsv
rm rfam.raw.txt rfam.txt Rfam.cm rfam2go


# download NCBI Taxonomy DB
printf "\n3/X: download NCBI Taxonomy DB ...\n"
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
printf "\n4/X: setup SQLite Bakta db ...\n"
python3 ${BAKTA_DB_SCRIPTS}/init-db.py --db bakta.db


############################################################################
# Build unique protein sequences (UPSs) based on UniRef100 entries
# - download UniProt UniRef100
# - read, filter and transform UniRef100 entries and store to ups.db
############################################################################
printf "\n5/X: download UniProt UniRef100 ...\n"
wget ftp://ftp.expasy.org/databases/uniprot/current_release/uniref/uniref100/uniref100.xml.gz
printf "\n5/X: read, filter and store UniRef100 entrie ...:\n"
python3 ${BAKTA_DB_SCRIPTS}/init-ups.py --taxonomy nodes.dmp --xml uniref100.xml.gz --db bakta.db
rm uniref100.xml.gz


############################################################################
# Build protein sequence clusters (PSCs) based on UniRef90 entries
# - download UniProt UniRef90
# - read and transform UniRef90 XML file to DB and Fasta file
# - build PSC Diamond db
############################################################################
printf "\n6/X: download UniProt UniRef90 ...\n"
wget ftp://ftp.expasy.org/databases/uniprot/current_release/uniref/uniref90/uniref90.xml.gz
printf "\n6/X: read UniRef90 entries and build Protein Sequence Cluster sequence and information databases:\n"
python3 ${BAKTA_DB_SCRIPTS}/init-psc.py --taxonomy nodes.dmp --xml uniref90.xml.gz --db bakta.db --fasta psc.faa
printf "\n5/X: build PSC Diamond db ...\n"
diamond makedb --in psc.faa --db psc
rm uniref90.xml.gz node.dmp psc.faa


############################################################################
# Integrate NCBI nonredundant protein identifiers and PCLA cluster information
# - download bacterial RefSeq nonredundant proteins and cluster files
# - annotate UPSs with NCBI nrp IDs (WP_*)
# - annotate PSCs with NCBI gene names (WP_* -> hash -> UniRef100 -> UniRef90 -> PSC)
############################################################################
printf "\n7/10: download RefSeq nonredundant proteins and clusters ...\n"
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/CLUSTERS/PCLA_proteins.txt
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/CLUSTERS/PCLA_clusters.txt
mkdir refseq-bacteria
cd refseq-bacteria
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/bacteria.nonredundant_protein.*.protein.faa.gz
cd ..
pigz -dc refseq-bacteria/bacteria.nonredundant_protein.* | seqtk seq -CU > refseq-bacteria-nrp.trimmed.faa
rm -r refseq-bacteria
printf "\n7/10: annotate UPSs and PSCs ...\n"
python3 ${BAKTA_DB_SCRIPTS}/annotate-ncbi-nrp.py --db bakta.db --nrp refseq-bacteria-nrp.trimmed.faa --pcla-proteins PCLA_proteins.txt --pcla-clusters PCLA_clusters.txt
rm refseq-bacteria-nrp.trimmed.faa PCLA_proteins.txt PCLA_clusters.txt


############################################################################
# Integrate NCBI COG db
# - download NCBI COG db
# - align UniRef90 proteins to COG protein sequences
# - annotate PSCs with COG info
############################################################################
printf "\n8/X: download COG db ...\n"
wget ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/cognames2003-2014.tab  # COG IDs and functional class
wget ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/cog2003-2014.csv # Mapping GenBank IDs -> COG IDs
wget ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/prot2003-2014.fa.gz  # annotated protein sequences
pigz -dc prot2003-2014.fa.gz | seqtk seq -CU > cog.faa
printf "\n8/X: annotate PSCs ...\n"
diamond makedb --in cog.faa --db cog.dmnd
diamond blastp --query psc.faa --db cog.dmnd --id 90 --query-cover 80 --subject-cover 80 --max-target-seqs 1 --out diamond.tsv --outfmt 6 qseqid sseqid length pident qlen slen evalue
python3 ${BAKTA_DB_SCRIPTS}/annotate-cog.py --db bakta.db --alignments diamond.tsv --cog-ids cognames2003-2014.tab --id-mapping cog2003-2014.csv
rm cognames2003-2014.tab cog2003-2014.csv prot2003-2014.fa.gz diamond.tsv cog.*


############################################################################
# Integrate NCBI Pathogen AMR db
# - download AMR gene WP_* annotations from NCBI Pathogen AMR db
# - annotate UPSs with AMR info
############################################################################
printf "\n9/X: download AMR gene WP_* annotations from NCBI Pathogen AMR db ...\n"
wget ftp://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/Data/latest/ReferenceGeneCatalog.txt
printf "\n9/10: annotate UPSs ...\n"
python3 ${BAKTA_DB_SCRIPTS}/annotate-ncbi-amr.py --db bakta.db --amr ReferenceGeneCatalog.txt
rm ReferenceGeneCatalog.txt


############################################################################
# Integrate ISfinder db
# - download IS protein sequences from GitHub (thanhleviet/ISfinder-sequences)
# - annotate UPSs with IS info
############################################################################
printf "\n10/X: download ISfinder protein sequences ...\n"
wget https://raw.githubusercontent.com/thanhleviet/ISfinder-sequences/master/IS.faa
printf "\n10/X: annotate UPSs ...\n"
diamond makedb --in IS.faa --db is
diamond blastp --query psc.faa --db is.dmnd --id 98 --query-cover 99 --subject-cover 99 --max-target-seqs 1 --out diamond.tsv --outfmt 6 qseqid sseqid stitle length pident qlen slen evalue
python3 ${BAKTA_DB_SCRIPTS}/annotate-is.py --db bakta.db --alignments diamond.tsv
rm IS.faa is.dmnd


# Cleanup
rm psc.faa
