#!/bin/bash


THREADS="1"

while getopts "t" options; do
    case "${options}" in
        t)
            if [[ "$OPTARG" == -* ]]; then
                echo "Error: -p requires an argument!"
                exit 1
            else
                THREADS=${OPTARG}
            fi
            ;;
        esac
    done

mkdir db
cd db

printf "Create Bakta database\n"

# download rRNA covariance models from Rfam
printf "\n1/X: download rRNA covariance models from Rfam:\n"
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
printf "\n2/X: download ncRNA covariance models from Rfam:\n"
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
printf "\n3/X: download NCBI Taxonomy DB:\n"
mkdir taxonomy
cd taxonomy
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -I pigz -xf taxdump.tar.gz
cd ..
mv taxdump/nodes.dmp .
rm -rf taxdump


# download UniProt UniRef100
############################################################################
# Build unique protein sequences (UPSs) based on UniRef100 entries
# - download UniProt UniRef100
# - read, filter and transform UniRef100 entries and store to ups.db
############################################################################
printf "\n4/X: download UniProt UniRef100:\n"
wget ftp://ftp.expasy.org/databases/uniprot/current_release/uniref/uniref100/uniref100.xml.gz
printf "\n4/X: read, filter and store UniRef100 entries:\n"
python3 init-ups.py --taxonomy node.dmp --xml uniref100.xml.gz --db ups.db
rm uniref100.xml.gz


############################################################################
# Build protein sequence clusters (PSCs) based on UniRef90 entries
# - download UniProt UniRef90
# - read and transform UniRef90 XML file to DB and Fasta file
# - build PSC GhostZ db
############################################################################
printf "\n5/X: download UniProt UniRef90:\n"
wget ftp://ftp.expasy.org/databases/uniprot/current_release/uniref/uniref90/uniref90.xml.gz
printf "\n5/X: read UniRef90 entries and build Protein Sequence Cluster sequence and information databases:\n"
python3 init-psc.py --taxonomy node.dmp --xml uniref90.xml.gz --db psc.db --fasta psc.faa
printf "\n5/X: build PSC GhostZ db...\n"
ghostz db -i psc.faa -o psc -L 6
rm uniref90.xml.gz node.dmp psc.faa


############################################################################
# Integrate NCBI nonredundant protein identifiers and PCLA cluster information
# - download bacterial RefSeq nonredundant proteins and cluster files
# - annotate UPSs with NCBI nrp IDs (WP_*)
# - annotate PSCs with NCBI gene names (WP_* -> hash -> UniRef100 -> UniRef90 -> PSC)
############################################################################
printf "\n6/10: download RefSeq nonredundant proteins and clusters:\n"
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/CLUSTERS/PCLA_proteins.txt
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/CLUSTERS/PCLA_clusters.txt
mkdir refseq-bacteria
cd refseq-bacteria
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/bacteria.nonredundant_protein.*.protein.faa.gz
cd ..
pigz -dc refseq-bacteria/bacteria.nonredundant_protein.* | seqtk seq -CU > refseq-bacteria-nrp.trimmed.faa
rm -r refseq-bacteria
printf "\n6/10: annotate UPSs and PSCs:\n"
python3 annotate-ups-psc-ncbi-nrp.py --db-ups ups.db --db-psc psc.db --nrp refseq-bacteria-nrp.trimmed.faa --pcla-proteins PCLA_proteins.txt  --pcla-clusters PCLA_clusters.txt
rm refseq-bacteria-nrp.trimmed.faa


# download AMR gene WP_* annotations from NCBI Pathogen AMR db
printf "\n7/X: download AMR gene WP_* annotations from NCBI Pathogen AMR db...\n"
wget ftp://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/Data/latest/ReferenceGeneCatalog.txt
python3 annotate-ups-ncbi-amr.py --db ups.db --amr ReferenceGeneCatalog.txt
rm ReferenceGeneCatalog.txt


# download COG db and annotate PSCs
printf "\n8/X: download COG db and annotate PSCs...\n"
wget ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/cognames2003-2014.tab  # COG IDs and functional class
wget ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/prot2003-2014.tab  # GeneID -> RefSeq protein ID (old)
wget ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/prot2003-2014.fa.gz  # annotated protein sequences
pigz -dc prot2003-2014.fa.gz | seqtk seq -CU > cog.faa
ghostz db -i cog.faa -o cog
ghostz aln -i psc.faa -d cog -o ghostz.tsv -b 1 -a ${THREADS}
python3 annotate-psc-cog.py --db psc.db --cog-ids cognames2003-2014.tab --id-mapping prot2003-2014.tab
rm cognames2003-2014.tab prot2003-2014.tab prot2003-2014.fa.gz ghostz.tsv cog.*
