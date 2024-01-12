#!/bin/bash

set -e

mkdir db
cd db

printf "Create Bakta database\n"

# download rRNA covariance models from Rfam
printf "\n1/20: download rRNA covariance models from Rfam ...\n"
wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
pigz -d Rfam.cm.gz
cmfetch Rfam.cm RF00001 >  rRNA
cmfetch Rfam.cm RF00177 >> rRNA
cmfetch Rfam.cm RF02541 >> rRNA
cmpress rRNA
rm rRNA


# download and extract ncRNA gene covariance models from Rfam
printf "\n2/20: download ncRNA gene covariance models from Rfam ...\n"
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
printf "\n3/20: download ncRNA region covariance models from Rfam ...\n"
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
printf "\n4/20: download and extract spurious ORF HMMs from AntiFam ...\n"
mkdir antifam-dir
cd antifam-dir
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/AntiFam/current/Antifam.tar.gz
tar -xzf Antifam.tar.gz
cd ..
mv antifam-dir/AntiFam_Bacteria.hmm antifam
hmmpress antifam
rm -r antifam antifam-dir/


# download & extract oriT sequences
printf "\n5/20: download and extract oriT sequences from Mob-suite ...\n"
wget https://zenodo.org/records/10304948/files/data.tar.gz
tar -xvzf data.tar.gz
mv data/orit.fas ./orit.fna
rm -r data/ data.tar.gz
printf "\n5/20: download oriC/V sequences from DoriC ...\n"
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


# download NCBI Taxonomy DB
printf "\n6/20: download NCBI Taxonomy DB ...\n"
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
printf "\n7/20: setup SQLite Bakta db ...\n"
python3 ${BAKTA_DB_SCRIPTS}/init-db.py --db bakta.db


############################################################################
# Build protein sequence clusters (PSCCs) based on UniRef50 entries
# - download UniProt UniRef50
# - read and transform UniRef50 XML file to DB and Fasta file
# - build PSCC Diamond db
############################################################################
printf "\n8/20: download UniProt UniRef50 ...\n"
wget https://ftp.expasy.org/databases/uniprot/current_release/uniref/uniref50/uniref50.xml.gz
for i in {1..200}; do
    wget https://ftp.expasy.org/databases/uniprot/current_release/uniparc/fasta/active/uniparc_active_p${i}.fasta.gz
    pigz -dc uniparc_active_p${i}.fasta.gz >> uniparc_active.fasta
    rm uniparc_active_p${i}.fasta.gz
done
printf "\n8/20: read UniRef90 entries and build Protein Sequence Cluster sequence and information databases:\n"
python3 ${BAKTA_DB_SCRIPTS}/init-pscc.py --taxonomy nodes.dmp --uniref50 uniref50.xml.gz --uniparc uniparc_active.fasta --db bakta.db --pscc pscc.faa --sorf pscc_sorf.faa
printf "\n8/20: build PSCC Diamond db ...\n"
diamond makedb --in pscc.faa --db pscc
diamond makedb --in pscc_sorf.faa --db sorf
mkdir db-lite
cp bakta.db db-lite
mv pscc.dmnd sorf.dmnd db-lite
cd db-lite
python3 ${BAKTA_DB_SCRIPTS}/optimize-db.py --db bakta.db
cd ..
rm uniref50.xml.gz


############################################################################
# Build protein sequence clusters (PSCs) based on UniRef90 entries
# - download UniProt UniRef90
# - read and transform UniRef90 XML file to DB and Fasta file
# - build PSC Diamond db
############################################################################
printf "\n9/20: download UniProt UniRef90 ...\n"
wget https://ftp.expasy.org/databases/uniprot/current_release/uniref/uniref90/uniref90.xml.gz
printf "\n9/20: read UniRef90 entries and build Protein Sequence Cluster sequence and information databases:\n"
python3 ${BAKTA_DB_SCRIPTS}/init-psc.py --taxonomy nodes.dmp --uniref90 uniref90.xml.gz --uniparc uniparc_active.fasta --db bakta.db --psc psc.faa --sorf sorf.faa
printf "\n9/20: build PSC Diamond db ...\n"
diamond makedb --in psc.faa --db psc
diamond makedb --in sorf.faa --db sorf
rm uniref90.xml.gz


############################################################################
# Build unique protein sequences (IPSs) based on UniRef100 entries
# - download UniProt UniRef100
# - read, filter and transform UniRef100 entries and store to ips.db
############################################################################
printf "\n10/20: download UniProt UniRef100 ...\n"
wget https://ftp.expasy.org/databases/uniprot/current_release/uniref/uniref100/uniref100.xml.gz
printf "\n10/20: read, filter and store UniRef100 entries ...:\n"
python3 ${BAKTA_DB_SCRIPTS}/init-ups-ips.py --taxonomy nodes.dmp --uniref100 uniref100.xml.gz --uniparc uniparc_active.fasta.gz --db bakta.db --ips ips.faa
rm uniref100.xml.gz uniparc_active.fasta.gz


############################################################################
# Integrate NCBI COG db
# - download NCBI COG db
# - align UniRef90 proteins to COG protein sequences
# - annotate PSCs with COG info
############################################################################
printf "\n11/20: download COG db ...\n"
wget https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.def.tab  # COG IDs and functional class
wget https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.cog.csv # Mapping GenBank IDs -> COG IDs
for i in $(seq -f "%04g" 1 5950)
do
    wget https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/fasta/COG${i}.fa.gz
    pigz -dc COG${i}.fa.gz | seqtk seq -CU >> cog.faa
    rm COG${i}.fa.gz
done
printf "\n11/20: annotate PSCs ...\n"
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
printf "\n12/20: download KEGG kofams HMM models...\n"
wget https://www.genome.jp/ftp/db/kofam/ko_list.gz
wget https://www.genome.jp/ftp/db/kofam/profiles.tar.gz
zcat ko_list.gz | grep full | awk '{ if($5>=0.77) print $0}' > hmms.kofam.selected.tsv
cut -f1 hmms.kofam.selected.tsv > hmms.ids.txt
tar -I pigz -xf profiles.tar.gz
for kofam in `cat profiles/prokaryote.hal`; do cat profiles/$kofam >> kofam-prok; done
hmmfetch -f -o kofams kofam-prok hmms.ids.txt
hmmpress kofams
printf "\n12/20: annotate PSCs...\n"
mkdir -p work/tblout work/domtblout
nextflow run ${BAKTA_DB_SCRIPTS}/hmmsearch.nf --in psc.faa --db kofams --no_tc --out hmmsearch.kofam.tblout
python3 ${BAKTA_DB_SCRIPTS}/annotate-kofams.py --db bakta.db --hmms hmms.kofam.selected.tsv --hmm-results hmmsearch.kofam.tblout
rm -rf profiles ko_list.gz kofam* hmmsearch.kofam.* hmms*


############################################################################
# Integrate NCBI nonredundant protein identifiers and PCLA cluster information
# - download bacterial RefSeq nonredundant proteins and cluster files
# - annotate UPSs with NCBI nrp IDs (WP_*)
# - annotate IPSs/PSCs with NCBI gene names (WP_* -> hash -> UniRef100 -> UniRef90 -> PSC)
############################################################################
printf "\n13/20: download RefSeq nonredundant proteins and clusters ...\n"
wget https://ftp.ncbi.nlm.nih.gov/genomes/CLUSTERS/PCLA_proteins.txt
wget https://ftp.ncbi.nlm.nih.gov/genomes/CLUSTERS/PCLA_clusters.txt
for i in {1..360}; do
    wget https://ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/bacteria.nonredundant_protein.${i}.protein.faa.gz
    pigz -dc bacteria.nonredundant_protein.${i}.protein.faa.gz | seqtk seq -CU >> refseq-bacteria-nrp.trimmed.faa
    rm bacteria.nonredundant_protein.${i}.protein.faa.gz
done
printf "\n13/20: annotate IPSs and PSCs ...\n"
python3 ${BAKTA_DB_SCRIPTS}/annotate-ncbi-nrp.py --db bakta.db --nrp refseq-bacteria-nrp.trimmed.faa --pcla-proteins PCLA_proteins.txt --pcla-clusters PCLA_clusters.txt
rm refseq-bacteria-nrp.trimmed.faa PCLA_proteins.txt PCLA_clusters.txt


############################################################################
# Integrate UniProt Swissprot information
# - download SwissProt annotation xml file
# - annotate PSCs if IPS have PSC UniRef90 identifier (seq -> hash -> UPS -> IPS -> PSC)
# - annotate IPSs if IPS have no PSC UniRef90 identifier (seq -> hash -> UPS -> IPS)
############################################################################
printf "\n14/20: download UniProt/SwissProt ...\n"
wget https://ftp.expasy.org/databases/swiss-prot/release/uniprot_sprot.xml.gz
printf "\n14/20: annotate IPSs and PSCs ...\n"
python3 ${BAKTA_DB_SCRIPTS}/annotate-swissprot.py --taxonomy nodes.dmp --xml uniprot_sprot.xml.gz --db bakta.db
rm uniprot_sprot.xml.gz


############################################################################
# Integrate NCBIfams HMM models
# - download NCBIfams HMM models
# - annotate PSCs
############################################################################
printf "\n15/20: download NCBIfams HMM models...\n"
wget https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.LIB
wget https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.tsv
grep -v "(Provisional)" hmm_PGAP.tsv > hmms.non-prov.tsv
grep exception hmms.non-prov.tsv > hmms.ncbi.selected.tsv
grep equivalog hmms.non-prov.tsv >> hmms.ncbi.selected.tsv
cut -f1 hmms.ncbi.selected.tsv > hmms.ids.txt
hmmfetch -f -o ncbifams hmm_PGAP.LIB hmms.ids.txt
hmmpress ncbifams
printf "\n15/20: annotate PSCs...\n"
mkdir -p work/tblout work/domtblout
nextflow run ${BAKTA_DB_SCRIPTS}/hmmsearch.nf --in psc.faa --db ncbifams --out hmmsearch.ncbifams.tblout
python3 ${BAKTA_DB_SCRIPTS}/annotate-ncbi-fams.py --db bakta.db --hmms hmms.ncbi.selected.tsv --hmm-results hmmsearch.ncbifams.tblout
rm ncbifams* hmms.* hmm_PGAP.* hmmsearch.ncbifams.tblout


############################################################################
# Integrate PHROG DB of phage orthologous
# - download PHROG protein sequences
# - filter unannotated PHROGs
# - annotate PSCs
############################################################################
printf "\n16/20: download PHROGs ...\n"
wget https://phrogs.lmge.uca.fr/downloads_from_website/FAA_phrog.tar.gz
wget https://phrogs.lmge.uca.fr/downloads_from_website/phrog_annot_v4.tsv
tar -xzf FAA_phrog.tar.gz
cat FAA_phrog/*.faa >> phrogs-raw.faa
python3 ${BAKTA_DB_SCRIPTS}/extract-phrogs.py --annotation phrog_annot_v4.tsv --proteins phrogs-raw.faa --filtered-proteins phrogs.faa
diamond makedb --in phrogs.faa --db phrog
printf "\n16/20: annotate PSCs...\n"
python3 ${BAKTA_DB_SCRIPTS}/extract-hypotheticals.py --psc psc.faa --db bakta.db --hypotheticals hypotheticals.faa
nextflow run ${BAKTA_DB_SCRIPTS}/diamond.nf --in hypotheticals.faa --db phrog.dmnd --block 1000000 --id 90 --qcov 80 --scov 80 --out diamond.phrog.psc.tsv
python3 ${BAKTA_DB_SCRIPTS}/annotate-phrogs.py --db bakta.db --annotation phrog_annot_v4.tsv --psc-alignments diamond.phrog.psc.tsv
rm -r FAA_phrog.tar.gz phrog_annot_v4.tsv FAA_phrog phrogs-raw.faa phrogs.faa phrog.dmnd hypotheticals.faa


############################################################################
# Integrate NCBI Pathogen AMR db
# - download AMR gene WP_* annotations from NCBI Pathogen ReferenceGeneCatalog
# - annotate IPSs with AMR info
############################################################################
printf "\n17/20: download AMR gene WP_* annotations from NCBI Pathogen AMR db ...\n"
wget https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/ReferenceGeneCatalog.txt
printf "\n17/20: annotate PSCs...\n"
python3 ${BAKTA_DB_SCRIPTS}/annotate-ncbi-amr.py --db bakta.db --genes ReferenceGeneCatalog.txt
rm ReferenceGeneCatalog.txt


############################################################################
# Integrate ISfinder db
# - download IS protein sequences from GitHub (oschwengers/ISfinder-sequences)
# - extract IS transposase sequences and mark ORF A/B transposases
# - annotate IPSs/PCSs with IS info
############################################################################
printf "\n18/20: download & extract ISfinder protein sequences ...\n"
wget https://github.com/oschwengers/ISfinder-sequences/raw/2e9162bd5e3448c86ec1549a55315e498bef72fc/IS.faa
python3 ${BAKTA_DB_SCRIPTS}/extract-is.py --input IS.faa --output is.transposase.faa
printf "\n18/20: annotate IPSs/PCSs ...\n"
diamond makedb --in is.transposase.faa --db is
nextflow run ${BAKTA_DB_SCRIPTS}/diamond.nf --in ips.faa --db is.dmnd --block 1000000 --id 95 --qcov 90 --scov 90 --out diamond.is.ips.tsv
nextflow run ${BAKTA_DB_SCRIPTS}/diamond.nf --in psc.faa --db is.dmnd --block 1000000 --id 90 --qcov 80 --scov 80 --out diamond.is.psc.tsv
python3 ${BAKTA_DB_SCRIPTS}/annotate-is.py --db bakta.db --ips-alignments diamond.is.ips.tsv --psc-alignments diamond.is.psc.tsv
rm is.transposase.faa is.dmnd diamond.is.ips.tsv diamond.is.psc.tsv


############################################################################
# Integrate Pfam A
# - download all Pfam A HMM models
# - extract families & domains
# - compress HMM models
# - annotate hypothetical PSC via Pfam families
############################################################################
printf "\n19/20: download HMM models from Pfam ...\n"
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz
python3 ${BAKTA_DB_SCRIPTS}/extract-pfam.py --pfam Pfam-A.hmm.dat.gz --family pfam.families.tsv --non-family pfam.non-families.tsv
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
pigz -d Pfam-A.hmm.gz
hmmfetch -o pfam-families -f Pfam-A.hmm pfam.families.tsv
hmmpress pfam-families
hmmfetch -o pfam -f Pfam-A.hmm pfam.non-families.tsv
hmmpress pfam
python3 ${BAKTA_DB_SCRIPTS}/extract-hypotheticals.py --psc psc.faa --db bakta.db --hypotheticals hypotheticals.faa
mkdir -p work/tblout work/domtblout
nextflow run ${BAKTA_DB_SCRIPTS}/hmmsearch.nf --in hypotheticals.faa --db pfam-families --out hmmsearch.pfam-families.tblout
python3 ${BAKTA_DB_SCRIPTS}/annotate-pfam.py --db bakta.db --hmms pfam-families --hmm-results hmmsearch.pfam-families.tblout
rm pfam-families* pfam *.tsv Pfam* hmmsearch.pfam-families.tblout


############################################################################
# Setup expert protein sequences
# - import IS sequences
# - import NCBI BlastRules models
# - import VFDB sequences
############################################################################
printf "\n20/20: download AA sequences for expert annotation system ...\n"
wget https://ftp.ncbi.nlm.nih.gov/pub/blastrules/4.2.2.tgz
tar -xzf 4.2.2.tgz
wget http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz
gunzip VFDB_setA_pro.fas.gz
python3 ${BAKTA_DB_SCRIPTS}/expert/setup-is.py --expert-sequence expert-protein-sequences.faa --proteins IS.faa
python3 ${BAKTA_DB_SCRIPTS}/expert/setup-ncbiblastrules.py --expert-sequence expert-protein-sequences.faa --ncbi-blastrule-tsv 4.2.2/data/blast-rules_4.2.2.tsv --proteins 4.2.2/data/proteins.fasta
python3 ${BAKTA_DB_SCRIPTS}/expert/setup-vfdb.py --expert-sequence expert-protein-sequences.faa --proteins VFDB_setA_pro.fas
diamond makedb --in expert-protein-sequences.faa --db expert-protein-sequences
rm -r 4.2.2/ 4.2.2.tgz IS.faa VFDB_setA_pro.fas expert-protein-sequences.faa

# Cleanup
ls -l bakta.db
python3 ${BAKTA_DB_SCRIPTS}/optimize-db.py --db bakta.db
ls -l bakta.db
rm psc.faa sorf.faa node.dmp
