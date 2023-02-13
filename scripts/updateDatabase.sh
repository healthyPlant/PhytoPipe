#!/usr/bin/env bash

#************************************************************
#This script is for building and updating databases
#************************************************************
echo "please find RVDB version from https://rvdb-prot.pasteur.fr/, ex. v25.0"
#check argument number
if [[ $# -ne 3 ]]; then
    echo "Usage: bash /path/to/phytopipe/scripts/updateDatabase.sh [/path/to/phytopipe] [/path/my/database]  [rvdb_version, ex. v25.0]" >&2
    exit 2
fi

phytopipe_dir=$1
mydb=$2  #/ppq/data2/resources
rvdb_version=$3 #v25.0 #https://rvdb-prot.pasteur.fr/

ncbi=$mydb/ncbi
ncbi_nt=$mydb/ncbi_nt
krakendb=$mydb/kraken_db
kaijudb=$mydb/kaiju_db
taxondb=$mydb/ncbi/taxonomy
rrnadb=$mydb/rRNA
pv_taxonFile=$phytopipe_dir/db/plantvirus_taxon.txt
excludeSeqId=$phytopipe_dir/db/unwantedSeqId.txt

krakendb0=$mydb/kraken_db

if [ -d $krakendb ]; then 
    echo "$krakendb exists, we will build Kraken2 database under a new directory $krakendb.1"
    echo "After you test the new $krakendb.1, you can rename it by running \"mv $krakendb.1 $krakendb\" "
fi


echo "#*****************************"
echo "Check programs"
programs=("kraken2-build" "kaiju-makedb" "update_blastdb.pl" "diamond" "updateTaxonomy.sh" "updateAccessions.sh" "parallel" "blastdbcmd" "makeblastdb" "kaiju-addTaxonNames" "get_species_taxids.sh" "esearch" "filterbyname.sh")
for prog in ${programs[@]}
do
    if ! command -v $prog &> /dev/null
    then
        echo "$prog could not be found. Please intall it and/or add it in $PATH."
        exit
    else
        echo "$prog is ok"
    fi
done


#****************************************************
#Check database folders, if not exsit, create them
#my database directory
if [ ! -d $mydb ]; then 
    mkdir -p $mydb  
fi
#nr and reference sequence directory 
if [ ! -d $ncbi ]; then 
    mkdir -p $ncbi  
fi
#NCBI nt directory
if [ ! -d $ncbi_nt ]; then 
    mkdir -p $ncbi_nt  
fi
#Kraken2 db directory
if [ ! -d $krakendb ]; then 
    mkdir -p $krakendb 
else
    mkdir -p $krakendb.1
    krakendb=$krakendb.1
fi
#Kaiju db directory
if [ ! -d $kaijudb ]; then 
    mkdir -p $kaijudb 
fi
#NCBI taxonomy db directory 
if [ ! -d $taxondb ]; then 
    mkdir -p $taxondb  
fi
#rRNA db directory 
if [ ! -d $rrnadb ]; then 
    mkdir -p $rrnadb  
fi
#*****************************************
#Build or update databases
echo "#*****************************"
echo "#1. build Kraken2 db"
cd $krakendb
echo "kraken2-build --download-taxonomy --threads 16 --db $krakendb"
#kraken2-build --download-taxonomy --threads 16 --db $krakendb
echo "kraken2-build --download-library nt --db $krakendb"
kraken2-build --download-library nt --db $krakendb
echo "kraken2-build --build --threads 16 --db $krakendb &"
kraken2-build --build --threads 16 --db $krakendb &
echo "kraken2-build is running in the background and may take several days."
PID_kraken2=$!

#**********************************************
echo "#2. build Kaiju db"
cd $kaijudb
echo "kaiju-makedb -t 16 -s nr_euk &"
echo "kaiju-makedb is running in the background and may take several days."
kaiju-makedb -t 16 -s nr_euk &
PID_kaiju=$!

#**********************************************
echo "#3. update NCBI_nt"
cd $ncbi_nt
echo "update_blastdb.pl --decompress --force nt"
update_blastdb.pl --decompress --force nt

#*****************************************
echo "#4. update NCBI nr diamond format"
cd $ncbi
wget -q --show-progress ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
echo "diamond makedb --in nr.gz -d nr"
diamond makedb --in nr.gz -d nr

#*****************************************
echo "#5. update Krona taxonomy db"
echo "updateTaxonomy.sh && updateAccessions.sh"
updateTaxonomy.sh
updateAccessions.sh
#the database files all.accession2taxid.sorted and taxonomy.tab are in the folder [Krona software folder]/taxonomy

#*****************************************
echo "#6. download NCBI taxonomy db"
cd $taxondb

echo "cp $kaijudb/nr_euk/prot.accession2taxid $taxondb/"
if [ -f $kaijudb/nr_euk/prot.accession2taxid ]; then
    cp $kaijudb/nr_euk/prot.accession2taxid $taxondb/
else
    wget -q --show-progress https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
    gunzip prot.accession2taxid.gz
fi

echo "cp $krakendb/taxonomy/*.dmp $taxondb/"
if [[ -f $krakendb/taxonomy/nodes.dmp && -f $krakendb/taxonomy/names.dmp ]]; then
    cp $krakendb/taxonomy/nodes.dmp $taxondb/
    cp $krakendb/taxonomy/names.dmp $taxondb/
else
    wget -q --show-progress https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
    tar -xvzf taxdump.tar.gz
fi

echo "cp $krakendb/taxonomy/nucl_gb.accession2taxid $taxondb/"
if [ -f $krakendb/taxonomy/nucl_gb.accession2taxid ]; then
    cp $krakendb/taxonomy/nucl_gb.accession2taxid $taxondb/
else
    wget -q --show-progress https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
    gunzip nucl_gb.accession2taxid.gz
fi

#*****************************************
echo "#7. build viral reference blastn db"
cd $ncbi
wget -q --show-progress -r -nd -np -A '*.genomic.fna.gz' ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral
zcat viral.*.genomic.fna.gz > refseq_viral_genomic.fa
echo "makeblastdb -in refseq_viral_genomic.fa -dbtype nucl"
makeblastdb -in refseq_viral_genomic.fa -dbtype nucl

#<<comment
#get viral reference taxonomy using kaiju-addTaxonNames
cd $taxondb
grep ">" $ncbi/refseq_viral_genomic.fa | cut -d " " -f 1 | sed 's/>//g' > refseq_viral_genomic.id.txt
parallel -k --pipepart -a nucl_gb.accession2taxid --block 100M fgrep -F -f refseq_viral_genomic.id.txt | cut -f 2,3 > refseq_viral_genomic.taxonId.txt
awk -v FS='\t' -v OFS='\t'  '{print "C",$1,$2}'  refseq_viral_genomic.taxonId.txt > refseq_viral_genomic.kaiju.out 
kaiju-addTaxonNames -t nodes.dmp -n names.dmp -i refseq_viral_genomic.kaiju.out -o refseq_viral_genomic.gb_taxon1.txt -p
cut -f 2- refseq_viral_genomic.gb_taxon1.txt > refseq_viral.gb_taxon.txt
rm -f refseq_viral_genomic.id.txt
rm -f refseq_viral_genomic.taxonId.txt
rm -f refseq_viral_genomic.gb_taxon1.txt
rm -f refseq_viral_genomic.kaiju.out
#comment

#*****************************************
echo "#8. build plant virus blastn db"
cd $ncbi
#get virus taxon id
cut -f1 $pv_taxonFile | sed '/TaxonId/d'> $taxondb/taxIDs
#extract virus sequences from NCBI nt database
blastdbcmd -db $ncbi_nt/nt -dbtype nucl -taxidlist $taxondb/taxIDs -out plantvirus.0.fa

#remove unwanted sequences
filterbyname.sh ow=t include=f in=plantvirus.0.fa out=plantvirus.fa names=$excludeSeqId
#filterbyname.sh is in bbmp

#make blastn database
echo "makeblastdb -in plantvirus.fa -dbtype nucl"
makeblastdb -in plantvirus.fa -dbtype nucl

#<<comment
#get sequence taxonomy information
#get viral reference taxonomy using kaiju-addTaxonNames
cd $taxondb
grep ">" $ncbi/plantvirus.fa | cut -d " " -f 1 | sed 's/>//g' > seqid
parallel -k --pipepart -a nucl_gb.accession2taxid --block 100M fgrep -F -f seqid | cut -f 2,3 > taxonId.txt
awk -v FS='\t' -v OFS='\t'  '{print "C",$1,$2}'  taxonId.txt > kaiju.out 
kaiju-addTaxonNames -t nodes.dmp -n names.dmp -i kaiju.out -o gb_taxon1.txt -p
cut -f 2- gb_taxon1.txt > plantvirus.gb_taxon.txt

#remove temp files
rm -rf $ncbi/plantvirus.0.fa
rm -rf taxIDs
rm -f seqid
rm -f taxonId.txt
rm -f gb_taxon1.txt
rm -f kaiju.out
#comment

#*****************************************
echo "#9. Download Reference Viral Database (RVDB-prot) and build diamond format db"
cd $ncbi
wget -q --show-progress https://rvdb-prot.pasteur.fr/files/U-RVDB"$rvdb_version"-prot.fasta.xz  #v25.0
xz --decompress U-RVDB"$rvdb_version"-prot.fasta.xz
if [ -s RVDB.fasta ]; then
    rm -rf RVDB.fasta
fi
ln -s U-RVDB"$rvdb_version"-prot.fasta RVDB.fasta
echo "diamond makedb --in RVDB.fasta -d rvdb"
diamond makedb --in RVDB.fasta -d rvdb


cd $taxondb
#get RVDB accession taxonomy using kaiju-addTaxonNames
grep ">" $ncbi/RVDB.fasta | cut -d "|" -f 3 > RVDB-prot.id.txt
parallel -k --pipepart -a prot.accession2taxid --block 100M fgrep -F -f RVDB-prot.id.txt | cut -f 2,3 > RVDB-prot.taxonId.txt
awk -v FS='\t' -v OFS='\t'  '{print "C",$1,$2}'  RVDB-prot.taxonId.txt > RVDB-prot.kaiju.out 
kaiju-addTaxonNames -t nodes.dmp -n names.dmp -i RVDB-prot.kaiju.out -o RVDB-prot.gb_taxon1.txt -p
cut -f 2- RVDB-prot.gb_taxon1.txt > rvdb.gb_taxon.txt
rm -f RVDB-prot.id.txt
rm -f RVDB-prot.taxonId.txt
rm -f RVDB-prot.gb_taxon1.txt
rm -f RVDB-prot.kaiju.out

#**********************************************
echo "#10. download Eukaryote ribosomal RNA database"
cd $rrnadb
if [ ! -s silva-euk_combined_rRNA.fasta ]; then
    wget -q --show-progress https://github.com/biocore/sortmerna/raw/master/data/rRNA_databases/silva-euk-18s-id95.fasta
    wget -q --show-progress https://github.com/biocore/sortmerna/raw/master/data/rRNA_databases/silva-euk-28s-id98.fasta
    cat silva-euk-* > silva-euk_combined_rRNA.fasta
fi

#**********************************************
echo "#11. get micobial taxon ids"
cd $ncbi
get_species_taxids.sh -t 2 > bacteria.tids
get_species_taxids.sh -t 10239 > viruses.tids
get_species_taxids.sh -t 4751 > fungi.tids
get_species_taxids.sh -t 4762 > oomycetes.tids
cat fungi.tids bacteria.tids viruses.tids oomycetes.tids > microbial.tids

wait $PID_kraken2
echo "Kraken2 database building has finished."
wait $PID_kaiju
echo "Kaiju database building has finished."

#Clean up
rm -rf $ncbi/*.gz
rm -rf $taxondb/*.gz
if [ -s $krakendb/hash.k2d ]; then
    echo "Kraken2 database building succeeded."
    rm -rf $krakendb/library
    rm -rf $krakendb/taxonomy
else
    echo "Kraken2 database building failed, please build it following the Kraken2 manual."
fi
if [ -s $kaijudb/nr_euk/kaiju_db_nr_euk.fmi ]; then
    echo "Kaiju database building succeeded."
    mv $kaijudb/nr_euk/kaiju_db_nr_euk.fmi $kaijudb/
    rm -rf $kaijudb/nr_euk
    rm -rf $kaijudb/taxdump.tar.gz
else
    echo "Kaiju database building failed, please build it following the Kaiju manual."
fi

echo "Database building has finished."
echo "#*****************************************************"
echo "Please update database paths in the config.yaml"
echo "krakenDb: $krakendb"
echo "kaijuDb: $kaijudb/kaiju_db_nr_euk.fmi"
echo "blastnDb: $ncbi_nt/nt"
echo "blastxDb: $ncbi/nr.dmnd"
echo "#blastnViralDb: $ncbi/refseq_viral_genomic.fa"
echo "blastnViralDb: $ncbi/plantvirus.fa"
echo "#blastnViralTaxonDb: $taxondb/refseq_viral.gb_taxon.txt"
echo "blastnViralTaxonDb: $taxondb/plantvirus.gb_taxon.txt"
echo "blastxViralDb: $ncbi/rvdb.dmnd"
echo "blastxViralTaxonDb: $taxondb/rvdb.gb_taxon.txt"
echo "euk_rRNA: $rrnadb/silva-euk_combined_rRNA.fasta"
echo "taxDb: $taxondb"
echo "microbialTaxon: $ncbi/microbial.tids"
echo "acronymDb: $phytopipe_dir/db/ICTV_virus_acronym2019.txt"
echo "monitorPathogen: $phytopipe_dir/db//monitorPathogen.txt"
echo "control: $phytopipe_dir/db/phi-X174.fasta"
