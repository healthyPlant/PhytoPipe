#update taxonomy
#************************************************************
taxondb=/ppq/data2/resources/ncbi/taxonomy
krakendb=/ppq/data2/resources/kraken_db
kaijudb=/ppq/data2/resources/kaiju_db
blastxTaxondb=$taxondb/rvdb.gb_taxon.txt
blastnTaxondb=$taxondb/refseq_viral.gb_taxon.txt
blastdb=/ppq/data2/resources/ncbi
#rvdb=/ppq/data2/resources/ncbi/RVDB.fa
rvdbFile=U-RVDBv18.0-prot.fasta  #get name from RVDB website https://rvdb-prot.pasteur.fr/

#1. Update krona taxonomy
conda activate ppq_pipeline
ktUpdateTaxonomy.sh
#/ppq/data1/software/anaconda3/envs/ppq_pipeline/bin/ktUpdateTaxonomy.sh
#location: /ppq/data1/software/anaconda3/envs/ppq_pipeline/opt/krona/taxonomy
#file: /ppq/data1/software/anaconda3/envs/ppq_pipeline/opt/krona/taxonomy/taxonomy.tab
#Fetching taxdump.tar.gz from ftp://ftp.ncbi.nih.gov/pub/taxonomy

#2. Kraken
# https://ccb.jhu.edu/software/kraken2/index.shtml?t=manual#kraken-2-databases

#It requires hash.k2d, opts.k2d, and taxo.k2d files. They are in /ppq/data2/resources/kraken_db.

#1). Install a taxonomy.
if [ -d $krakendb/taxonomy ]; then
	rm -rf $krakendb/taxonomy
fi
kraken2-build --download-taxonomy --threads 36 --db $krakendb
#This will download the accession number to taxon maps, as well as the taxonomic name and tree information from NCBI. These files can be found in $DBNAME/taxonomy/

#2). Install one or more reference libraries
kraken2-build --download-library nt --threads 36 --db $krakendb
# nt: NCBI non-redundant nucleotide database

#build the database
if [ -f $krakendb/seqid2taxid.map ]; then
	rm -rf $krakendb/seqid2taxid.map
	rm -rf $krakendb/*.k2d
fi

kraken2-build --build --threads 36 --db $krakendb

#3. Kaiju taxonomy
#download nr_euk http://kaiju.binf.ku.dk/database/kaiju_db_nr_euk_2019-06-25.tgz in the Kaiju server http://kaiju.binf.ku.dk/server
#nr_euk: Subset of NCBI BLAST nr database containing all proteins belonging to Archaea, Bacteria and Viruses, fungi and microbial eukaryotes
#kaiju-makedb -s nr_euk
#cd kaijudb
#wget http://kaiju.binf.ku.dk/database/kaiju_db_nr_euk_2019-06-25.tgz
#tar zxvf kaiju_db_nr_euk_2019-06-25.tgz
#nr_euk folder and nodes.dmp, names.dmp are extracted

#4. build viral Blast taxononmy
# RVDB for blastx

#get taxonomy information for refseq_viral
cd $blastdb
#download refseq_viral from genBank refSeq
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.1.genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.3.1.genomic.fna.gz
zcat viral.1.1.genomic.fna.gz viral.2.1.genomic.fna.gz viral.3.1.genomic.fna.gz  > refseq_viral_genomic.fa
makeblastdb -in refseq_viral_genomic.fa -input_type fasta -dbtype nucl 

#download the latest RVDB from website https://rvdb.dbi.udel.edu/  #RVDB nt, don't use it for diamond
#wget https://rvdb.dbi.udel.edu/download/$rvdbFile.gz
#gunzip $rvdbFile.gz
#wget https://rvdb-prot.pasteur.fr/files/U-RVDBv18.0-prot.fasta.bz2 #RVDB nr
wget https://rvdb-prot.pasteur.fr/files/$rvdbFile.bz2
bzip2 -d $rvdbFile.bz2
ln -s $rvdbFile RVDB.fa
#make diamond format db
diamond makedb --in RVDB.fa -d rvdb

#get taxonomy files
cd $taxondb
if [ -f $taxondb/nodes.dmp ]; then
	rm $taxondb/nodes.dmp
	ln -s $krakendb/taxonomy/nodes.dmp $taxondb/nodes.dmp
fi
if [ -f $taxondb/names.dmp ]; then
	rm $taxondb/names.dmp
	ln -s $krakendb/taxonomy/names.dmp $taxondb/names.dmp
fi
if [ -f $taxondb/nucl_gb.accession2taxid ]; then
	rm $taxondb/nucl_gb.accession2taxid
	ln -s $krakendb/taxonomy/nucl_gb.accession2taxid $taxondb/nucl_gb.accession2taxid
fi

if [ -f $taxondb/prot.accession2taxid ]; then
	rm $taxondb/prot.accession2taxid
fi
#download prot_accession2taxid
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
gunzip prot.accession2taxid.gz #-c

grep ">" $blastdb/refseq_viral_genomic.fa | cut -d " " -f 1 | sed 's/>//g' > refseq_viral_genomic.id.txt

#get viral refSeq taxon id
parallel -k --pipepart -a nucl_gb.accession2taxid --block 100M fgrep -F -f refseq_viral_genomic.id.txt | cut -f 2,3 > refseq_viral_genomic.taxonId.txt
#use kaiju-addTaxonNames to add taxon to refseq_viral_genomic.taxonId.txt
awk -v FS='\t' -v OFS='\t'  '{print "C",$1,$2}'  refseq_viral_genomic.taxonId.txt > refseq_viral_genomic.kaiju.out 
kaiju-addTaxonNames -t nodes.dmp -n names.dmp -i refseq_viral_genomic.kaiju.out -o refseq_viral_genomic.gb_taxon1.txt -p
#-r species,genus,family,order,class,phylum,superkingdom
#-p full taxonomy path
cut -f 2- refseq_viral_genomic.gb_taxon1.txt > refseq_viral.gb_taxon.txt
#delete intermedia files
rm refseq_viral_genomic.gb_taxon1.txt
rm refseq_viral_genomic.kaiju.out


#get RVDB protein accession id
grep ">" $blastdb/RVDB.fa | cut -d "|" -f 3 > RVDB-prot.id.txt
#get taxon id
parallel -k --pipepart -a prot.accession2taxid --block 100M fgrep -F -f RVDB-prot.id.txt | cut -f 2,3 > RVDB-prot.taxonId.txt

#use kaiju-addTaxonNames to add taxon to RVDB-prot.taxonId.txt
awk -v FS='\t' -v OFS='\t'  '{print "C",$1,$2}'  RVDB-prot.taxonId.txt > RVDB-prot.kaiju.out 
kaiju-addTaxonNames -t nodes.dmp -n names.dmp -i RVDB-prot.kaiju.out -o RVDB-prot.gb_taxon1.txt -p
#-r species,genus,family,order,class,phylum,superkingdom
# -p            Print full taxon path.

#check missing taxon
#grep -w '0$' RVDB-prot.taxonId.txt | head

#change ; to \t
#cut -f 2- U-RVDBv17.0-prot.gb_taxon1.txt | tr ';' \\t  > U-RVDBv17.0-prot.gb_taxon.txt
cut -f 2- RVDB-prot.gb_taxon1.txt > rvdb.gb_taxon.txt
#delete intermedia files
rm RVDB-prot.gb_taxon1.txt
rm RVDB-prot.kaiju.out
