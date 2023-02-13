#/usr/bin/env bash
#please put all fastq files in $workDir/raw folder
#and modify config.docker.yaml
workDir=/home/alex/Dataset10  #Dataset8  #fastq file must be in $workDir/raw folder
config=/ppq/data0/phytopipe_docker/config.docker.yaml

#databases
#please change path to your database path
#get your local databases
#kraken2 database folder
krakenDb_dir=/my/work/directory/phytopipe_db/kraken_db #hash.k2d  opts.k2d  seqid2taxid.map  taxo.k2d 
#kaiju database file
kaijuDb_dir=/my/work/directory/phytopipe_db/kaiju_db  #/kaiju_db_nr_euk.fmi
#NCBI blastn nt database downloaded from NCBI
ncbi_nt_dir=/my/work/directory/phytopipe_db/ncbi_nt  #/nt
#SILVA 18S and 28S Eukaryote ribosomal RNA database
euk_rRNA_dir=/my/work/directory/phytopipe_db/rRNA #/silva-euk_combined_rRNA.fasta
#NCBI db
ncbiDb_dir=/my/work/directory/phytopipe_db/ncbi #nr.dmnd refseq_viral_genomic.fa  rvdb.dmnd
#NCBI taxonomy database
taxDb_dir=/my/work/directory/phytopipe_db/ncbi/taxonomy  #refseq_viral.gb_taxon.txt rvdb.gb_taxon.txt nodes.dmp names.dmp
#Krona taxonomy
kronaTaxDb_dir=/my/work/directory/phytopipe_db/krona_db  #all.accession2taxid.sorted taxonomy.tab

#use -d to run docker in background
# -u $UID:$GROUPS \
docker run -it --rm --name phytopipe \
           -v $workDir:/data \
           -v $krakenDb_dir:/opt/phytopipe/db/kraken_db \
           -v $kaijuDb_dir:/opt/phytopipe/db/kaiju_db \
           -v $ncbi_nt_dir:/opt/phytopipe/db/ncbi_nt \
           -v $euk_rRNA_dir:/opt/phytopipe/db/rRNA \
           -v $ncbiDb_dir:/opt/phytopipe/db/ncbi \
           -v $taxDb_dir:/opt/phytopipe/db/ncbi/taxonomy \
           -v $kronaTaxDb_dir:/opt/KronaTools-2.8.1/taxonomy \
           -v $config:/opt/phytopipe/config.yaml \
           xhu556/phytopipe \
           snakemake --configfile /opt/phytopipe/config.yaml \
           -s /opt/phytopipe/Snakefile \
           --config workDir=/data \
           --cores 32  
