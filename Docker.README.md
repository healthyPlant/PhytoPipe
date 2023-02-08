# PhytoPipe Docker

This repository is designed to show examples of ways in which you can utilize
the Docker images that come prepackaged with
[PhytoPipe](https://github.com/healthyPlant/PhytoPipe).


## Bioinformatics tools in the image

|Tool|Version| 
--- | --- | 
Python3|v3.8.10|
Java|openjdk 11.0.17 2022-10-18|
BBtools|v39.01|
FastQC|v0.11.9|
Trimmomatic|v0.39|
QUAST|v5.2.0|
snakemake|v7.18.0|
SPAdes|v3.13.1|
multiqc|v1.14|
Kraken2|v2.0.8-beta|
Kaiju|v1.9.2|
MUMmer4|v4.0.0rc1|
Blast|v2.13.0+|
Diamond|v0.9.30|
BWA|v0.7.17-r1188|
Bowtie2|v2.3.5.1|
bcftools|v1.16|
samtools|v1.16.1|
bedtools|v2.27.1|
Seqtk|v1.3-r106|
Trinity|v2.15.0|
KronaTools|v2.8.1|
 | | 

## Examples to use the image (Linux)
### 1. use the interactive bash
```shell
workDir=/my/work/directory
docker run -it --rm -v $workDir:/data xhu556/phytopipe /bin/bash
#/data is a directory in the image and links to your local workDir (/my/work/directory) through the parameter "-v". You can run any Linux commands in the docker bash, for example list all files in /my/work/directory.
ls -al /data
```
Run all PhytoPipe commands in [PhytoPipe wiki](https://github.com/healthyPlant/PhytoPipe/wiki).
```shell
#for example, build nr Diamond database
mkdir -p /data/phytopipe_db/ncbi
cd /data/phytopipe_db/ncbi
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
diamond makedb --in nr.gz -d nr
ls -ltr
``` 
You can see the file nr.dmnd in the local folder /my/work/directory.

PhytoPipe and other tools are installed in the /opt folder. You can run PhytoPipe after you build all required [databases](https://github.com/healthyPlant/PhytoPipe/wiki#databases).
Suppose your databases are in the folder /my/work/directory/phytopipe_db and your sequence files are in the folder /my/work/directory/phytopipe_run/raw. Remember these two folders in the docker container are under /data folder (we use -v to build a link between the local folder /my/work/directory and the docker folder /data). The folder phytopipe_run is the phytopipe work directory. Now you can modify the config file /opt/phytopipe/config.yaml.
```shell
nano /opt/phytopipe/config.yaml
#Select sequence type (pe or se); change database paths to /data/phytopipe_db and make sure that files are available. 

#dry run phytopipe
snakemake  --configfile /opt/phytopipe/config.yaml -s /opt/phytopipe/Snakefile --config workDir=/data/phytopipe_run --cores 32 -n

#if dry run success, run 
snakemake  --configfile /opt/phytopipe/config.yaml -s /opt/phytopipe/Snakefile --config workDir=/data/phytopipe_run --cores 32
```
The final results are in the local folder /my/work/directory/phytopipe_run.  Please see the more information in [PhytoPipe wirki](https://github.com/healthyPlant/PhytoPipe/wiki#output). 

### 2. Use the docker container as a second system
You can run the interactive bash in a docker container. Every time you start this docker container, your previous works (for example, installed tools) are still there.
Here are some steps to use our phytopipe docker image.
#### i). Create a phytopipe docker container
```
workDir=/my/work/directory  #local folder
docker run -it -v $workDir:/data xhu556/phytopipe /bin/bash 
#do something
touch testmycontainer
exit 
```
#### ii). find your docker container id
```
docker ps -a
```
Suppose your "CONTAINER ID" is a2b0f6be4c72
#### iii). start your docker container 
```
docker start a2b0f6be4c72
docker exec -it a2b0f6be4c72 /bin/bash
```
#### iv). work on the docker container 
```
ls -al testmycontainer
#You should see your previous work (testmycontainer)
cd /data
touch test
#please work on /data folder in the container. So you can access all your works on your local folder /my/work/directory. 
#You should see the file "test" in your local folder /my/work/directory.
```
Start here, you can build databases and run PhytoPipe in this container following [PhytoPipe wiki](https://github.com/healthyPlant/PhytoPipe/wiki).

#### v). run PhytoPipe on the docker container 
Before running PhytoPipe, please put your fastq files in /my/work/directory/phytopipe_run/raw (don not change 'raw', which is the default fastq file folder name PhytoPipe use ). After you modify /opt/phytopipe/config.yaml, you can run PhytoPipe using the following commands
```shell
#dry run
snakemake  --configfile /opt/phytopipe/config.yaml -s /opt/phytopipe/Snakefile --config workDir=/data/phytopipe_run --cores 32 -n
#if dry run success, please run
snakemake  --configfile /opt/phytopipe/config.yaml -s /opt/phytopipe/Snakefile --config workDir=/data/phytopipe_run --cores 32 
```
After successfully running, you should see the report in the folder /my/work/directory/phytopipe_run. **/data is a folder name in the container for the local folder  /my/work/directory.**

### 3. Directly run a command in the image
This is example to run KronaTools to generate a Krona pie chart for taxonomy.
```shell
docker run --rm -v /my/local/software/KronaTools-2.8.1/taxonomy:/opt/KronaTools-2.8.1/taxonomy -v /my/workDir/Dataset1:/data xhu556/phytopipe ktImportTaxonomy -m 3 -t 5 -o /data/Dataset1.kraken2.krona.html /data/Dataset1.kraken2.report.txt
```
"-v /my/local/software/KronaTools-2.8.1/taxonomy:/opt/KronaTools-2.8.1/taxonomy" builds a bridge between the local folder /my/local/software/KronaTools-2.8.1/taxonomy and the folder /opt/KronaTools-2.8.1/taxonomy in the image.  Any files in the local folder can be accessed by the docker container. Doing this allows that the Krona taxonomy database in the local  machine can be used by the docker. 

"-v /my/workDir/Dataset1:/data" links your local work directory (Dataset1) to the folder /data in the docker container.

"xhu556/phytopipe" is the docker image name.

"ktImportTaxonomy" is the command to run.

" -m 3 -t 5 " is the parameters for the command.

"-o /data/Dataset1.kraken2.krona.html" is the output. It is in the folder /data in the docker. But you can see it in your local work directory because of the bridge we built.

"/data/Dataset1.kraken2.report.txt" is the input from Kraken2 classification.  The file "Dataset1.kraken2.report.txt" is in the local work directory. 


### 4. Run PhytoPipe 
#### i). Create your config.yaml file
You can download "config.docker.yaml" from [PhytoPipe](https://github.com/healthyPlant/PhytoPipe) to edit or copy the following code to edit
```shell
#For single-read file
#seq_type: 'se'  #'se' for single-readd, 'pe' for paired-end  
#strand1: ''  #'' for single-read file name without R1, 'R1' for single-read file name end with R1 or R1_001
#strand2: 'R2'  # keep it even no pair R2 file

#For paired-end file
seq_type: 'pe'  #'se' for single-read, 'pe' for paired-end  
strand1: 'R1'  # for pair-end file name with R1 and R2 or R1_001 and R2_001
strand2: 'R2'   #for R2 or R2_001

input_format: 'fastq.gz'  #read file format. 'fq.gz' is an option 
number_of_threads: 16 #number of computer cores or threads

#please check wiki to see how to build the databases
#please do not change database paths here. They are the paths in the docker image. You can use "-v" to link the local databases 
#spike-in or contorl's sequence file 
control: /opt/phytopipe/db/phi-X174.fasta
#kraken2 database folder
krakenDb: /opt/phytopipe/db/kraken_db
#kaiju database file
kaijuDb: /opt/phytopipe/db/kaiju_db/kaiju_db_nr_euk.fmi
#NCBI blastn nt database downloaded from NCBI
blastnDb: /opt/phytopipe/db/ncbi_nt/nt
#NCBI nr database made by Diamond
blastxDb: /opt/phytopipe/db/ncbi/nr.dmnd
#viral reference database
blastnViralDb: /opt/phytopipe/db/ncbi/refseq_viral_genomic.fa
#blastnViralDb:/ppq/data2/resources/ncbi/plantvirus.fa
#viral reference taxonomy onformation
blastnViralTaxonDb: /opt/phytopipe/db/ncbi/taxonomy/refseq_viral.gb_taxon.txt
#blastnViralTaxonDb: /ppq/data2/resources/ncbi/taxonomy/plantvirus.gb_taxon.txt
#RVDB viral protein database made by Diamond
blastxViralDb: /opt/phytopipe/db/ncbi/rvdb.dmnd
#RVDB viral protein sequences taxonomy information
blastxViralTaxonDb: /opt/phytopipe/db/ncbi/taxonomy/rvdb.gb_taxon.txt
#SILVA 18S and 28S Eukaryote ribosomal RNA database
euk_rRNA: /opt/phytopipe/db/rRNA/silva-euk_combined_rRNA.fasta
#NCBI taxonomy database
taxDb: /opt/phytopipe/db/ncbi/taxonomy

#Virus acronym from ICTV
acronymDb: /opt/phytopipe/db/ICTV_virus_acronym.txt

#microbial taxon ids
microbialTaxon: /opt/phytopipe/db/ncbi/microbial.tids

#your monitored pathogens, if you don't have them, leave the file empty
#the file is a tab delimited
#column titles are:TaxonId Species
#each line have a specices taxonomy ID and name
monitorPathogen: /opt/phytopipe/db/monitorPathogen.txt

#Blast databases
#'virus' is for blasting against RVDB and viral references.
#'all' is for blasting NCBI nt and nr. It could take a long time to finish blastn
blastDbType: virus # or 'all' 
#Blast E-value cutoff 
blastEvalue: 1e-10

#reads used for assembly or mapping to a reference
#'trimmed' is for all trimmed reads including host reads
#'clean' is for all possible pathogen reads. host reads are supposed to be removed
mapReadType: 'clean' #'trimmed'

#trimmomatic, its adapters and parameters
trimmomatic: /opt/Trimmomatic-0.39/trimmomatic-0.39.jar
adapters: /opt/Trimmomatic-0.39/adapters/TruSeq3-SE.fa
trimmomatic_param: " LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 "

#de novo genome assembly tool
#'Trinity' or 'Spades'
assembler: Trinity # 'Spades'
#Spades parameters. Different kmers (k) may generate different results. 
spades_param: " --only-assembler --phred-offset 33 -k 21,51,71,91 "


#mapping tool
#'bwa' or 'bowtie'
mappingTool: bwa # 'bowtie2'
#mapping tool parameters
bwa_param: " -B 4 "  #default bwa setting #" -k 12 -A 1 -B 3 -O 1 -E 1 " #bwa loose mapping; 
bowtie2_param: " --very-sensitive-local -k 100 --score-min L,20,1.0 " #bowtie2 sensitive mapping; " --mp 20  --score-min L,-0.1,-0.1 " #bowtie2 very stringent mapping

#clumpify parameters
clumpify_param: "dedupe=t subs=0 passes=2 " #dupedist=40 optical=t # These parameters identify reads as duplicated only if they are an exact match (i.e., no substitution allowed).

#gobal variables, please do not change
#they can be changed in the command line
fastqDir: ""
flowCellDir: ""
workDir: ""
samples: ""

run_info:
  raw: raw                                #raw fastq file directoy
  log: logs                               #log files diretory
  qc: qc                                  #fastqc and multiqc directory
  trim: trimmed                           #trimmed directory
  clean: cleaned                          #clean reads after host, control ex. Phx174, duplicates removal
  assemble: assembly                      #assemble directory
  classify: classification                #read classification directory
  annotate: annotation                    #blast annotation directory  
  map: mapping                            #mapping reads directory   
  report: report                          #report directory    
  novel: novelVirus                       #novelVirus directory      
```
**Important:** If you use single-reads, please use the single-read setting and comment paired-end setting
```shell
seq_type: 'se' 
strand1: ''  
strand2: 'R2'
```
If you use paired-end reads, please use the paired-end setting and comment single-read setting
```shell
seq_type: 'pe'  
strand1: 'R1'  
strand2: 'R2' 
```

#### ii). Get your local databases
```shell
#kraken2 database folder
krakenDb_dir=/my/db/path/kraken_db #hash.k2d  opts.k2d  seqid2taxid.map  taxo.k2d 
#kaiju database file
kaijuDb_dir=/my/db/path/kaiju_db  #/kaiju_db_nr_euk.fmi
#NCBI blastn nt database downloaded from NCBI
ncbi_nt_dir=/my/db/path/ncbi_nt  #/nt
#SILVA 18S and 28S Eukaryote ribosomal RNA database
euk_rRNA_dir=/my/db/path/ #/silva-euk_combined_rRNA.fasta
#NCBI db
ncbiDb_dir=/my/db/path/ncbi #nr.dmnd refseq_viral_genomic.fa  rvdb.dmnd
#NCBI taxonomy database
taxDb_dir=/my/db/path/taxonomy  #refseq_viral.gb_taxon.txt rvdb.gb_taxon.txt nodes.dmp names.dmp
#Krona taxonomy
kronaTaxDb_dir=/my/software/KronaTools-2.8.1/taxonomy  #all.accession2taxid.sorted taxonomy.tab
```

#### iii). Run PhytoPipe
Your sequence fastq.gz files are in the folder $workDir/raw. Here we use VIROMOCKchallenge [Dataset8](https://gitlab.com/ilvo/VIROMOCKchallenge)(paried-end reads) and [Dataset10](https://gitlab.com/ilvo/VIROMOCKchallenge) (single-reads, novel virus) to test the pipeline.
``` shell
workDir=/home/myname/Dataset8   
config=/home/myname/Dataset8/config.docker.yaml

#for Dataset10
#workDir=/home/myname/Dataset10   
#config=/home/myname/Dataset10/config.docker.yaml

#fastq file must be in $workDir/raw folder
ls -al $workDir/raw
#run docker
#use "-d" for runing in the background
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
```
After successfully running, you should see the report in the $workDir.

You can download "runDocker.sh" from Github [PhytoPipe](https://github.com/healthyPlant/PhytoPipe) to edit and use it.
```
bash ./runDocker.sh
```
