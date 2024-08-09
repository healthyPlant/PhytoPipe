# PhytoPipe: A Phytosanitary Pipeline for plant pathogen detection and diagnosis using RNA-seq data

PhytoPipe is an open-source bioinformatics pipeline for plant pathogen detection using RNA-seq data. The pipeline is written in [Snakemake](https://snakemake.readthedocs.io) (Köster and Rahmann 2018), a workflow management system for the development of data analysis workflows. PhytoPipe is divided into read quality control, read taxonomic classification, *de novo* assembly, contig blast annotation, virus species selection, viral reference mapping, and summary report. PhytoPipe can detect plant pathogens: bacteria (in addition to low expressed phytoplasma), fungi, viruses, and viroids. PhytoPipe has been tested on RNA-seq data (Illumina single-end and paired-end reads of different lengths) from different quarantined plant species such as apple, pear, peach, potato, sweet potato, cassava, rice, sugarcane, and bamboo. PhytoPipe is being used by the USDA Plant Germplasm Quarantine Program (PGQP) as the cornerstone of HTS-based diagnostics.  

# Workflow

![scheme of workflow](doc/VDflow_scheme.png?raw=true)

## Tools in pipeline

1. Convert Blc to FastQ files ([`bcl2fastq`](https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html), *optional*)
2. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
3. Read preprocessing
	1. Host ribosomal RNA removal (*bbduk.sh* in [`BBTools`](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/))
	2. PCR duplicates removal (*clumpify.sh* in [`BBTools`](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/))
	3. Control/spike-in removal (*bbsplit.sh* in [`BBTools`](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/))
	4. Adaptor and read trimming ([`Trimmomatic`](http://www.usadellab.org/cms/?page=trimmomatic))
4. Read classification
	1. Nucleotide level ([`Kraken2`](http://ccb.jhu.edu/software/kraken2/))
	2. Protein level ([`Kaiju`](https://kaiju.binf.ku.dk/))
	3. View result with Krona pie chart ([`KronaTools`](https://github.com/marbl/Krona/tree/master/KronaTools))
5. Assembly-based annotation
   1. Removal of host reads (modified [`KrakenTools`](https://ccb.jhu.edu/software/krakentools/) (included) and [`seqtk`](https://github.com/lh3/seqtk))
   2. *De novo* assembly tools ([`Trinity`](https://github.com/trinityrnaseq/trinityrnaseq/wiki) *||* [`SPAdes`](http://cab.spbu.ru/software/spades/))
   3. Assembly assessment report ([`QUAST`](http://quast.sourceforge.net/quast))
   4. Blastn to [`viral reference genome`](https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/) or plant viral database or [`NCBI nt database`](https://ftp.ncbi.nlm.nih.gov/blast/db/) ([`blastn`](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch))
   5. Diamond blastx to [`RVDB`](https://rvdb-prot.pasteur.fr/) or [`NCBI nr database`](https://ftp.ncbi.nlm.nih.gov/blast/db/) ([`Diamond`](https://github.com/bbuchfink/diamond))
   6. View conig classification with Krona pie chart ([`KronaTools`](https://github.com/marbl/Krona/tree/master/KronaTools))
6. Reference-based mapping
   1. Read alignment ([`BWA`](https://github.com/lh3/bwa) *||* [`Bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml))
   2. Sort and index alignments ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
   3. Consensus callers ([`BCFTools`](http://samtools.github.io/bcftools/bcftools.html) and [`BEDTools`](https://github.com/arq5x/bedtools2/))
   4. Make pseudo contig (modified *order_and_orient* with [`MUMmer4`](https://mummer4.github.io/) in [`viral-NGS`](https://github.com/broadinstitute/viral-ngs) (included))
7. Present multisample QC and visualisation for raw read and assembly ([`MultiQC`](http://multiqc.info/))


## Quick start
### Installation

For a Windows system, please use [PhytoPipe docker image](https://hub.docker.com/r/healthyplant/phytopipe).

For a Linux or Mac system, you can download our code or clone the repository using Git:
```
cd /path/to/software
git clone https://github.com/healthyPlant/PhytoPipe.git
```

Then install dependencies. For an Ubuntu system, simply run
```
sudo bash /path/to/PhytoPipe/scripts/installTools.sh /path/to/software
```
For other Linux or Mac systems, please follow [PhytoPipe wiki](https://github.com/healthyPlant/PhytoPipe/wiki#dependencies) to install dependencies.

***Since some tools conflict in a conda environment and are hard to update, we do not recommend installing all dependencies using conda.*** For a conda user, please check [PhytoPipe wiki](https://github.com/healthyPlant/PhytoPipe/wiki#dependencies) to see how to create a PhytoPipe conda environment.

### Build databases
To install/update all databases, simply run
```
nohup bash /path/to/PhytoPipe/scripts/updateDatabase.sh /path/to/software/PhytoPipe /path/to/database v25.0 &
```
Please find the newest RVDB version from [RVDB database, protein version](https://rvdb-prot.pasteur.fr/) and use the full path in the command. v25.0 in the above command is the RVDB version.  

### Set up configuration
Customize the workflow based on your needs in `config.yaml`, such as your sequence file extension and database paths. Please see the details in [wiki](https://github.com/healthyPlant/PhytoPipe/wiki).

### Run PhytoPipe
Please check [dependencies and databases requirements](https://github.com/healthyPlant/PhytoPipe/wiki) first using a dry-run (-n flag). If they are met, then you can run PhytoPipe

#### 1. Run whole PhytoPipe
For fastq.gz reads input dry-run:
```shell
$ snakemake  --configfile /path/to/PhytoPipe/config.yaml -s /path/to/PhytoPipe/Snakefile --config workDir=/path/to/output/folder fastqDir=/path/to/input/fastq/folder --cores [number of cores ex. 16] -n 
```

For raw Illumina flowcell's data input dry-run:
```shell
$ snakemake  --configfile /path/to/PhytoPipe/config.yaml -s /path/to/PhytoPipe/Snakefile --config workDir=/path/to/output/folder flowCellDir=/path/to/Illumina/run --cores [number of cores ex. 16] -n
```

If the dry-run succeeds, please remove '-n' parameter to run the pipeline. If you'd like to run it in the background, please use 'nohup'. For example:
```shell
$ nohup snakemake  --configfile /path/to/PhytoPipe/config.yaml -s /path/to/PhytoPipe/Snakefile --config workDir=/path/to/output/folder fastqDir=/path/to/input/fastq/folder --cores [number of cores ex. 16] &
```
**Important:** For workDir and fastqDir paths, full paths must be used.

You can view progress or errors in the file 'nohup.out' using the command

`more nohup.out`

#### 2. Run partial PhytoPipe: NCBI nt and nr blast after the above viral pipeline 
After finishing the virus detection using the virus database and RVDB, if you want to check NCBI nt and nr databases for other pathogens, you can try   
```shell
$ snakemake  --configfile your/config.yaml -s /path/to/PhytoPipe/Snakefile_runBlast --config workDir=/path/to/output/folder blastDbType=all --cores [number of cores ex. 16] -n 
```

#### 3. Run partial PhytoPipe: trim-assemble-annotation pipe
If your computer memory is limited (<256Gb), you can skip the classification steps (Kraken2 and Kaiju) to directly run trim-assemble-annotation, then check your results in the annotation folder (sample.blast.nx.txt, sample.blastnt0.txt, sample.blastnr0.txt). Here is a dry-run command example for running Trimmomatic trim, Spades assemble, and Blast against viral databases and NCBI nt/nr database annotation:
```shell
$ snakemake  --configfile /path/to/PhytoPipe/config.yaml -s /path/to/PhytoPipe/Snakefile_runAnnotation --config workDir=/path/to/output/folder fastqDir=/path/to/input/fastq/folder mapReadType=trimmed assembler=Spades blastDbType=all --cores [number of cores ex. 16] -n 
```
Note: seq_type (se/pe), mapReadType, assembler, blastDbType ban be changed in config.yaml

#### 4. Run partial PhytoPipe: trim-assemble-annotation pipe with only NCBI nt and nr databases
If you just want annotate the contigs using NCBI nt (directly download) and nr databases (Diamond built), you can run another trim-assemble-annotation pipeline, then check your results in the annotation folder (sample.blastnt0.txt, sample.blastnr0.txt). Here is a dry-run command example for running Trimmomatic trim, Spades assemble, and Blast against NCBI nt/nr database annotation:
```shell
$ snakemake  --configfile /path/to/PhytoPipe/config.yaml -s /path/to/PhytoPipe/Snakefile_runAnnotationNt --config workDir=/path/to/output/folder fastqDir=/path/to/input/fastq/folder mapReadType=trimmed assembler=Spades --cores [number of cores ex. 16] -n 
```

### PhytoPipe quick test
After the software and databases are ready, you can run a quick test using [the VIROMOCKchallenge Dataset8](https://gitlab.com/ilvo/VIROMOCKchallenge/-/blob/master/Datasets/Dataset8.md), which is in the test/data folder. The steps can be found in [test notebook](test/quick_test.ipynb).  

# PhytoPipe docker image
[The PhytoPipe docker image](https://hub.docker.com/r/healthyplant/phytopipe) can be pulled 
```
docker pull healthyplant/phytopipe
```
Docker can help you avoid manually installing the software. You can use the docker image on many systems (Linux, Mac, Windows). PhytoPipe docker image usage is in its [docker README](https://hub.docker.com/r/healthyplant/phytopipe). PhytoPipe docker tutorial can be found in [test docker](test/PhytoPipe_docker_tutorial.ipynb).

# Documentation

More information on input/output, dependencies, and databases are found in the [wiki](https://github.com/healthyPlant/PhytoPipe/wiki)



# Citation

Hu, X., Hurtado-Gonzales, O.P., Adhikari, B.N. et al. PhytoPipe: a phytosanitary pipeline for plant pathogen detection and diagnosis using RNA-seq data. BMC Bioinformatics 24, 470 (2023). https://doi.org/10.1186/s12859-023-05589-2

Contact
------------
Alex Hu (xiaojun.hu (at) usda.gov)
