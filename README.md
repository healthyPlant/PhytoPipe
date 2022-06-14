# PhytoPipe: A Phytosanitary Pipeline for plant pathogen detection and diagnosis using RNA-seq data

PhytoPipe is an open-source bioinformatics pipeline for plant pathogen detection using RNA-seq data. The pipeline is written in [Snakemake](https://snakemake.readthedocs.io) (Köster and Rahmann 2018), a workflow management system for the development of data analysis workflows. PhytoPipe is divided into read quality control, read taxonomic classification, *de novo* assembly, contig blast annotation, virus species selection, viral reference mapping, and summary report. PhytoPipe can detect plant pathogens: bacteria(in addition to low expressed phytoplasma), fungi, viruses, and viroids. PhytoPipe has been tested on RNA-seq data (Illumina single-end and paired-end reads of different lengths) from different quarantined plant species such as apple, pear, peach, potato, sweet potato, cassava, rice, sugarcane, and bamboo. PhytoPipe is being used by the USDA Plant Germplasm Quarantine Program (PGQP) as the cornerstone of HTS-based diagnostics.  

# Workflow

![scheme of workflow](doc/VDflow_scheme.png?raw=true)

## Tools in pipeline

1. Convert Blc to FastQ files ([`bcl2fastq`](https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html), *optional*)
2. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
3. Read precprocessing
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
   4. Blastn to [`viral reference genome`](https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/) or [`NCBI nt database`](https://ftp.ncbi.nlm.nih.gov/blast/db/) ([`blastn`](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch))
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

Check [dependencies and databases requirements](https://github.com/healthyPlant/PhytoPipe/wiki), if they are met, then

Clone the repository:

`git clone https://github.com/healthyPlant/PhytoPipe.git`

Create the environment without tools if they are already installed:

`conda create -n phytopipe python=3.7`

Or create the environment with installing tools by conda:

`conda env create -f [PhytoPipe]/environment.yaml`

Activate the environment:

`conda activate phytopipe`

If you want to remove PhytoPipe, simply type the follwowing command,

`conda env remove --name phytopipe`

### Set up configuration
Customize the workflow based on your need in `config.yaml`, such as your sequence file extension and database paths. Please see the details in [wiki](https://github.com/healthyPlant/PhytoPipe/wiki).

### Run PhytoPipe

For fastq.gz reads input dry-run:
```shell
$ snakemake  --configfile /path/to/PhytoPipe/config.yaml -s /path/to/PhytoPipe/Snakefile --config workDir=/path/to/output/folder fastqDir=/path/to/input/fastq/folder --cores [number of cores ex. 16] -n 
```

For raw Illumina flowcell's data input dry-run:
```shell
$ snakemake  --configfile /path/to/PhytoPipe/config.yaml -s /path/to/PhytoPipe/Snakefile --config workDir=/path/to/output/folder flowCellDir=/path/to/Illumina/run --cores [number of cores ex. 16] -n
```

If dry-run (-n) succeeds, please remove '-n' parameter to run the pipeline. If you'd like to run it in the background, please use 'nohup' and '&'. For example:
```shell
$ nohup snakemake  --configfile /path/to/PhytoPipe/config.yaml -s /path/to/PhytoPipe/Snakefile --config workDir=/path/to/output/folder fastqDir=/path/to/input/fastq/folder --cores [number of cores ex. 16] &
```
You can view progress or errors in the file 'nohup.out'

`more nohup.out`

# Documentation

More information on input/output, dependencies, and databases are found in the [wiki](https://github.com/healthyPlant/PhytoPipe/wiki)

# Citation

> The manuscript is on the way.
> Will update later. 

Contact
------------
Alex Hu (xiaojun.hu (at) usda.gov)
