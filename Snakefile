"""
    Snakemake main file 
"""
__author__ = 'Xiaojun Hu <xiaojun.hu@usda.gov>'

import os,subprocess,itertools,glob,re,sys
import datetime
from snakemake.utils import report

#SETUP PROCESS-RELATED CONFIGURATION 
try:
	CONFIGFILE = str(workflow.overwrite_configfiles[0])
except:
	CONFIGFILE = str(workflow.overwrite_configfile[0])

#config["workDir"] = os.getcwd()
#setup workdir dynamically
#rules will then be expressed relative to the workdir:
workdir: config["workDir"]  #can be changed at the snakemake command line
fastqDir = config["fastqDir"]  #can be changed at the snakemake command line
flowCellDir = config["flowCellDir"]  #can be changed at the snakemake command line

#Get relative position of Snakefile from wd
SNAKEFILE = workflow.snakefile
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)
CORES = workflow.cores #note workflow.cores work only on version 5.11.2 above
#print(CORES)
config["number_of_threads"] = CORES  #change threads number by command paramter --cores/-j 16
config["snakefile"] = SNAKEFILE
config["snakefile_dir"] = SNAKEFILE_DIR
#print(SNAKEFILE)

#print(workflow.__dict__) #print all attributes
#Establish snakefile and environment dictionaries
snakefiles_dir = os.path.abspath(os.path.join(SNAKEFILE_DIR, "rules"))		#directory for additional snakefiles
#environments_dir = os.path.abspath(os.path.join(SNAKEFILE_DIR, "envs"))	#directory for conda environment yaml files
scripts_dir = os.path.abspath(os.path.join(SNAKEFILE_DIR, "scripts")) 			#directory for extra scripts used in workflow

#running folders
seq_type = config["seq_type"] 
rawReadDir = config["run_info"]["raw"]
qcDir = config["run_info"]["qc"]
logDir = config["run_info"]["log"]
trimDir = config["run_info"]["trim"]
assembleDir = config["run_info"]["assemble"]
classifyDir = config["run_info"]["classify"]
annotateDir = config["run_info"]["annotate"]
mapDir = config["run_info"]["map"]
reportDir = config["run_info"]["report"]
novelDir = config["run_info"]["novel"]
blastDbType = config["blastDbType"]
reportFile = reportDir + "/report.txt" 
htmlReport = reportFile.replace(".txt", ".html")

if os.path.isfile(reportFile):
    os.remove(reportFile)

## Check whether input fastq files are compressed or not
if  "input_format" in config.keys():
    input_format = config["input_format"]
else:
    input_format = "fastq.gz"

#Run bcl2fastq
#bclOption = " --barcode-mismatches 0 --no-lane-splitting "  #for one index  
bclOption = " --barcode-mismatches 1 --no-lane-splitting "  #for dual index 

#setup soft link or copy fastq files to rawReadDir
if fastqDir:
    #set up soft link for fastq files if fastq folder is surpplied
    if not os.path.exists(rawReadDir):
        os.system("ln -sf " + fastqDir + " " + rawReadDir)
    # Check if given path is link
    if os.path.exists(rawReadDir) and os.path.islink(rawReadDir):
        print(fastqDir, " soft link is built." )
    else:
        # soft link broken, copy files
        if not os.path.isdir(rawReadDir): #check a folder exists
            os.mkdir(rawReadDir) #make a folder
        os.system("cp " + fastqDir + "/*." + input_format + " " + rawReadDir + "/")  
			
elif flowCellDir and len(os.listdir(flowCellDir)) != 0 :
    print("Runing bcl2fastq")
    if not os.path.isdir(rawReadDir): #check a folder exists
        os.mkdir(rawReadDir) #make a folder
    os.system("bcl2fastq -R " + flowCellDir + " -o " + rawReadDir + bclOption + ">> " + rawReadDir + "/bcl2fastq.log 2>&1")
    #--sample-sheet SampleSheet.csv

if not os.listdir(rawReadDir): 
    print(rawReadDir + " is empty. Exit!") 
    sys.exit()

#get sample names from rawReadDir
#print(config['workdir'])
#print(os.path.basename(config['run_info']['raw_fastq']))
#SAMPLES, = glob_wildcards(rawReadDir + "/{sample}_R1_001.fastq.gz")
SAMPLES = [os.path.basename(x) for x in glob.glob(rawReadDir + '/*.' +  input_format)]
#remove 'Undetermined_*'
SAMPLES = [x for x in SAMPLES if not x.startswith('Undetermined_')]
#print(SAMPLES)
#remove '_R1_001.fastq.gz'
SAMPLES = [x.replace("." + input_format,"") for x in SAMPLES]
if config["strand1"] != '':
	SAMPLES = [re.sub(r"(.*)_%s.*" % config["strand1"], "\\1", x) for x in SAMPLES]

if seq_type == 'pe':
	SAMPLES = [re.sub(r"(.*)_%s.*" % config["strand2"], "\\1", x) for x in SAMPLES]
SAMPLES = list(set(SAMPLES))  #remove duplicated from a list using set

print("Samples: %s"  % ", ".join(SAMPLES))
config["samples"] = SAMPLES  #save sample names in a global viriable
#print("config samples: ", config['samples'])
STRANDS = [config["strand1"], config["strand2"]]
strand2 = config["strand2"]
if len(SAMPLES) == 0:
	sys.exit("No sample files found. Exit!")

# load rules 
include: os.path.join(snakefiles_dir, "cleanReads.smk")
include: os.path.join(snakefiles_dir, "classifyReads.smk")
include: os.path.join(snakefiles_dir, "assemble.smk")
include: os.path.join(snakefiles_dir, "annotate.smk")
include: os.path.join(snakefiles_dir, "mappingReads.smk")
include: os.path.join(snakefiles_dir, "makeReport.smk")

# Single-end
if (seq_type == "se"):
    #QC raw reads
	rawFastQC = expand(qcDir + "/raw_fastqc/{sample}_fastqc.zip", sample=SAMPLES)
	#check host ribosomal RNA
	rRNACheck = expand(trimDir + "/{sample}.filtRNA.fastq.gz", sample=SAMPLES)
	#remove duplicate reads
	rmDup = expand(trimDir + "/{sample}.rmdup.fastq.gz", sample=SAMPLES) #run clumpify to remove duplicate reads
	#remove PhiX174 contaminant
	rmCtm = expand(trimDir + "/{sample}.rmdup_ctm.fastq.gz", sample=SAMPLES) #run clumpify to remove duplicate reads
	#Trim reads
	trim = expand(trimDir + "/{sample}.trimmed.fastq.gz", sample=SAMPLES), #run trimmomatic 
	#QC trimmed reads
	trimmedFastQC = expand(qcDir + "/trimmed_fastqc/{sample}.trimmed_fastqc.zip",sample=SAMPLES), #run fastqc for trimmed reads
	#extract pathogen reads
	extractPathRead = expand(cleanDir + "/{sample}.pathogen.fastq.gz", sample=SAMPLES),

# Paired-ends
elif (seq_type == "pe"):
    #QC raw reads
    rawFastQC = expand(qcDir + "/raw_fastqc/{sample}_{strand}_fastqc.zip", sample=SAMPLES, strand=STRANDS)
    #check host ribosomal RNA
    rRNACheck = expand(trimDir + "/{sample}_" + strand2 + ".filtRNA.fastq.gz", sample=SAMPLES)
	#remove duplicate reads	
    rmDup = expand(trimDir + "/{sample}_{strand}.rmdup.fastq.gz",sample=SAMPLES, strand=STRANDS)
    #remove PhiX174 contaminant
    rmCtm = expand(trimDir + "/{sample}_{strand}.rmdup_ctm.fastq.gz", sample=SAMPLES, strand=STRANDS) #run clumpify to remove duplicate reads
    #Trim reads
    trim = expand(trimDir + "/{sample}_{strand}.trimmed.fastq.gz", sample=SAMPLES, strand=STRANDS), #run trimmomatic 
	#QC trimmed reads
    trimmedFastQC = expand(qcDir + "/trimmed_fastqc/{sample}_{strand}.trimmed_fastqc.zip", sample=SAMPLES, strand=STRANDS), #run fastqc for trimmed reads
    #extract pathogen reads
    extractPathRead = expand(cleanDir + "/{sample}_{strand}.pathogen.fastq.gz", sample=SAMPLES, strand=STRANDS),

else:
    sys.exit("Error: invalid 'seq_type parameter'. Must be 'se' or 'pe'")

#Give the user options to blast against only viral sequences or all NCBI nt and nr databases
runBlast = list()
#Run Blastn againt NCBI viral ref
viralblastn=expand(annotateDir + "/{sample}.blastn.txt", sample=SAMPLES)
#Run Diamond/Blastx againt RVDB		
viralblastx=expand(annotateDir + "/{sample}.blastx.txt", sample=SAMPLES)
runBlast.append(viralblastn)
runBlast.append(viralblastx)
if (blastDbType == "all"):
	#Run Blastn againt NCBI nt		
	allBlastnt = expand(annotateDir + "/{sample}.blastnt.txt", sample=SAMPLES)
	#Run Diamond/Blastx againt NCBI nr		
	allBlastnr = expand(annotateDir + "/{sample}.blastnr.txt", sample=SAMPLES)
	runBlast.append(allBlastnt)
	runBlast.append(allBlastnr)

rule all:
	input: #targets
		#QC raw reads
		rawFastQC,

		#remove host ribosomal RNAs
		rRNACheck,

		#remove duplicate reads
		rmDup,
		
		#remove PhiX174 contaminant
		rmCtm,
		
		#Trim reads
		trim,
		
		#QC trimmed reads
		trimmedFastQC,

		#summary fastqc by multiqc
		qcDir + "/multiqc/raw_multiqc.html", #run multiqc for raw reads fastqc
		qcDir + "/multiqc/trimmed_multiqc.html", #run multiqc for trimmed reads fastqc

		#kraken-krona
		expand(classifyDir + "/{sample}.kraken2.report.html", sample=SAMPLES),
		#extract pathogen reads
		extractPathRead,
		
		#Kaiju-krona
		expand(classifyDir + "/{sample}.kaiju.table.txt", sample=SAMPLES),
		expand(classifyDir + "/{sample}.kaiju_krona.html", sample=SAMPLES),

		#assemble pathogen reads
		expand(assembleDir + "/{sample}/contigs.fasta", sample=SAMPLES),
		expand(qcDir + "/quast/{sample}.quast", sample=SAMPLES),
		#quast multiqc
		qcDir + "/multiqc/quast_multiqc.html", #run multiqc for quast

		#annotate contigs using blastn and blastx
		runBlast,
		#reteive references 
		expand(annotateDir + "/{sample}.refNamen.txt", sample=SAMPLES),
		
		#mapping reads to references
		expand(logDir + "/checkPoint/{sample}.retrieveRef.done", sample=SAMPLES),
		
		#mapping reads to contigs, possible novel virus
		expand(logDir + '/checkPoint/{sample}.map2Contig.done', sample=SAMPLES),

		#blastn consensus against NCBI, if failed, run blastx
		expand(logDir + "/checkPoint/{sample}.blastn.done", sample=SAMPLES),
		
		reportFile,
		htmlReport

	message: "Rule all"
	shell: "echo Job done    `date '+%Y-%m-%d %H:%M'`"			
