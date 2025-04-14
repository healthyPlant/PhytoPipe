"""
    Snakemake main file 
"""
__author__ = 'Xiaojun Hu <xiaojun.hu@usda.gov>'

import os,subprocess,itertools,glob,re,sys
import datetime, time
from snakemake.utils import report

start = time.time()

#SETUP PROCESS-RELATED CONFIGURATION 
try:
	CONFIGFILE = str(workflow.overwrite_configfiles[0])
except:
	CONFIGFILE = str(workflow.overwrite_configfile[0])

#setup workdir dynamically
#rules will then be expressed relative to the workdir
workdir: config["workDir"]  #can be changed at the snakemake command line
fastqDir = config["fastqDir"]  #can be changed at the snakemake command line
flowCellDir = config["flowCellDir"]  #can be changed at the snakemake command line

#Get relative position of Snakefile from wd
SNAKEFILE = workflow.snakefile
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)
CORES = workflow.cores #note workflow.cores work only on version 5.11.2 above
config["number_of_threads"] = CORES  #change threads number by command paramter --cores/-j 16
config["snakefile"] = SNAKEFILE
config["snakefile_dir"] = SNAKEFILE_DIR

#print(workflow.__dict__) #print all attributes
#Establish snakefile and environment dictionaries
snakefiles_dir = os.path.abspath(os.path.join(SNAKEFILE_DIR, "rules"))		#directory for additional snakefiles
environments_dir = os.path.abspath(os.path.join(SNAKEFILE_DIR, "envs"))	#directory for conda environment yaml files
scripts_dir = os.path.abspath(os.path.join(SNAKEFILE_DIR, "scripts")) 			#directory for extra scripts used in workflow

#running folders
rawReadDir = config["workDir"] + "/" + config["run_info"]["raw"]
qcDir = config["workDir"] + "/" + config["run_info"]["qc"]
logDir = config["workDir"] + "/" + config["run_info"]["log"]
trimDir = config["workDir"] + "/" + config["run_info"]["trim"]
assembleDir = config["workDir"] + "/" + config["run_info"]["assemble"]
classifyDir = config["workDir"] + "/" + config["run_info"]["classify"]
annotateDir = config["workDir"] + "/" + config["run_info"]["annotate"]
mapDir = config["workDir"] + "/" + config["run_info"]["map"]
reportDir = config["workDir"] + "/" + config["run_info"]["report"]
novelDir = config["workDir"] + "/" + config["run_info"]["novel"]
blastDbType = config["blastDbType"]
seq_type = config["seq_type"]
reportFile = reportDir + "/report.txt" 
htmlReport = reportFile.replace(".txt", ".html")
strand1 = config["strand1"]
strand2 = config["strand2"]
STRANDS = [config["strand1"], config["strand2"]]

## Check whether input fastq files are compressed or not
if  "input_format" in config.keys():
    input_format = config["input_format"]
else:
    input_format = "fastq.gz"

#set up soft link for fastq files if fastq folder is surpplied
if fastqDir:
	if not os.path.exists(rawReadDir):
		os.system("ln -sf " + fastqDir + " " + rawReadDir)
	# Check if given path is link
	if os.path.exists(rawReadDir) and os.path.islink(rawReadDir):
		print(fastqDir, " soft link is built." )
	else:
		# soft link broken, copy files
		if not os.path.isdir(rawReadDir): #check a folder exists
			os.mkdir(rawReadDir) #make a folder
		print("cp " + fastqDir + "/*." + input_format + " " + rawReadDir + "/")  
			
if not os.listdir(rawReadDir): 
    print(rawReadDir + " is empty. Exit!") 
    sys.exit()

#get sample names from rawReadDir
SAMPLES = [os.path.basename(x) for x in glob.glob(rawReadDir + '/*.' +  input_format)]
#remove 'Undetermined_*'
SAMPLES = [x for x in SAMPLES if not x.startswith('Undetermined_')]

#remove '_R1_001.fastq.gz'
SAMPLES = [x.replace("." + input_format,"") for x in SAMPLES]
if seq_type == 'se' and config["strand1"] != '':
	SAMPLES = [re.sub(r"(.*)_%s" % config["strand1"], "\\1", x) for x in SAMPLES]

if seq_type == 'pe':
	SAMPLES = [re.sub(rf"_[{strand1}|{strand2}]$", "", x) for x in SAMPLES]

SAMPLES = list(set(SAMPLES))  #remove duplicated from a list using set
print("Samples: %s"  % ", ".join(SAMPLES))
config["samples"] = SAMPLES  #save sample names in a global viriable

if len(SAMPLES) == 0:
	sys.exit("No sample files found. Exit!")

#Run rules
include: os.path.join(snakefiles_dir, "cleanReads.smk")
include: os.path.join(snakefiles_dir, "classifyReads.smk")
include: os.path.join(snakefiles_dir, "assemble.smk")
include: os.path.join(snakefiles_dir, "annotate.smk")
include: os.path.join(snakefiles_dir, "mappingReads.smk")
include: os.path.join(snakefiles_dir, "makeReport.smk")

#the following two functions: gather_samples and gather_blast_input for handling assembly failed samples
#If no contigs, the downstream rules will not be executed for the sample
def gather_samples(wildcards):
	"""
	Aggregate samples having a successful assembly based on checkpoint extract_contig in annotate.smk
	"""
	samples = []
	for sample in SAMPLES:
		fn = checkpoints.extract_contig.get(sample=sample).output[0]
		if os.path.exists(fn) and os.stat(fn).st_size > 0:
			samples.append(sample)
	return samples

def gather_blast_input(wildcards):
	"""
	Aggregate samples for blast 
	"""
	samples = gather_samples(wildcards)
	if len(samples) == 0:
		sys.exit("No samples are sucessfully assembled. Exit!")
	else:
		print("Sucessfully assembled samples: ", samples)

	#Give the user options to blast against only viral sequences or all NCBI nt and nr databases
	runBlast = list()
	#Run Blastn againt NCBI viral ref
	viralblastn = expand(annotateDir + "/{sample}.blastn.txt", sample=samples)
	#Run Diamond/Blastx againt RVDB
	viralblastx = expand(annotateDir + "/{sample}.blastx.txt", sample=samples)
	runBlast.append(viralblastn)
	runBlast.append(viralblastx)
	if (blastDbType == "all"):
		#Run Blastn againt NCBI nt
		allBlastnt = expand(annotateDir + "/{sample}.blastnt.krona.html", sample=samples)
		#Run Diamond/Blastx againt NCBI nr
		allBlastnr = expand(annotateDir + "/{sample}.blastnr.krona.html", sample=samples)
		runBlast.append(allBlastnt)
		runBlast.append(allBlastnr)

	return runBlast

# a target rule to define the desired final output
rule all:
	input: #targets
		htmlReport,
		unpack(gather_blast_input),

	message: "Rule all"
	shell: "echo Job done    `date '+%Y-%m-%d %H:%M'`"

# The onsuccess handler is executed if the workflow finished without error.
onsuccess:
	#check contigs, if no conitgs, report assembly fail
	assembleFail = []
	for sample in SAMPLES:
		filename = assembleDir + "/" + sample + "/contigs.fasta"
		if not os.path.exists(filename) or os.stat(filename).st_size == 0:
			assembleFail.append(sample)
	failNum = len(assembleFail)

	print("\n############################################################################")
	if failNum > 0:
		print("# Aassembly failed sample(s): ", ", ".join(assembleFail))
	else:
		print("# PhytoPipe finished without errors ")
	print("# Plase check the result in the folder " + reportDir)
	print("############################################################################")
	print("\nRunning time in minutes: %s\n" % round((time.time() - start)/60,1))

# Else, the onerror handler is executed.
onerror:
	print("\n\n#####################\n# An error occurred #\n#####################\n\n")
	command = "grep -r 'Error:' " + logDir  
	try:
		output = subprocess.check_output(command, shell=True, text=True, stderr=subprocess.STDOUT)
		print(output)
	except subprocess.CalledProcessError as e:
		print(f"Error executing command: {e}")
		print("No errors found in the log files.")
	print("############################################################################")
	print("\nRunning time in minutes: %s\n" % round((time.time() - start)/60,1))
	#shell("mail -s "an error occurred" youremail@provider.com < {log}")

# onstart handler will be executed before the workflow starts. Note that dry-runs do not trigger any of the handlers
onstart:
	print("Running PhytoPipe for samples : %s"  % ", ".join(SAMPLES))
