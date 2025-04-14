"""
    Raw reads are cleaned by removing duplicates and control PhiX174, and trimmed using trimmomatic
    Raw reads fastQC and trimmed read fastQC are reported and summrized by multiQC
"""
__author__ = 'Xiaojun Hu <xiaojun.hu@usda.gov>'

workDir = config["workDir"]
rawReadDir = config["workDir"] + "/" + config["run_info"]["raw"]
qcDir = config["workDir"] + "/" + config["run_info"]["qc"]
logDir = config["workDir"] + "/" + config["run_info"]["log"]
trimDir = config["workDir"] + "/" + config["run_info"]["trim"]
reportDir = config["workDir"] + "/" + config["run_info"]["report"]
adapters = config["adapters"]
trimmomatic = config["trimmomatic"]
control = config["control"]
input_format = config["input_format"]
rRNA_db = config["euk_rRNA"]
seq_type = config["seq_type"]
clumpify_param = config["clumpify_param"]
trimmomatic_param = config["trimmomatic_param"]
strand1 = config["strand1"]
strand2 = config["strand2"]
samples = config["samples"]
#print(samples)
rawFiles = glob.glob(rawReadDir + "/*." + input_format)
#true001 = 0


#to prevent snakemake from interpreting wildcard strings as regular expressions. For this, escape each wildcard value using wildcard_constraints.
wildcard_constraints:
sample= '|'.join([re.escape(x) for x in samples]) #add escape (\) before special character and join sample names using "|"

# Define input files
def raw_inputs(wildcards):
	#print({wildcards.sample}, {wildcards.strand})
	if (seq_type == "pe"):
		return expand(rawReadDir + "/{sample}_{strand}." + input_format,  strand=[strand1,strand2], sample=wildcards.sample) 
	elif (seq_type == "se" ):
		if strand1:
			return rawReadDir + "/{sample}_" + strand1 + "." + input_format
		else:
			return rawReadDir + "/{sample}." + input_format
	else:
		sys.exit("Error: invalid sequencing type parameter. Must be 'se' or 'pe'")

def fastqc_raw_inputs(wildcards): 
	if (seq_type == "pe"):
		return rawReadDir + "/{sample}_{strand}." + input_format #,  strand=wildcards.strand, sample=wildcards.sample
	elif (seq_type == "se" ):
		if strand1:
			return rawReadDir + "/{sample}_" + strand1 + "." + input_format
		else:
			return rawReadDir + "/{sample}." + input_format
	else:
		sys.exit("Error: invalid sequencing type parameter. Must be 'se' or 'pe'")

# Define raw multiQC input files
def raw_multiqc_inputs(wildcards):
	if (seq_type == "pe"):
		qcOut = expand(qcDir + "/raw_fastqc/{sample}_{strand}_fastqc.zip", strand=[strand1,strand2], sample=samples)
	elif (seq_type == "se"):
		qcOut = expand(rules.raw_fastqc_se.output, sample=samples) #ensure the rule can only be started after fastqc processes all samples
	else:
		sys.exit("Error: invalid sequencing type parameter. Must be 'se' or 'pe'")
	return qcOut

# Define trimmed multiQC input files
def trimmed_multiqc_inputs(wildcards):
	if (seq_type == "pe"):
		qcOut = expand(qcDir + "/trimmed_fastqc/{sample}_{strand}.trimmed_fastqc.zip", strand=['R1','R2'], sample=samples)
	elif (seq_type == "se"):
		qcOut = expand(rules.trim_fastqc_se.output, sample=samples) #ensure the rule can only be started after fastqc processes all samples
	else:
		sys.exit("Error: invalid sequencing type parameter. Must be 'se' or 'pe'")
	return qcOut

if (seq_type == "se"):
	rule raw_fastqc_se:
		"""
		Perform a single-end raw reads quality check using fastqc
		"""
		input:
			fastqc_raw_inputs   
		output:
			#FastQC names output files by taking the original filename, removing the .fastq extension, and replacing the extension with _fastqc.html or _fastqc.zip.
			html = qcDir + "/raw_fastqc/{sample}_fastqc.html",
			zip = qcDir + "/raw_fastqc/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
		message:
			'''--- {wildcards.sample} raw reads quality check with fastqc.'''
		log:
			logDir + "/raw_fastqc/{sample}.log"
		threads: 4
		params: 
			strand1 = config["strand1"]
		shell:
			"""
			fastqc -t {threads} --quiet --outdir {qcDir}/raw_fastqc {input} >> {log} 2>&1
			if [[ -s qc/raw_fastqc/{wildcards.sample}_{params.strand1}_001_fastqc.html ]]; then
				mv qc/raw_fastqc/{wildcards.sample}_{params.strand1}_001_fastqc.html  {output.html}  #rename output file to fit multiQC input
				mv qc/raw_fastqc/{wildcards.sample}_{params.strand1}_001_fastqc.zip  {output.zip}
			fi
			if [[ -s qc/raw_fastqc/{wildcards.sample}_{params.strand1}_fastqc.html ]]; then
				mv qc/raw_fastqc/{wildcards.sample}_{params.strand1}_fastqc.html  {output.html}  #rename output file to fit multiQC input
				mv qc/raw_fastqc/{wildcards.sample}_{params.strand1}_fastqc.zip  {output.zip}
			fi
			"""

	rule check_rRNA_se:
		"""
		Perform a raw single reads ribosomal RNA check using bbduk.sh
		"""
		input:
			reads = raw_inputs
		output:
			log = logDir + "/rRNA_qc/{sample}.log",
			filtRNA = temp(trimDir + "/{sample}.filtRNA.fastq.gz"),
			rRNA = temp(trimDir + "/{sample}.ribo.fastq.gz")
		message:
			'''--- {wildcards.sample} ribosomal RNAs check with bbduk.sh.'''
		threads:
			config.get("number_of_threads", 1)
		resources:
			mem=12 #mem_mb=4000
		params: 
			kmer = "k=31"
		shell:
			"""
			bbduk.sh -Xmx10g in={input.reads} out={output.filtRNA} outm={output.rRNA} {params.kmer} ref={rRNA_db} 2>{output.log}
			"""
	rule remove_duplicate_se:
		"""
		Remove (PCR)-duplicated single-end reads using clumpify.sh in BBMap
		"""
		input:
			read = trimDir + "/{sample}.filtRNA.fastq.gz"
		output:
			temp(trimDir + "/{sample}.rmdup.fastq.gz")
		message:
			'''--- Remove {wildcards.sample} duplicate reads with clumpify.'''
		log:
			logDir + "/removeDuplicate/{sample}.log"
		resources:
			mem = 12,  # its 85% = Xmx10g
			threads = config.get("number_of_threads", 1)
		params:
			param = clumpify_param  
		group: "rmDup"
		shell:
			"""
			clumpify.sh -Xmx10g overwrite=true in={input.read} out={output} {params.param} threads={threads} 2> {log}
			"""

	rule remove_control_se:
		"""
		Remove Illimina control ex. PhiX174
		"""
		input:
			trimDir + "/{sample}.rmdup.fastq.gz"
		output:
			cleaned = temp(trimDir + "/{sample}.rmdup_ctm.fastq.gz"),
			removed = temp(trimDir + "/{sample}.ctm.fastq")  #if you want keep control reads, remove temp() 
		message:
			'''--- Remove {wildcards.sample} control, like PhiX174 with bbsplit.'''
		resources:
			mem = 12 #mem_mb=4000
		log:
			logDir + "/removeControl/{sample}.log"
		shell:
			"""
			#To index:  bbsplit.sh ref={control} path=/ppq/data2/resources/phiX174
			bbsplit.sh -Xmx10g ref={control} in={input} out_phi-X174={output.removed} outu={output.cleaned} 2> {log}  #automatically generate reference (folder name: ref) in the current work directory
			"""

	rule trimmomatic_se:
		"""
		Trim adaptors and low quality reads by Trimmomatic
		"""
		input:
			trimDir + "/{sample}.rmdup_ctm.fastq.gz"
		output:
			trimDir + "/{sample}.trimmed.fastq.gz"
		message:
			'''--- Trim {wildcards.sample} reads with trimmomatic.'''
		log:
			logDir + "/trimmomatic/{sample}.log"
		resources:
			mem = 20 
		params: 
			trimmer= trimmomatic_param,
			adapters = adapters
		threads: 4
		shell:
			"""
			java -jar {trimmomatic} SE -threads {threads} -phred33 {input} {output} ILLUMINACLIP:{params.adapters}:2:30:10 {params.trimmer} 2> {log}
			"""

	rule trim_fastqc_se:
		"""
		Perform a trimmed single-end reads quality check using fastqc
		"""
		input:
			trimDir + "/{sample}.trimmed.fastq.gz"
		output:
			html = qcDir + "/trimmed_fastqc/{sample}.trimmed_fastqc.html",
			zip = qcDir + "/trimmed_fastqc/{sample}.trimmed_fastqc.zip" 
		message:
			'''--- {wildcards.sample} trimmed reads quality check with fastqc.'''
		log:
			logDir + "/trimmed_fastqc/{sample}.log"
		threads: 4
		shell:
			"""
			fastqc -t {threads} --quiet --outdir {qcDir}/trimmed_fastqc {input} >> {log} 2>&1
			"""


if (seq_type == "pe"):
	rule raw_fastqc_pe:
		"""
		Perform a paired-end raw reads quality check using fastqc
		"""
		input:
			fastqc_raw_inputs 
		output:
			#FastQC names output files by taking the original filename, removing the .fastq extension, and replacing the extension with _fastqc.html or _fastqc.zip.
			html = qcDir + "/raw_fastqc/{sample}_{strand}_fastqc.html",
			zip = qcDir + "/raw_fastqc/{sample}_{strand}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
		message:
			'''--- {wildcards.sample}_{wildcards.strand} raw reads quality check with fastqc.'''
		log:
			logDir + "/raw_fastqc/{sample}_{strand}.log"
		threads: 4
		shell:
			"""
			fastqc -t {threads} --quiet --outdir {qcDir}/raw_fastqc {input} >> {log} 2>&1
			echo {input}
			if [[ -s qc/raw_fastqc/{wildcards.sample}_{wildcards.strand}_001_fastqc.html ]]; then
				mv qc/raw_fastqc/{wildcards.sample}_{wildcards.strand}_001_fastqc.html  {output.html}  #rename output file to fit multiQC input
				mv qc/raw_fastqc/{wildcards.sample}_{wildcards.strand}_001_fastqc.zip  {output.zip}
			fi
			"""

	rule check_rRNA_pe:
		"""
		Perform a raw paired-reads ribosomal RNA check using bbduk.sh
		"""
		input:
			reads = raw_inputs
		output:
			log = logDir + "/rRNA_qc/{sample}.log",
			filtRNA1 = temp(trimDir + "/{sample}" + "_R1.filtRNA.fastq.gz"),
			rRNA1 = temp(trimDir + "/{sample}" + "_R1.ribo.fastq.gz"), 
			filtRNA2 = temp(trimDir + "/{sample}" + "_R2.filtRNA.fastq.gz"),
			rRNA2 = temp(trimDir + "/{sample}" + "_R2.ribo.fastq.gz")
		message:
			'''--- {wildcards.sample} ribosomal RNAs check with bbduk.sh.'''
		threads:
			config.get("number_of_threads", 1)
		resources:
			mem=12 #mem_mb=4000
		params: 
			kmer = "k=31"
		shell:
			"""
			bbduk.sh -Xmx10g in={input.reads[0]} in2={input.reads[1]} out={output.filtRNA1} out2={output.filtRNA2} outm={output.rRNA1} outm2={output.rRNA2} {params.kmer} ref={rRNA_db} 2>{output.log}
			"""

	rule remove_duplicate_pe:
		"""
		Remove (PCR)-duplicated paired-end reads using clumpify.sh in BBMap
		"""
		input:
			read1 = trimDir + "/{sample}" + "_R1.filtRNA.fastq.gz",
			read2 = trimDir + "/{sample}" + "_R2.filtRNA.fastq.gz"
		output:
			rmdup1 = temp(trimDir + "/{sample}" + "_R1.rmdup.fastq.gz"),
			rmdup2 = temp(trimDir + "/{sample}" + "_R2.rmdup.fastq.gz")
		message:
			'''--- Remove {wildcards.sample} duplicate reads with clumpify.'''
		log:
			logDir + "/removeDuplicate/{sample}.log"
		resources:
			mem = 12,  # its 85% = Xmx10g
			threads = config.get("number_of_threads", 1)
		params: 
			param = clumpify_param  
		group: "rmDup"
		shell:
			"""
			clumpify.sh -Xmx10g overwrite=true in={input.read1} in2={input.read2} out={output.rmdup1} out2={output.rmdup2} {params.param} threads={threads} 2> {log}
			"""

	rule remove_control_pe:
		"""
		Remove Illimina control ex. PhiX174, treated as a control
		"""
		input:
			read1 = trimDir + "/{sample}" + "_R1.rmdup.fastq.gz",
			read2 = trimDir + "/{sample}" + "_R2.rmdup.fastq.gz"
		output:
			cleaned1 = temp(trimDir + "/{sample}" + "_R1.rmdup_ctm.fastq.gz"),
			cleaned2 = temp(trimDir + "/{sample}" + "_R2.rmdup_ctm.fastq.gz"),
			removed = temp(trimDir + "/{sample}.ctm.fastq") #if you want keep control reads, remove temp() 
		message:
			'''--- Remove {wildcards.sample} control ex. PhiX174 with bbsplit.'''
		threads:
			config.get("number_of_threads", 1)
		resources:
			mem=12 #mem_mb=4000
		log:
			logDir + "/removeControl/{sample}.log"
		shell:
			"""
			#To index:  bbsplit.sh ref={control} path=/ppq/data2/resources/phiX174
			bbsplit.sh -Xmx10g ref={control} in={input.read1} in2={input.read2} out_phi-X174={output.removed} outu={output.cleaned1} outu2={output.cleaned2} 2> {log}  #automatically generate reference (folder name: ref) in the current work directory
			#basename="trimDir/%_R#.fastq.gz" 
			"""

	rule trimmomatic_pe:
		"""
		Trim adaptors and low quality reads by Trimmomatic
		"""
		input:
			read1 = trimDir + "/{sample}" + "_R1.rmdup_ctm.fastq.gz",
			read2 = trimDir + "/{sample}" + "_R2.rmdup_ctm.fastq.gz",
		output:
			trimmed1 = trimDir + "/{sample}" + "_R1.trimmed.fastq.gz",
			trimmed2 = trimDir + "/{sample}" + "_R2.trimmed.fastq.gz",
			unpaired1 = trimDir + "/{sample}" + "_R1.unpaired.fastq.gz",
			unpaired2 = trimDir + "/{sample}" + "_R2.unpaired.fastq.gz",
		message:
			'''--- Trim {wildcards.sample} reads with trimmomatic.'''
		log:
			logDir + "/trimmomatic/{sample}.log"
		resources:
			mem = 20 
		params: 
			trimmer= trimmomatic_param, 
			adapters = adapters
		threads: 4
		shell:
			"""
			java -jar {trimmomatic} PE -threads {threads}  -phred33 {input.read1} {input.read2} {output.trimmed1} {output.unpaired1} {output.trimmed2} {output.unpaired2} ILLUMINACLIP:{params.adapters}:2:30:10 {params.trimmer} 2> {log}
			"""

	rule trim_fastqc_pe:
		"""
		Perform a trimmed paired-end reads quality check using fastqc
		"""
		input:
			trimDir + "/{sample}_{strand}.trimmed.fastq.gz"
		output:
			html = qcDir + "/trimmed_fastqc/{sample}_{strand}.trimmed_fastqc.html",
			zip = qcDir + "/trimmed_fastqc/{sample}_{strand}.trimmed_fastqc.zip" 
		message:
			'''--- {wildcards.sample} {wildcards.strand} trimmed reads quality check with fastqc.'''
		log:
			logDir + "/trimmed_fastqc/{sample}_{strand}.log"
		threads: 4
		shell:
			"""
			fastqc -t {threads} --quiet --outdir {qcDir}/trimmed_fastqc {input} >> {log} 2>&1
			"""

rule raw_multiqc:
	"""
	Summary raw redas fastqc result by multiqc
	"""
	input:
		raw_multiqc_inputs
	output:
		html = qcDir + "/multiqc/raw_multiqc.html",
		stats = qcDir + "/multiqc/raw_multiqc_data/multiqc_general_stats.txt"
	message:
		'''--- Summary raw reads fastqc using multiqc.'''
	log:
		logDir + "/multiqc/raw_multiqc.log"
	params:
		name = "raw_multiqc.html",
		inputDir = qcDir + "/raw_fastqc",
		outputDir = qcDir + "/multiqc"
	priority: -10
	shell:
		"""
		multiqc -f -n {params.name} -o {params.outputDir} {params.inputDir} 2> {log}
		"""

rule trim_multiqc:
	"""
	Summary trimmed redas fastqc result by multiqc
	"""
	input:
		trimmed_multiqc_inputs 
	output:
		html = qcDir + "/multiqc/trimmed_multiqc.html",
		stats = qcDir + "/multiqc/trimmed_multiqc_data/multiqc_general_stats.txt"
	message:
		'''--- Summary trimmed reads fastqc with multiqc.'''
	log:
		logDir + "/multiqc/trimmed_multiqc.log"
	params:
		name = "trimmed_multiqc.html",
		inputDir = qcDir + "/trimmed_fastqc",
		outputDir = qcDir + "/multiqc"
	priority: -15
	shell:
		"""
		multiqc -f -n {params.name} -o {params.outputDir} {params.inputDir} &>> {log}
		"""
