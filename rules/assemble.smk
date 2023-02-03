"""
    Assemble reads using Spades or Trinity. Default assembler is Trinity. 
"""
__author__ = 'Xiaojun Hu <xiaojun.hu@usda.gov>'

trimDir = config["workDir"] + "/" + config["run_info"]["trim"]
cleanDir = config["workDir"] + "/" + config["run_info"]["clean"]
logDir = config["workDir"] + "/" + config["run_info"]["log"]
assembleDir = config["workDir"] + "/" + config["run_info"]["assemble"]
qcDir = config["workDir"] + "/" + config["run_info"]["qc"]
seq_type = config["seq_type"]
samples = config["samples"]
assembler = config["assembler"]
spades_param = config["spades_param"]
mapReadType = config["mapReadType"]

if (seq_type == "pe"):
	paired_string = '--paired'
else:
	paired_string = ''

#print("Assemble samples: ", samples)
# Define input files
def assemble_inputs(wildcards):
	if (mapReadType == 'clean'):
		if (seq_type == "pe"):
			reads = expand(cleanDir + "/{sample}_{strand}.pathogen.fastq.gz", strand=['R1','R2'], sample=wildcards.sample)
		elif (seq_type == "se"):
			reads = cleanDir + "/{sample}.pathogen.fastq.gz"
		else:
			sys.exit("Error: invalid sequencing type parameter. Must be 'se' or 'pe'")
	else:
		if (seq_type == "pe"):
			reads = expand(trimDir + "/{sample}_{strand}.trimmed.fastq.gz", strand=['R1','R2'], sample=wildcards.sample)
		elif (seq_type == "se"):
			reads = trimDir + "/{sample}.trimmed.fastq.gz"
		else:
			sys.exit("Error: invalid sequencing type parameter. Must be 'se' or 'pe'")
	#print(reads)
	return reads

if assembler == "Spades":
	rule run_spades:
		"""
		Assemble reads using SPAdes 
		"""
		input:
			read = assemble_inputs
		output:
			contigs = assembleDir + "/{sample}/contigs.fasta"
		message:
			"""--- assemble {wildcards.sample} genome using SPAdes."""
		params:
			param = spades_param,
			sampleDir = assembleDir + "/{sample}/spades",
			paired_string = paired_string,
		threads: config["number_of_threads"] 
		group: "Spades"
		log:
			assemlog = logDir + "/assemble/{sample}.spades.log",
			normlog = logDir + "/assemble/{sample}.bbnorm.log"         
		resources:
			mem = 200,
			threads = config["number_of_threads"] 			
		shell: #8 white space instead of tab to resolve "Command must be given as string after the shell keyword. "
			"""
			if [ ! -z {params.paired_string} ]; then
				#paired-end
				#bbout1={params.sampleDir}"_"{strand1}".bbnorm.fastq.gz"
				#bbout2={params.sampleDir}"_"{strand2}".bbnorm.fastq.gz"
				#bbnorm.sh in={input.read[0]} in2={input.read[1]} out=$bbout1  out2=$bbout2 threads={threads} target=80 &>> {log.normlog}  #run bbnorm to reduce high depth
				#spades.py --threads {threads} {params.param} -o {params.sampleDir} -1 $bbout1 -2 $bbout2  &>> {log.assemlog}
				spades.py --threads {threads} {params.param} -o {params.sampleDir} -1 {input.read[0]} -2 {input.read[1]}  &>> {log.assemlog}
				cp {params.sampleDir}/contigs.fasta {output.contigs}
				rm -rf {params.sampleDir}
			else
				#single-end
				#bbout={params.sampleDir}".bbnorm.fastq.gz"
				#bbnorm.sh in={input} out=$bbout threads={threads} target=80 &>> {log.normlog}  #run bbnorm to reduce high depth
				#spades.py --threads {threads} {params.param} -o {params.sampleDir} -s $bbout &>> {log.assemlog}
				spades.py --threads {threads} {params.param} -o {params.sampleDir} -s {input} &>> {log.assemlog}
				cp {params.sampleDir}/contigs.fasta {output.contigs}
				rm -rf {params.sampleDir}
			fi
			"""
else:
	#Trinity assembler is good for low titer viruses
	rule run_trinity:
		"""
		Assemble reads using Trinity 
		"""
		input:
			read = assemble_inputs 
		output:
			contigs = assembleDir + "/{sample}/contigs.fasta" 
		message:
			"""--- assemble {wildcards.sample} genome using Trinity."""
		params:
			param = "--seqType fq --full_cleanup", 
			outdir= assembleDir + "/{sample}/trinity",
			paired_string = paired_string
		threads: config["number_of_threads"] 
		group: "Trinity"
		log:
			assemlog = logDir + "/assemble/{sample}.trinity.log"
		resources:
			mem = 200,
			threads = config["number_of_threads"] 			
		shell:
			"""
			if [ ! -z {params.paired_string} ]; then
				#paired-end
				Trinity {params.param} --CPU {threads} --max_memory {resources.mem}G --left {input.read[0]} --right {input.read[1]} --output {params.outdir} > {log.assemlog}
			else
				#single-end		  
				Trinity {params.param} --CPU {threads} --single {input} --max_memory {resources.mem}G --output {params.outdir} > {log.assemlog}
			fi
			mv {params.outdir}.Trinity.fasta {output.contigs}
			rm {params.outdir}.Trinity.fasta.gene_trans_map
			"""

rule run_quast:
	"""
	Perform genome assembly quality assesment using QUAST
	"""
	input:
		assembleDir + "/{sample}/contigs.fasta"
	output:
		directory(qcDir + "/quast/{sample}.quast")
	message:
		'''--- {wildcards.sample} genome assembly quality assesment using QUAST.'''
	log:
		logDir + "/quast/{sample}.log"
	shell:
	  "quast.py {input} -m 100 -o {output} 2> {log} 1>&2"

rule multiqc_quast:
	"""
	Summary assembly quality result by multiqc
	"""
	input:
		expand(qcDir + "/quast/{smp}.quast", smp=samples)
	output:
		html = qcDir + "/multiqc/quast_multiqc.html",
		stats = qcDir + "/multiqc/quast_multiqc_data/multiqc_general_stats.txt"
	message:
		'''--- Summary quast with multiqc.'''
	log:
		logDir + "/multiqc/quast_multiqc.log"
	params:
		name = "quast_multiqc.html",
		inputDir = qcDir + "/quast",
		outputDir = qcDir + "/multiqc"    
	priority: -10
	shell:
		"""
		multiqc -f -n {params.name} -o {params.outputDir} {params.inputDir} 2> {log} 1>&2  
		"""
