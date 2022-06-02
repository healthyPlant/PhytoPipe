"""
    Run Kraken2 against NCBI nt and Kaiju against nr_euk pathogen proteins
"""
__author__ = 'Xiaojun Hu <xiaojun.hu@usda.gov>'

trimDir = config["run_info"]["trim"]
cleanDir = config["run_info"]["clean"]
logDir = config["run_info"]["log"]
krakenDir = config["run_info"]["classify"]
krakenDb = config["krakenDb"]
kaijuDb = config["kaijuDb"]
taxDb = config["taxDb"]
seq_type = config["seq_type"]
strand1 = config["strand1"]
strand2 = config["strand2"]

if (seq_type == "pe"):
    paired_string = '--paired'
else:
    paired_string = ''

#for ambiguous input
ruleorder: extractPathReads_pe > extractPathReads_se

# Define input files
def kraken_inputs(wildcards):
    if (seq_type == "pe"):
        reads = expand(trimDir + "/{sample}_{strand}.trimmed.fastq.gz", strand=[strand1,strand2], sample=wildcards.sample)
    elif (seq_type == "se"):
        reads = trimDir + "/{sample}.trimmed.fastq.gz"
    else:
        sys.exit("Error: invalid sequencing type parameter. Must be 'se' or 'pe'")
    return reads

def kaiju_inputs(wildcards):
    reads = trimDir + "/{sample}.pathogen.fastq.gz"
    if (seq_type == "pe"):
        reads = expand(trimDir + "/{sample}_{strand}.trimmed.fastq.gz", strand=[strand1,strand2], sample=wildcards.sample)
    elif (seq_type == "se"):
        reads = trimDir + "/{sample}.trimmed.fastq.gz"
    else:
        sys.exit("Error: invalid sequencing type parameter. Must be 'se' or 'pe'")
    return reads

#since kraken takes a large memory, it can't be used parallel
rule run_kraken:
	"""
	Reads taxonomy classification by Kraken2  
	"""
	input:
		kraken_inputs
	output:
		report = krakenDir + "/{sample}.kraken2.report.txt",
		readClass = temp(krakenDir + "/{sample}.kraken2.txt")
	message:
		'''--- {wildcards.sample} kraken2'''
	threads: config["number_of_threads"] 
	params:
		param = " --use-names --gzip-compressed --confidence 0.05 ", #--quick --fastq-input
		paired_string = paired_string,
		unclassifiedRead = krakenDir + "/{sample}.kraken2.unclassified"
	resources:
		mem = 200,
		runtime = lambda wildcards, attempt: attempt**2*60*10		
	priority: 5
	group: "classfication"
	log:
		logDir + "/kraken/{sample}.log"
	shell: 
		"""
		unclassifiedReadOut={params.unclassifiedRead}'.fastq'
		# the -z operator checks whether the string is null
		if [ ! -z {params.paired_string} ]; then
			unclassifiedReadOut={params.unclassifiedRead}".R#.fastq"
		fi
		#run Kraken2
		kraken2 --threads {threads} --db {krakenDb} {params.param} {params.paired_string} --output {output.readClass} --report {output.report} --unclassified-out $unclassifiedReadOut {input} &>> {log}
		
		if [ -f {params.unclassifiedRead}".R_2.fastq" ]; then
			mv {params.unclassifiedRead}".R_1.fastq" {krakenDir}"/"{wildcards.sample}"_"{strand1}".kraken2.unclassified.fastq"
			mv {params.unclassifiedRead}".R_2.fastq" {krakenDir}"/"{wildcards.sample}"_"{strand2}".kraken2.unclassified.fastq"
		fi
		"""

rule run_krona:
	"""
	View Kraken report using Krona 
	"""
	input:
		krakenDir + "/{sample}.kraken2.report.txt"
	output:
		reducedReport = temp(krakenDir + "/{sample}.kraken2.reduced.report.txt"),
		report = krakenDir + "/{sample}.kraken2.report.html"
	params:
		" -m 3 -t 5 "
	message:
		'''--- {wildcards.sample} kraken2 to krona'''
	shell:
		"""
		cat {input} | awk "{{if(\$1>=0.01) print}}"  > {output.reducedReport}  #reduce report size
		ktImportTaxonomy {params} -o {output.report} {output.reducedReport}
		"""   

rule extractPathReads_se:
	"""
	Extract pathogen reads and merge them with unclassified reads to form a new pathogen fastq file 
	"""
	input:
		classFile = krakenDir + "/{sample}.kraken2.txt",
		report = krakenDir + "/{sample}.kraken2.report.txt",
		reads = trimDir + "/{sample}.trimmed.fastq.gz",
	output:
		idFile = temp(krakenDir + "/{sample}.pathogen.readId.txt"),
		pathReadFile = temp(cleanDir + "/{sample}.pathogen.fastq"),
		gzFile = cleanDir + "/{sample}.pathogen.fastq.gz"
	params:
		param = " -t 2 10239 4751 4762 --include-children ", #bacteria, virus, fungi, oomycetes taxonomy id #archaea(2157),
		unclassifiedRead = krakenDir + "/{sample}.kraken2.unclassified.fastq"
	message:
		'''--- {wildcards.sample} pathogen reads extraction'''
	shell:
		"""
		python {scripts_dir}/getReadIdFromKraken.py -k {input.classFile} -r {input.report} {params.param} -o {output.idFile}
		seqtk subseq {input.reads} {output.idFile} | cat {params.unclassifiedRead} - > {output.pathReadFile}
		gzip -c {output.pathReadFile} > {output.gzFile} 
		rm {params.unclassifiedRead}  
		"""

rule extractPathReads_pe:
	"""
	Extract pathogen reads and merge them with unclassified reads to form a new pathogen fastq file 
	"""
	input:
		classFile = krakenDir + "/{sample}.kraken2.txt",
		report = krakenDir + "/{sample}.kraken2.report.txt",
		reads1 = trimDir + "/{sample}" + "_" + strand1 + ".trimmed.fastq.gz",
		reads2 = trimDir + "/{sample}" + "_" + strand2 + ".trimmed.fastq.gz",
	output:
		idFile = temp(krakenDir + "/{sample}.pathogen.readId.txt"),
		pathReadFile1 = temp(cleanDir + "/{sample}" + "_" + strand1 + ".pathogen.fastq"),
		pathReadFile2 = temp(cleanDir + "/{sample}" + "_" + strand2 + ".pathogen.fastq"),
		gzFile1 = cleanDir + "/{sample}" + "_" + strand1 + ".pathogen.fastq.gz",
		gzFile2 = cleanDir + "/{sample}" + "_" + strand2 + ".pathogen.fastq.gz"
	params:
		param = " -t 2 10239 4751 4762 --include-children ", #bacteria, virus, fungi, oomycetes taxonomy id #archaea(2157),
		unclassifiedRead1 = krakenDir + "/{sample}" + "_" + strand1 + ".kraken2.unclassified.fastq",
		unclassifiedRead2 = krakenDir + "/{sample}" + "_" + strand2 + ".kraken2.unclassified.fastq"
	message:
		'''--- {wildcards.sample} pathogen reads extraction'''
	shell:
		"""
		python {scripts_dir}/getReadIdFromKraken.py -k {input.classFile} -r {input.report} {params.param} -o {output.idFile}
		seqtk subseq {input.reads1} {output.idFile} | cat {params.unclassifiedRead1} - > {output.pathReadFile1}
		seqtk subseq {input.reads2} {output.idFile} | cat {params.unclassifiedRead2} - > {output.pathReadFile2}
		gzip -c {output.pathReadFile1} > {output.gzFile1}   
		gzip -c {output.pathReadFile2} > {output.gzFile2}
		rm {params.unclassifiedRead1}
		rm {params.unclassifiedRead2}
		"""

rule run_kaiju:
	"""
	Reads taxonomy classification by Kaiju based on NCBI protein database  
	"""
	input:
		reads = kaiju_inputs
	output:
		temp(krakenDir + "/{sample}.kaiju.out")
	params:
		paired_string = paired_string
	message:
		'''--- {wildcards.sample} kaiju'''
	threads: config["number_of_threads"] 
	resources:
		mem = 200,
		runtime = lambda wildcards, attempt: attempt**2*60*10
	priority: 3
	group: "classfication"
	log:
		logDir + "/kaiju/{sample}.log"
	shell:
		"""
		if [ ! -z {params.paired_string} ]; then
			kaiju -z {threads} -t {taxDb}/nodes.dmp -f {kaijuDb} -i {input.reads[0]} -j {input.reads[1]} -o {output} &>> {log}
		else
			kaiju -z {threads} -t {taxDb}/nodes.dmp -f {kaijuDb} -i {input.reads} -o {output} &>> {log} 
		fi
		"""

rule run_kaijuKrona:
	"""
	View Kaiju report using Krona 
	"""
	input:
		krakenDir + "/{sample}.kaiju.out"
	output:
		kaijuTable = krakenDir + "/{sample}.kaiju.table.txt",
		kaijuKrona = krakenDir + "/{sample}.kaiju_krona",
		kaijuKronaHtml = krakenDir + "/{sample}.kaiju_krona.html"
	params:
		" -r species -c 100 -p -e "  #-l superkingdom,phylum,class,order,family,genus,species; -c 100 the minimum required number of reads for the taxon: 100 (except viruses) 
	message:
		'''--- {wildcards.sample} Kaiju to Krona'''
	shell:
		"""
		kaiju2table -t {taxDb}/nodes.dmp -n {taxDb}/names.dmp {params} -o {output.kaijuTable} {input}
		kaiju2krona -t {taxDb}/nodes.dmp -n {taxDb}/names.dmp -i {input} -o {output.kaijuKrona}
		ktImportText -o {output.kaijuKronaHtml} {output.kaijuKrona}
		"""   
