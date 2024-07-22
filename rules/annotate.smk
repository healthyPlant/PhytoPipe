"""
	Annotate contigs by Blastn against NCBI viral refSeq and Diamond blastx against RVDB (only for viruses)
	Blastn against NCBI nt and Diamond blastx against NCBI nr are also supported.
"""
__author__ = 'Xiaojun Hu <xiaojun.hu@usda.gov>'

annotateDir = config["workDir"] + "/" + config["run_info"]["annotate"]
krakenDir = config["workDir"] + "/" + config["run_info"]["classify"]
assembleDir = config["workDir"] + "/" + config["run_info"]["assemble"]
logDir = config["workDir"] + "/" + config["run_info"]["log"]
blastndb = config["blastnViralDb"]
blastxdb = config["blastxViralDb"]
blastxTaxondb = config["blastxViralTaxonDb"] 
blastnTaxondb = config["blastnViralTaxonDb"] 
nt = config["blastnDb"]  #NCBI blastdb nt
nr = config["blastxDb"]  #diamond format NCBI nr
microbialTaxonIds = config["microbialTaxon"]
blastEvalue = config["blastEvalue"]

#http://www.metagenomics.wiki/tools/blast/blastn-output-format-6
#download ncbi nr fasta sequences
#wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
#make diamond nr database
#diamond makedb --in nr.gz -d nr
#download ncbi blastn db
#ftp://ftp.ncbi.nlm.nih.gov/blast/db/${db}.tar.gz

rule run_blastn:
	"""
	Annotate contigs by blastn
	"""
	input:
		assembleDir + "/{sample}/contigs.fasta"
	output:
		contigL200 = annotateDir + "/{sample}.contigs.l200.fasta",
		blastnOut0 = temp(annotateDir + "/{sample}.blastn0.txt"),
		blastnOut = annotateDir + "/{sample}.blastn.txt",  #blast the viral reference
		blastnKrona = annotateDir + "/{sample}.blastn.krona.html"
	message:
		'''--- annotatate {wildcards.sample} contigs using blastn.'''
	params:
		" -word_size 11 -max_target_seqs 5 -evalue " + blastEvalue + " -outfmt \"6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore stitle \" " #scomnames staxids sscinames sskingdoms
	threads: 8  #threads > 4, performance do not improve lots
	log:
		logDir + "/annotate/{sample}.blastn.log",
	shell:
		"""
		seqtk seq -L 200 {input} > {output.contigL200}   #only use length >200bp contigs
		blastn -db {blastndb} -num_threads {threads} -query {output.contigL200} -out {output.blastnOut0} {params} &>> {log} 
		python {scripts_dir}/setBlastTaxon.py -b {output.blastnOut0} -t {blastnTaxondb} -o {output.blastnOut} 
		ktImportBLAST -o {output.blastnKrona} {output.blastnOut0}
		"""

rule run_blastx:
	"""
	Annotate contigs by blastx through diamond
	"""
	input:
		annotateDir + "/{sample}.contigs.l200.fasta"
	output:
		blastxOut0 = temp(annotateDir + "/{sample}.blastx0.txt"),
		blastxOut1 = temp(annotateDir + "/{sample}.blastx1.txt"),
		blastxOut = annotateDir + "/{sample}.blastx.txt",
		blastxKrona = annotateDir + "/{sample}.blastx.krona.html"
	message:
		'''--- annotatate {wildcards.sample} contigs using blastx.'''
	params:
		" --max-target-seqs 1 -e  " + blastEvalue + "  --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore stitle "
	threads: 8
	log:
		logDir + "/annotate/{sample}.blastx.log",
	shell:
		"""
		diamond blastx -d {blastxdb} -q {input} -p {threads} -o {output.blastxOut0} {params} &>> {log} 
		python {scripts_dir}/setBlastTaxon.py -b {output.blastxOut0} -t {blastxTaxondb} -o {output.blastxOut} 
		#limit Blast out to ktImportBLAST
		awk -F "\t" 'OFS="\t" {{ split($2,a,"|");$2=a[5]; if($11 < 1e-100) print }}' {output.blastxOut0} | sort -k 11,11g | head -10000 > {output.blastxOut1}
		ktImportBLAST -o {output.blastxKrona} {output.blastxOut1}
		"""

rule run_blastnt:
	"""
	Annotate contigs by blastn againt NCBI nt
	"""
	input:
		contigs = annotateDir + "/{sample}.contigs.l200.fasta",
		microbialTids = microbialTaxonIds
	output:
		blastntOut = annotateDir + "/{sample}.blastnt0.txt"  #blast NCBI all nt
	message:
		'''--- annotatate {wildcards.sample} contigs using blastn against NCBI nt.'''
	params:
		" -max_target_seqs 1 -evalue  " + blastEvalue + " -outfmt \"6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore stitle \" " #-task megablast, megablast faster than blast for intraspecies comparison as it uses a large word-size of 28bp, but it misses something
	threads: 8 
	log:
		logDir + "/annotate/{sample}.blastnt.log"
	shell:
		"""
		blastn -db {nt} -num_threads {threads} -query {input.contigs} -out {output.blastntOut} -taxidlist {input.microbialTids} {params} &>> {log} || true 
		"""

rule summarize_blastnt:
	"""
	Annotate contigs by blastn againt NCBI nt
	"""
	input:
		contigs = annotateDir + "/{sample}.contigs.l200.fasta",
		krakenReport = krakenDir + "/{sample}.kraken2.report.txt",
		blastntOut0 = annotateDir + "/{sample}.blastnt0.txt",
	output:
		blastntOut1 = temp(annotateDir + "/{sample}.blastnt1.txt"),
		acc_list = temp(annotateDir + "/{sample}.blastnt.acc_list.txt"),
		acc_ti = temp(annotateDir + "/{sample}.blastnt.acc_ti.txt"),
		blastntOut = annotateDir + "/{sample}.blastnt.txt",  #blast NCBI all nt
		blastntKrona = annotateDir + "/{sample}.blastnt.krona.html",
		blastntSummary = annotateDir + "/{sample}.blastnt.summary.txt"
	message:
		'''--- summarize {wildcards.sample} contigs blastnt result.'''
	log:
		logDir + "/annotate/{sample}.summarizeBlastnt.log"
	shell:
		"""
		#get taxon id (ti) for subjects in blastnt result
		cut -f2 {input.blastntOut0} | uniq > {output.acc_list}
		#Add " || true " to output errors to the output file, otherwise bash strict mode cause exit
		blastdbcmd -db {nt} -entry_batch {output.acc_list} -outfmt "%i %T"  > {output.acc_ti} || true 
		#add kraken2 report to the blastnt result
		python {scripts_dir}/mergeBlastnKraken.py -b {input.blastntOut0} -k {input.krakenReport} -a {output.acc_ti} -c {input.contigs} -o {output.blastntOut}
		python {scripts_dir}/summarizeBlast.py {output.blastntOut} {output.blastntSummary} 
		#classify blast
		awk -F "\\t" 'OFS="\\t" {{ if($11 < 1e-20) print }}' {input.blastntOut0} | sort -k 11,11g  | head -10000 > {output.blastntOut1} || true 
		ktImportBLAST -o {output.blastntKrona} {output.blastntOut1} &>> {log} || true 
		"""

rule run_blastnr:
	"""
	Annotate contigs by blastx through diamond againt NCBI nr
	"""
	input:
		contigs = annotateDir + "/{sample}.contigs.l200.fasta",
	output:
		blastnrOut = annotateDir + "/{sample}.blastnr0.txt",
	message:
		'''--- annotatate {wildcards.sample} contigs using blastx against NCBI nr'''
	params:
		" --max-target-seqs 1 -e " + blastEvalue + " --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore stitle "
	threads: workflow.cores
	log:
		logDir + "/annotate/{sample}.blastnr.log"
	shell:
		"""
		diamond blastx -d {nr} -q {input.contigs} -p {threads} -o {output.blastnrOut} {params} &>> {log}
		"""

rule summarize_blastnr:
	"""
	Annotate contigs by blastx through diamond againt NCBI nr
	"""
	input:
		contigs = annotateDir + "/{sample}.contigs.l200.fasta",
		kaijuTable = krakenDir + "/{sample}.kaiju.table.txt",
		blastnrOut0 = annotateDir + "/{sample}.blastnr0.txt",
	output:
		blastnrOut1 = temp(annotateDir + "/{sample}.blastnr1.txt"),
		blastnrOut = annotateDir + "/{sample}.blastnr.txt",  #blast NCBI all nr
		blastnrKrona = annotateDir + "/{sample}.blastnr.krona.html",
		blastnrSummary = annotateDir + "/{sample}.blastnr.summary.txt"
	message:
		'''--- summarize {wildcards.sample} contigs blastnr results'''
	log:
		logDir + "/annotate/{sample}.summarizeBlastnr.log"
	shell:
		"""
		#add kaiju report to the blastnr result
		python {scripts_dir}/mergeBlastxKaiju.py -b {input.blastnrOut0} -k {input.kaijuTable} -c {input.contigs} -o {output.blastnrOut}
		python {scripts_dir}/summarizeBlast.py {output.blastnrOut} {output.blastnrSummary} 
		#classify blast
		awk -F "\\t" 'OFS="\\t" {{ if($11 < 1e-20) print }}' {input.blastnrOut0} | sort -k 11,11g | head -10000 > {output.blastnrOut1} || true
		ktImportBLAST -o {output.blastnrKrona} {output.blastnrOut1} &>> {log} || true 
		"""

rule merge_blastnx:
	"""
	Merge blastn and blastx results for viruses
	"""
	input:
		blastxOut = annotateDir + "/{sample}.blastx.txt",
		blastnOut = annotateDir + "/{sample}.blastn.txt",
	output:
		blastSummary = annotateDir + "/{sample}.blast.nx.txt"
	params:
		contigs = annotateDir + "/{sample}.contigs.l200.fasta",
	message:
		'''--- summarize {wildcards.sample} contigs blastn + blastx results'''
	shell:
		"""
		#merge blastn and blastx results
		python {scripts_dir}/mergeBlastnx.py -n {input.blastnOut} -x {input.blastxOut} -c {params.contigs} -o {output.blastSummary}  
		"""
