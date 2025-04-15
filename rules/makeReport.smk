"""
    Make a final report including table report and html report
"""
__author__ = 'Xiaojun Hu <xiaojun.hu@usda.gov>'

workDir = config["workDir"]
mapDir = config["workDir"] + "/" + config["run_info"]["map"]
reportDir = config["workDir"] + "/" + config["run_info"]["report"]
novelDir = config["workDir"] + "/" + config["run_info"]["novel"]
qcDir = config["workDir"] + "/" + config["run_info"]["qc"]
annotateDir = config["workDir"] + "/" + config["run_info"]["annotate"]
logDir = config["workDir"] + "/" + config["run_info"]["log"]
acronymDb = config["acronymDb"]
reportFile = reportDir + "/report.txt"
seq_type = config["seq_type"]
strand1 = config["strand1"]
strand2 = config["strand2"]
fastq_format = config["input_format"]
simReportFile = reportFile.replace(".txt", ".simple.txt")
htmlReport = reportFile.replace(".txt", ".html")
nt = config["blastnDb"]  #NCBI blastdb nt
nr = config["blastxDb"]  #diamond format NCBI nr
samples = config["samples"]    
monitorPathogen = config["monitorPathogen"]
excludeTaxids = config["excludeTaxids"]

filterKeyFile = ""
try:
    filterKeyFile = config["filterKeys"]
except KeyError:
    filterKeyFile = ""

samples = config["samples"]    

blastnEvalue = "1e-40" 
blastxEvalue = "1e-20"

def aggregate_input(wildcards):
	"""
	Blast rule depends on the contigs from assembly, if no contigs, skip blast; refNamex is from blastx; refNamen is from combine blastn and mapping
	'aggregate_input' function return values are the rule aggregate input, which are another two rules outputs
	"""
	fn = checkpoints.extract_contig.get(sample=wildcards.sample).output[0]

	if os.stat(fn).st_size > 0: #if having contigs, run blast and mapping
		chosenfile = reportDir + "/{sample}.consensus.blastn.txt"  #rule aggregate_blast output
	else:
		chosenfile = logDir + "/assemble/{sample}.assembly.failed"  #rule check_contig output

	return chosenfile

rule check_contig:
	"""
	Alternative rule for checkpoint extract_contig. Mark assembly failed samples 
	"""
	input:
		assembleDir + "/{sample}/contigs.fasta"
	output:
		logDir + "/assemble/{sample}.assembly.failed"
	message:
		'''--- Check contig file size from assembly.'''
	params:
		annotateDir + "/{sample}.selectedRef.txt",
	shell:
		"""
		if [ ! -s {input} ]; then   #check if the contig file is not empty
			touch {output}
			touch {params}
		fi
		"""


rule aggregate_blast:
	"""
	Aggregate consensus/contig blast result for checkpoint extract_contig
	"""
	input:
		aggregate_input
	output:
		logDir + "/checkPoint/{sample}.blastn.done"
	message:
		'''--- Aggregate blastn/x results for final consensus/contigs.'''
	shell:
		"""
		touch {output}
		#get novel virus mapping information
		if [ -d {novelDir}/map2Contig ]; then
			found=$(find {novelDir}/map2Contig -name '*.fasta' | wc -l)
			if (( $found > 0 )); then
				python {scripts_dir}/getContigMapInfo.py -w {workDir} -o {reportDir}"/novelVirusMapping.txt"
			fi
		fi
		"""

rule blastn_local:
	"""
	Blastn againt local NCBI nt
	"""
	input:
		consensus = mapDir + "/map2Ref/{sample}",
		contig = novelDir + "/finalAssembly/{sample}"
	output:
		#done = touch(logDir + "/checkPoint/{sample}.blastn.done"),
		consensusFile = temp(reportDir + "/{sample}.consensus.fa"),  
		blastnResult = temp(reportDir + "/{sample}.consensus.blastn.txt"), 
	message:
		"""--- {wildcards.sample} local Blast consensus/contigs against NCBI nt/nr."""
	priority: -10
	threads: 16
	log:
		logDir + "/annotate/{sample}.finalConsensus.blastn.log"
	params:
		blastn = " -word_size 11 -max_target_seqs 1 -evalue " + blastnEvalue + " -outfmt \"6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore stitle qcovs \" ", #qcovs:Query Coverage Per Subject
		blastx = " --max-target-seqs 1 -e  " + blastxEvalue + "  --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore stitle scovhsp ", #scovhsp:Subject Coverage Per HSP*
		blastnFailed = reportDir + "/{sample}.consensus.blastn.failed.fasta",
		blastxResult = reportDir + "/{sample}.consensus.blastx.txt",
		out = reportDir + "/ncbiBlast",
		smpName = "{sample}"  #Note: {sample} can't be directly used in shell 
	shell:
		"""
		#set -x  #for shell debug
		#create empty files in case no virus found in sample
		touch {output.consensusFile}
		touch {output.blastnResult}
		
		#combine all final consensus/contig to a file
		if [ ! -d  {reportDir}/ncbiBlast ]; then
			mkdir -p {reportDir}/ncbiBlast
		fi

		#for consensus
		inputPath={mapDir}/map2Ref 
		fileNum=$(find $inputPath -name "{params.smpName}.*.consensus.N.fasta" | wc -l)
		if (( $fileNum > 0 )); then
			for consensus in `ls $inputPath/{params.smpName}.*.consensus.N.fasta`  #file name has sample name and refernce name
			do
				consensus0=$(basename $consensus '.consensus.N.fasta' )
				refName="$(cut -d'.' -f2- <<<"$consensus0")"
				
				#in a re-run, the old reference is not in the selected file (.selectedRef.txt), check whehter the reference ($refName) is available, if yes, run compareSeq.py
				#grep -q argument does not output the matched text, but only return the exit status code
				if grep -q $refName {annotateDir}/{params.smpName}.selectedRef.txt; then
					#compare a contig and a consensus sequence, if a contig is longer than 90% of a consensus (without N), use it
					python {scripts_dir}/compareSeq.py -r $refName -s {params.smpName} -b {annotateDir}/{params.smpName}.selectedRef.txt -g {assembleDir}/{params.smpName}/contigs.fasta -c $inputPath/{params.smpName}.$refName.consensus.N.fasta -o $inputPath/{params.smpName}.$refName.contig.fasta 
					if [ -f $inputPath/{params.smpName}.$refName.contig.fasta ]; then
						cat $inputPath/{params.smpName}.$refName.contig.fasta >> {output.consensusFile}
					else
						#change sequence title
						sed "1s/^.*$/>{params.smpName}.$refName.consensus/" "$consensus" >> {output.consensusFile}  #give the sequence new title (sampleName.reference name)
					fi
				fi
			done
		fi

		#for contigs
		inputPath={novelDir}/finalAssembly/{params.smpName}
		fileNum=$(find $inputPath -name '*.contig.fasta' | wc -l)
		#echo $fileNum
		if (( $fileNum > 0 )); then
			for contig in `ls $inputPath/*.contig.fasta` 
			do
				refName=$(basename $contig '.contig.fasta')  #remove '.contig.fasta'
				#if the contig is in the file, replace it using the contig from the final assembly
				repName=$(cat {output.consensusFile} | grep $refName | cut -d' ' -f1 | sed 's/>//') || true
				echo $repName
				if [ ! -z $repName ]; then 
					filterbyname.sh in={output.consensusFile} out=temp.fa names=$repName include=f  #remove the contig by its name
					mv temp.fa {output.consensusFile}
					sed "1s/^.*$/>{params.smpName}.$refName.contig/" "$contig" >> {output.consensusFile}  #give the sequence new title (sampleName.reference name)
					rm -rf temp.fa
				fi
			done
		fi

		#run blastn against local NCBI nt database
		if [ -f {output.consensusFile} ] && [ -s  {output.consensusFile} ]; then
			#make soft link for the file taxonomy4blast.sqlite3, which is required for blastn with the parameter -taxids 
			from=$(dirname {nt})
			to=`pwd`
			if [[ -e $from/taxonomy4blast.sqlite3 && ! -e $to/taxonomy4blast.sqlite3 ]]; then
				ln -s $from/taxonomy4blast.sqlite3 $to/
			fi
			if [ -e $to/taxonomy4blast.sqlite3 ]; then
				#run parallel blastn with excluded taxons, each parallel 4 sequences using 4 threads
				cat {output.consensusFile} | parallel -q -j {threads} -N 4 --recstart '>' --pipe blastn -task megablast -db {nt} -num_threads 4 -negative_taxids {excludeTaxids} {params.blastn} > {output.blastnResult}  2>> {log} 
				#rm -f $to/taxonomy4blast.sqlite3
			else
				blastn -task megablast -db {nt} -num_threads {threads} -query {output.consensusFile} -out {output.blastnResult} {params.blastn} &>> {log}
			fi

			#split blastn result to different files, like remote blastn results
			python {scripts_dir}/splitBlastResult.py -s {output.consensusFile} -b {output.blastnResult} -f {params.blastnFailed} -t blastn -o {params.out} 
			#if blastn failed, try blastx
			if [ -f {params.blastnFailed} ] && [ -s  {params.blastnFailed} ]; then
				diamond blastx -d {nr} -q {params.blastnFailed} -p {threads} -o {params.blastxResult} {params.blastx} &>> {log}
				python {scripts_dir}/splitBlastResult.py -s {params.blastnFailed} -b {params.blastxResult} -t blastx -o {params.out}
				rm -f {params.blastnFailed}
				rm -f {params.blastxResult}
			fi
		fi                       
		
		"""

rule generate_report:
	"""
	Generate the final report
	"""
	input:
		#rawQC = qcDir + "/multiqc/raw_multiqc.html",
		#trimmedQC = qcDir + "/multiqc/trimmed_multiqc.html",
		kaijuDone = expand(classifyDir + "/{sample}.kaiju_krona.html", sample=samples),
		krakenDone = expand(classifyDir + "/{sample}.kraken2.report.html", sample=samples),
		blastDone = expand(logDir + "/checkPoint/{sample}.blastn.done",sample=samples),
		#refSeqDone = expand(logDir + "/checkPoint/{sample}.retrieveRef.done", sample=samples),
		#map2RefDone = expand(logDir + "/checkPoint/{sample}.map2Ref.done", sample=samples),
		#map2ContigDone = expand(logDir + "/checkPoint/{sample}.map2Contig.done", sample=samples),
	output:
		rpt = reportFile,
		log = logDir + "/report/report.log"
	message:
		'''--- Generate the final report.'''
	priority: -100
	shell:
		"""
		if [ {seq_type}=='se' ]; then
			python {scripts_dir}/getReport.py -a {acronymDb} -o {output.rpt} -w {workDir} > {output.log}
		else
			python {scripts_dir}/getReport.py -f {strand1} -r {strand2} -a {acronymDb} -o {output.rpt} -w {workDir} > {output.log}
		fi

		rm -rf rm {reportDir}/*.failed.fasta

		if [ -e `pwd`/taxonomy4blast.sqlite3 ]; then
			rm -rf `pwd`/taxonomy4blast.sqlite3
		fi

		"""

rule sum_readStats:
	"""
	Summary raw read stats 
	"""
	input:
		contigs = expand(assembleDir + "/{sample}/contigs.fasta", sample=samples),
		raw_multiqc = qcDir + "/multiqc/raw_multiqc.html",
		trimmed_multiqc = qcDir + "/multiqc/trimmed_multiqc.html",
		quast_multiqc = qcDir + "/multiqc/quast_multiqc.html",
	output:
		reportDir + "/qcReadNumber.txt"
	message:
		'''--- Summary raw read stats.'''
	shell:
		"""
		#get QC read number
		if [ ! {strand1} ]; then
			python {scripts_dir}/getReadNumber.py -w {workDir} -f "" -m {fastq_format} -o {output}
		else
			python {scripts_dir}/getReadNumber.py -w {workDir} -f {strand1} -r {strand2} -m {fastq_format} -o {output}
		fi
		"""


rule generate_htmlReport:
	"""
	Generate the final html report
	"""
	input:
		rpt = rules.generate_report.output.rpt,   #reportFile (report.txt),
		qcrn = reportDir + "/qcReadNumber.txt",
	output:
		log = logDir + "/report/report.log",
		html = htmlReport
	message:
		'''--- Generate the final html report.'''
	params:
		sampleInfo = keyExcelFile,
		filterKey = filterKeyFile,
		filteredRpt = reportDir + "/report.filtered.txt",
		pgqpReport = reportDir + "/report.pgqp.txt",
		simpleReport = reportDir + "/report.simple.txt"
	priority: -200
	shell:
		"""
		#get novel virus mapping information
		found=$(find {novelDir}/map2Contig -name '*.fasta' | wc -l)
		if (( $found > 0 )); then
			python {scripts_dir}/getContigMapInfo.py -w {workDir} -o {reportDir}"/novelVirusMapping.txt"
		fi

		#copy blast result to report folder
		if [ ! -d {reportDir}/blastnx ]; then
			mkdir {reportDir}/blastnx
		fi
		cp {workDir}/annotation/*.blast.nx.txt {reportDir}/blastnx/

		found=$(find {workDir}/annotation -name '*.blastnt.txt' | wc -l)
		if [[ $found >0 ]]; then
			cp {workDir}/annotation/*.blastnt.txt {reportDir}/blastnx/
			cp {workDir}/annotation/*.blastnr.txt {reportDir}/blastnx/
			cp {workDir}/annotation/*.blastnt.summary.txt {reportDir}/blastnx/
			cp {workDir}/annotation/*.blastnr.summary.txt {reportDir}/blastnx/
		fi

		#Filter report using filter keys
		if [ ! -z {params.filterKey} ] && [ -f  {params.filterKey} ]; then
			python {scripts_dir}/filterReport.py {input.rpt} {params.filterKey} {params.filteredRpt} >> {output.log}
		fi

		#generate a html report
		python {scripts_dir}/getHtmlReport.py {workDir} {monitorPathogen} {reportDir} {input.rpt} {output.html}

		"""
