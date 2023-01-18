"""
    Make a final report including table report and html report
"""
__author__ = 'Xiaojun Hu <xiaojun.hu@usda.gov>'

workDir = config["workDir"]
mapDir = config["run_info"]["map"]
reportDir = config["run_info"]["report"]
novelDir = config["run_info"]["novel"]
acronymDb = config["acronymDb"]
reportFile = reportDir + "/report.txt"
qcDir = config["run_info"]["qc"]
annotateDir = config["run_info"]["annotate"]
logDir = config["run_info"]["log"]
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
blastEvalue = config["blastEvalue"]

rule blastn_local:
    """
    Blastn againt local NCBI nt
    """
    input:
	    consensus = mapDir + "/map2Ref/{sample}",
	    contig = novelDir + "/finalAssembly/{sample}"
    output:
	    done = touch(logDir + "/checkPoint/{sample}.blastn.done"),
	    consensusFile = temp(reportDir + "/{sample}.consensus.fa"),  
	    blastnResult = temp(reportDir + "/{sample}.consensus.blastn.txt"),
    message:
	    """--- {wildcards.sample} local Blast consensus/contigs against NCBI nt/nr."""
    priority: -10
    threads: 8
    log:
	    logDir + "/annotate/{sample}.finalConsensus.blastn.log"
    params:
	    blastn = " -word_size 11 -max_target_seqs 1 -evalue 1e-100 -outfmt \"6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore stitle qcovs \" ", #qcovs:Query Coverage Per Subject
	    blastx = " --max-target-seqs 1 -e  " + blastEvalue + "  --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore stitle scovhsp ", #scovhsp:Subject Coverage Per HSP*
	    blastnFailed = reportDir + "/{sample}.consensus.blastn.failed.fasta",
	    blastxResult = reportDir + "/{sample}.consensus.blastx.txt",
	    out = reportDir + "/ncbiBlast",
	    smpName = "{sample}"  #Note: {sample} can't be directly used in shell 
    shell:
        """
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
                #change sequence title
                sed "1s/^.*$/>{params.smpName}.$refName.consensus/" "$consensus" >> {output.consensusFile}  #give the sequence new title (sampleName.reference name)
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
                #change sequence title
                sed "1s/^.*$/>{params.smpName}.$refName.contig/" "$contig" >> {output.consensusFile}  #give the sequence new title (sampleName.reference name)
            done
        fi

        #run blastn against local NCBI nt database
        if [ -f {output.consensusFile} ] && [ -s  {output.consensusFile} ]; then
            blastn -db {nt} -num_threads {threads} -query {output.consensusFile} -out {output.blastnResult} {params.blastn} &>> {log};
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
	    blastDone = expand(logDir + "/checkPoint/{smp}.blastn.done", smp=samples),
	    rawQC = qcDir + "/multiqc/raw_multiqc.html",
	    trimmedQC = qcDir + "/multiqc/trimmed_multiqc.html",
	    refSeqDone = expand(logDir + "/checkPoint/{smp}.retrieveRef.done", smp=samples),
	    map2RefDone = expand(logDir + "/checkPoint/{smp}.map2Ref.done", smp=samples),
	    map2ContigDone = expand(logDir + "/checkPoint/{smp}.map2Contig.done", smp=samples),
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
        """
        
rule generate_htmlReport:
    """
    Generate the final html report
    """
    input:
	    rpt = rules.generate_report.output.rpt   #reportFile,
    output:
	    qcReadNumber = reportDir + "/qcReadNumber.txt",
	    html = htmlReport
    message:
	    '''--- Generate the final html report.'''
    priority: -200
    shell:
	    """
	    #get QC read number
        if [ ! {strand1} ]; then
            python {scripts_dir}/getReadNumber.py -w {workDir} -m {fastq_format} -o {output.qcReadNumber}
        else
            python {scripts_dir}/getReadNumber.py -w {workDir} -f {strand1} -r {strand2} -m {fastq_format} -o {output.qcReadNumber}
        fi
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
        fi

	    #generate a html report
        python {scripts_dir}/getHtmlReport.py {workDir} {monitorPathogen} {workDir}/{reportDir} {input.rpt} {output.html}
		rm -rf {reportDir}/*.failed.fasta
	    """
