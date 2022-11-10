"""
    Retrieve selected reference sequences from NCBI, then map reads to the reference and plot the coverage graph
    Known viruses are selected by BLASTn and kraken2. Novel viruses are seelcted by BLASTx and Kaiju
"""
__author__ = 'Xiaojun Hu <xiaojun.hu@usda.gov>'

trimDir = config["run_info"]["trim"]
annotateDir = config["run_info"]["annotate"]
krakenDir = config["run_info"]["classify"]
logDir = config["run_info"]["log"]
mapDir = config["run_info"]["map"]
cleanDir = config["run_info"]["clean"]
reportDir = config["run_info"]["report"]
novelDir = config["run_info"]["novel"]
seq_type = config["seq_type"]
nt = config["blastnDb"]  #NCBI blastdb nt
viralRefDb = config["blastnViralDb"]
mappingTool = config["mappingTool"]
bwa_param = config["bwa_param"]
bowtie2_param = config["bowtie2_param"]
mapReadType = config["mapReadType"]

# Define input files
def map_inputs(wildcards):
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

    return reads

rule retrieve_references:
    """
    Retrieve reference sequeces from NCBI
    """
    input:
	    blastnOut = annotateDir + "/{sample}.blastn.txt",
	    blastxOut = annotateDir + "/{sample}.blastx.txt",
	    contigs = annotateDir + "/{sample}.contigs.l200.fasta",
	    kaijuTable = krakenDir + "/{sample}.kaiju.table.txt",
	    krakenReport = krakenDir + "/{sample}.kraken2.report.txt"     
    output:
	    selectedRef = annotateDir + "/{sample}.selectedRef.txt",
	    blastnx = annotateDir + "/{sample}.blast.nx.txt",
	    refNamen = annotateDir + "/{sample}.refNamen.txt",
	    refNamex = annotateDir + "/{sample}.refNamex.txt",
	    refSeq = directory(mapDir + "/ref/{sample}"),
	    done = touch(logDir + "/checkPoint/{sample}.retrieveRef.done")
    message:
	    '''--- retrieve {wildcards.sample} reference sequeces from NCBI.'''
    shell:
      """
      if [ ! -d {output.refSeq} ]; then 
        mkdir {output.refSeq}  
      fi    
      #sort blastn.txt by bitsocre
      sort -nk12,12 -r {input.blastnOut} > {input.blastnOut}.sorted    
      python {scripts_dir}/selectCandidate.py -n {input.blastnOut}.sorted -x {input.blastxOut} -k {input.krakenReport} -j {input.kaijuTable} -c {input.contigs} -s {output.selectedRef} -a {output.blastnx} -l {output.refNamen} -r {output.refNamex}
      if [[ -s {output.refNamen} ]]; then   #check if file is not empty
        bash {scripts_dir}/getSeq.sh {nt} {viralRefDb} {scripts_dir} {output.refNamen} {output.refSeq}
      fi
      if [[ -s {output.refNamex} ]]; then   #check if file is not empty
        bash {scripts_dir}/getSeq.sh {nt} {viralRefDb} {scripts_dir} {output.refNamex} {output.refSeq}
      fi
      """

rule map_ref:
    """
    mapping reads to references, plot coverage graph and make consensus
    """
    input:
	    reads = map_inputs,
	    refSeq = mapDir + "/ref/{sample}"
    output:
	    map2Ref = directory(mapDir + "/map2Ref/{sample}"),
	    done = touch(logDir + '/checkPoint/{sample}.map2Ref.done'), 
    message:
	    '''--- mapping {wildcards.sample} reads to references.'''
    priority: -10
    params:
      mappingTool = mappingTool,
	    seq_type = seq_type
    threads: config["number_of_threads"]
    log:
	    logDir + "/mapping/{sample}.map2Ref.log"
    shell:
	    """
      if [ ! -d {output.map2Ref} ]; then 
        mkdir {output.map2Ref}  
      fi   
      if [ $(ls -l {input.refSeq} | wc -l) -gt 1 ]; then
        for ref in `ls {input.refSeq}/*.fasta` 
        do
          refName=$(basename $ref .fasta)
          if [[ "{params.mappingTool}" == "bowtie2" ]]; then
            bowtie2-build $ref $ref
            #mapping reads to a reference       
            if [[ "{params.seq_type}" == "pe" ]]; then
              bowtie2 {bowtie2_param} --threads {threads} -x $ref -1 {input.reads[0]} -2 {input.reads[0]} -S {output.map2Ref}/$refName.sam 2>> {log}
            elif [[ "{params.seq_type}" == "se" ]]; then
              bowtie2 {bowtie2_param} --threads {threads} -x $ref -U  {input.reads} -S {output.map2Ref}/$refName.sam 2>> {log}
            fi
          else
            bwa index $ref
            #mapping reads to a reference       
            if [[ "{params.seq_type}" == "pe" ]]; then
              bwa mem {bwa_param} -t {threads} $ref {input.reads[0]} {input.reads[1]} > {output.map2Ref}/$refName.sam 2>> {log}
            elif [[ "{params.seq_type}" == "se" ]]; then
              bwa mem {bwa_param} -t {threads} $ref {input.reads} > {output.map2Ref}/$refName.sam 2>> {log}
            fi
          fi

          samtools view -Sb -F 4 {output.map2Ref}/$refName.sam | samtools sort -o {output.map2Ref}/$refName.sorted.bam - #filter out unmapped reads -F 4
          #draw coverage graph
          rm {output.map2Ref}/$refName.sam
          samtools depth -aa {output.map2Ref}/$refName.sorted.bam > {output.map2Ref}/$refName.coverage.txt
          if [[ -s {output.map2Ref}/$refName.coverage.txt ]]; then
            python {scripts_dir}/plotViralCoverage.py -i {output.map2Ref}/$refName.coverage.txt -o {output.map2Ref}.$refName.coverage.png
            #make consensus
            bcftools mpileup -Ou -f $ref {output.map2Ref}/$refName.sorted.bam | bcftools call -M -Ou -mv --ploidy=1 | bcftools norm -c w -f $ref -Oz -o {output.map2Ref}/$refName.vcf.gz
            # -M, --keep-masked-ref           keep sites with masked reference allele (REF=N)
            #-c, --check-ref <e|w|x|s>         check REF alleles and exit (e), warn (w), exclude (x), or set (s) bad sites [e]
            bcftools index {output.map2Ref}/$refName.vcf.gz
            bcftools consensus -f $ref {output.map2Ref}/$refName.vcf.gz -o {output.map2Ref}/$refName.consensus.fasta --exclude 'REF=="N"'
            #--exclude 'REF=="N"' to solve the problem "The fasta sequence does not match the REF allele N"
            echo "map to "$refName
            bedtools genomecov -bga -ibam {output.map2Ref}/$refName.sorted.bam -g {output.map2Ref}/$refName.consensus.fasta > {output.map2Ref}/$refName.bed
            awk "{{if(\$4 == 0) print}}" {output.map2Ref}/$refName.bed > {output.map2Ref}/$refName.0.bed  #note {{ and \ in awk
            if [[ -s {output.map2Ref}/$refName.0.bed ]]; then  #if {output.map2Ref}/$refName.0.bed is not empty
              bcftools consensus -m {output.map2Ref}/$refName.0.bed -f $ref {output.map2Ref}/$refName.vcf.gz -o {output.map2Ref}.$refName.consensus.N.fasta --exclude 'REF=="N"'
            else
              cp {output.map2Ref}/$refName.consensus.fasta {output.map2Ref}.$refName.consensus.N.fasta
            fi
          fi
        done
      fi
      """

rule map_contig:
    """
    mapping reads to pseudo-contigs, plot coverage graph and make consensus. This is only for novel virus
    """
    input:
	    reads = map_inputs,
	    contigs = assembleDir + "/{sample}/contigs.fasta",
	    refSeq = mapDir + "/ref/{sample}",
	    selectedRef = annotateDir + "/{sample}.selectedRef.txt",
	    refNamex = annotateDir + "/{sample}.refNamex.txt",
    output:
      map2Contig = directory(novelDir + "/map2Contig/{sample}"),
      finalAssembly = directory(novelDir + "/finalAssembly/{sample}"),
      pseudoContig = directory(novelDir + "/pseudoContig/{sample}"),
      done = touch(logDir + '/checkPoint/{sample}.map2Contig.done'), 
    message:
	    '''--- mapping {wildcards.sample} reads to pseudo-contigs.'''
    priority: -20
    params:
      mappingTool = mappingTool,
      seq_type = seq_type
    threads: config["number_of_threads"] #8
    log:
	    logDir + "/mapping/{sample}.map2Contig.log"
    shell:
	    """
      if [ ! -d {output.map2Contig} ]; then 
        mkdir -p {output.map2Contig}  
      else
        rm {output.map2Contig}/*
      fi
      if [ ! -d {output.finalAssembly} ]; then 
        mkdir -p {output.finalAssembly}  
      else
        rm -r {output.finalAssembly}/*
      fi 
      if [ ! -d {output.pseudoContig} ]; then 
        mkdir -p {output.pseudoContig}  
      else
        rm -r {output.pseudoContig}/*
      fi       
      #if blastx find references      
      if [[ -s {input.refNamex} ]]; then   #check if file is not empty  
        #process each reference
        lines=`cat {input.refNamex}`
        for refName in $lines; do
          #echo "$refName"
          #Align contigs to the reference and generate a reference adjusted pseudo-contig
          ref={input.refSeq}/$refName.fasta
          #run Viral-NGS order and orientied to make a pseudoContig
          errormessage=$(python {scripts_dir}/orderOrient.py -i {input.contigs} -r $ref -o {output.pseudoContig}/$refName.pseudoContig.fasta 2>&1)
          #echo $errormessage

          #Select longer contigs from original contigs and Viral-NGS refined pseudo-contigs
          python {scripts_dir}/selectContigs.py -r {input.refNamex} -b {input.selectedRef} -c {input.contigs} -p {output.pseudoContig} -o {output.finalAssembly}

          #mapping reads to finalAssembly
          refName1=$refName".contig"
          ref1={output.finalAssembly}/$refName1.fasta
          if [[ "{params.mappingTool}" == "bowtie2" ]]
          then
            bowtie2-build $ref1 $ref1
            if [[ "{params.seq_type}" == "pe" ]]
            then
              bowtie2 {bowtie2_param} --threads {threads} -x $ref1 -1 {input.reads[0]} -2 {input.reads[0]} -S {output.map2Contig}/$refName1.sam 2>> {log}
            elif [[ "{params.seq_type}" == "se" ]]
            then
              bowtie2 {bowtie2_param} --threads {threads} -x $ref1 -U  {input.reads} -S {output.map2Contig}/$refName1.sam 2>> {log}
            fi
          else
            bwa index $ref1
            if [[ "{params.seq_type}" == "pe" ]]
            then
              bwa mem {bwa_param} -t {threads} $ref1 {input.reads[0]} {input.reads[1]} > {output.map2Contig}/$refName1.sam 2>> {log} 
            elif [[ "{params.seq_type}" == "se" ]]
            then
              bwa mem {bwa_param} -t {threads} $ref1 {input.reads} > {output.map2Contig}/$refName1.sam 2>> {log} 
            fi
          fi

          samtools view -Sb -F 4 {output.map2Contig}/$refName1.sam | samtools sort -o {output.map2Contig}/$refName1.sorted.bam -
          samtools depth -aa {output.map2Contig}/$refName1.sorted.bam > {output.map2Contig}/$refName1.coverage.txt
          python {scripts_dir}/plotViralCoverage.py -i {output.map2Contig}/$refName1.coverage.txt -o {output.map2Contig}.$refName1.coverage.png
          rm {output.map2Contig}/$refName1.sam  
        
          #generate consensus sequence 
          #normalize the vcf by bcftools norm
          bcftools mpileup -Ou -f $ref1 {output.map2Contig}/$refName1.sorted.bam | bcftools call -Ou -mv | bcftools norm -f $ref1 -Oz -o {output.map2Contig}/$refName1.vcf.gz
          bcftools index {output.map2Contig}/$refName1.vcf.gz
          bcftools consensus -f $ref1 {output.map2Contig}/$refName1.vcf.gz -o {output.map2Contig}.$refName1.consensus.fasta
         
        done 
      fi
	    """      

