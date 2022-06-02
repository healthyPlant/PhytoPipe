##Mapping reads to pseudo-contig
#!/usr/bin/bash

ref=$1  #reference file
contigs=$2 #contigs from de novo assembly
readName=$3 #sample name
read1=$4 #read1 (R1) file
read2=$5 #read1 (R2) file

refName=`basename $ref .fasta`
#readName=`basename $reads .fastq.gz`
outputName=$readName"."$refName".pseudoCcontig"
scripts_dir=/ppq/data1/software/ppq_pipeline/scripts

#Align contigs to the reference and generate a reference adjusted pseudo-contig

#run Viral-NGS order and orientied
errormessage=$(python $scripts_dir/orderOrient.py -i $contigs -r $ref -o $outputName.fasta 2>&1)
echo $errormessage
#if [[ -f {output.pseudoContig}/$refName.pseudoContig.fasta ]]; then
#Select longer contigs from Spades original contigs and Viral-NGS refined pseudo-contigs
#python $scripts_dir/selectContigs.py -r {input.refNamex} -b {input.selectedRef} -c {input.contigs} -p {output.pseudoContig} -o {output.finalAssembly}

#mapping reads to finalAssembly
ref1=$outputName.fasta
bwa index $ref1
#bowtie2-build $ref1 $ref1
#bwa mem {params.param} $ref1 {input.reads} | samtools view -Sb -F 4 - > {output.map2Contig}/$refName1.bam  2>> {log}   #filter out unmapped reads -F 4
#samtools sort -o {output.map2Contig}/$refName1.sorted.bam {output.map2Contig}/$refName1.bam

if [[ ! -f $ref1.bwt ]]; then
	bwa index $ref1
fi
if [ -z "$read2" ];then
	bwa mem -k 12 -A 1 -B 3 -O 1 -E 1 $ref1 $read1 > $outputName.sam
	#bowtie2 {params.param} -x $ref1 -U  {input.reads} -S {output.map2Contig}/$refName1.sam 2>> {log}
else
	bwa mem -k 12 -A 1 -B 3 -O 1 -E 1 $ref1 $read1 $read2 > $outputName.sam 
	##bowtie2 {params.param} -x $ref1 -1 {input.reads[0]} -2 {input.reads[0]} -S {output.map2Contig}/$refName1.sam 2>> {log}
fi

samtools view -Sb -F 4 $outputName.sam | samtools sort -o $outputName.sorted.bam -
samtools depth -aa $outputName.sorted.bam > $outputName.coverage.txt

#get mapped read number
echo -n "mapped read number:"
samtools view -c -F 260 $outputName.sorted.bam
count0=$(grep -w 0$ $outputName.coverage.txt | wc -l)
total=$(cat $outputName.coverage.txt | wc -l)
echo -n "Genome covered %:"
echo "($total-$count0)/$total*100" | bc -l
echo -n "Reference lenght: "
wc -l  $outputName.coverage.txt
echo -n "Mean coverage: "
awk '{ total += $3; count++ } END { print total/count }' $outputName.coverage.txt

echo "Plot coverage graph"
python $scripts_dir/plotViralCoverage.py -i $outputName.coverage.txt -o $outputName.coverage.png
rm $outputName.sam  

  #generate consensus sequence 
#normalize the vcf by bcftools norm
bcftools mpileup -Ou -f $ref1 $outputName.sorted.bam | bcftools call -Ou -mv | bcftools norm -f $ref1 -Oz -o $outputName.vcf.gz
bcftools index $outputName.vcf.gz
bcftools consensus -f $ref1 $outputName.vcf.gz -o $outputName.consensus.fasta
