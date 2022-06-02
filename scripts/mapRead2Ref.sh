#!/usr/bin/env bash
######################################################################
#Alex Hu <xiaojun.hu@usda.gov>
#Updated: 02/18/2021
#This program maps reads to the reference and draw coverage graph
#Usage: mapRead2Ref.sh [script directory] [viral reference fasta file]  [sample name] [read1] [read2]
######################################################################

scripts_dir=$1
ref=$2  #reference file
readName=$3 #sample name
read1=$4 #read1 (R1) file
read2=$5 #read1 (R2) file

#Please change mapping tools and their peremeters here
mappingTool="bwa" # or "bowtie2"
bwaOption=" -B 4 " #default bwa options. it's good for read length >=150, change mismatch panelty from 1 to 3 (-B)
#bwaOption=" -k 12 -A 1 -B 3 -O 1 -E 1 " # good for read length <=100, Viral-NGS use it.
#bwaOption=" -B 9 -O 16 " #BWA strict mapping

#uncomment the following lines to use bowtie2
#mappingTool="bowtie2"
#Bowtie2 parameters
#bowtieOption=" --mp 20  --score-min L,-0.1,-0.1 " #read almost exactly match to the reference. It's good for mapping very similar virus genomes
#bowtieOption=" --very-sensitive-local -k 100 --score-min L,20,1.0 " #mapping reads with mismatches. Pathoscope2 used it. 

refName=`basename $ref .fasta`
#readName=`basename $reads .fastq.gz`
outputName=$readName"."$refName

if [[ $mappingTool == "bowtie2" ]]; then
	bowtie2-build $ref $ref

	if [ -z "$read2" ];then
		bowtie2 $bowtieOption -x $ref -U  $read1 | samtools view -Sb -F 4 - | samtools sort - > $outputName.sorted.bam
	else
		bowtie2 $bowtieOption -x $ref -1 $read1 -2 $read2 | samtools view -Sb -F 4 - | samtools sort - > $outputName.sorted.bam
	fi
else
	bwa index $ref
	#bwa mem -t {threads} {params.param} $ref {input.reads} | samtools sort -@{threads} -o {params.outputName}/$refName.bam - 2>> {log} 1>&2
	if [ -z "$read2" ];then
		bwa mem $bwaOption  $ref $read1 | samtools view -Sb -F 4 - | samtools sort - > $outputName.sorted.bam  #filter out unmapped reads -F 4
	else
		bwa mem $bwaOption $ref $read1 $read2 | samtools view -Sb -F 4 - | samtools sort - > $outputName.sorted.bam  #filter out unmapped reads -F 4
	fi
fi

#samtools sort -o $outputName.sorted.bam $outputName.bam
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

echo "Make consensus"
bcftools mpileup -Ou -f $ref $outputName.sorted.bam | bcftools call -Ou -mv | bcftools norm -f $ref -Oz -o $outputName.vcf.gz
bcftools index $outputName.vcf.gz
bcftools consensus -f $ref $outputName.vcf.gz -o $outputName.consensus.fasta
bedtools genomecov -bga -ibam $outputName.sorted.bam -g $outputName.consensus.fasta > $outputName.bed
awk '{if($4 == 0) print}' $outputName.bed > $outputName.0.bed  #note {{ and \ in awk
if [[ -s $outputName.0.bed  ]]; then
  bcftools consensus -m $outputName.0.bed -f $ref $outputName.vcf.gz -o $outputName.consensus.N.fasta
fi