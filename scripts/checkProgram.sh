#!/usr/bin/env bash
# List of programs used in PhytoPipe
programs=("fastqc" "java" "multiqc" "kraken2" "ktImportTaxonomy" "seqtk" "kaiju" "kaiju2table" "kaiju2krona" "ktImportText" "quast.py" "blastn" "ktImportBLAST" "diamond" "blastdbcmd" "samtools" "bcftools" "nucmer" "parallel" "blastdbcmd" "kaiju-addTaxonNames")
nanopore=("flye" "minimap2" "miniasm" "NanoFilt" "seqkit")
illumina=("bbduk.sh" "clumpify.sh" "bbsplit.sh" "Trinity" "spades.py" "bowtie2-build" "bowtie2" "bwa")
# check programs
echo "########################################"
echo "Check the programs used in PhytoPipe"
echo "########################################"
for prog in ${programs[@]};
do
	if ! command -v $prog &> /dev/null
	then
		echo "$prog could not be found. Please intall it and/or put it in the path. "
	else
		echo "$prog is ok"
	fi
done

echo "########################################"
echo "If you analyze Illumina sequencing data, the following additional programs are needed."
echo "########################################"
for prog in ${illumina[@]};
do
	if ! command -v $prog &> /dev/null
	then
		echo "$prog could not be found. Please intall it and/or put it in the ~/bin path. "
	else
		echo "$prog is ok"
	fi
done
