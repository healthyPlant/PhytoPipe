#!/usr/bin/env bash
######################################################################
#Alex Hu <xiaojun.hu@usda.gov>
#Updated: 02/18/2021
#This program extracts sequences from multiple sources
#Usage: getSeq.sh [NCBI nt directory] [viral reference fasta file] [script directory] [sequence gb id file] [output directory]
######################################################################

nt_db="$1"  #/ppq/data2/resources/ncbi_nt/nt
viral_db="$2"  #/ppq/data2/resources/ncbi/refseq_viral_genomic.fa
scripts_dir="$3"  #/ppq/data1/software/ppq_pipeline/scripts
seqIdFile="$4"
outputFolder="$5"

#extract sequences from viral reference database
seqtk subseq $viral_db $seqIdFile > $outputFolder/temp.fasta
cat $outputFolder/temp.fasta | awk -v folder=$outputFolder '{ if (substr($0,1,1)==">") {split(substr($0,2),f, " ")} {filename=f[1] ".fasta"} print $0 > folder"/"filename}'
rm -rf $outputFolder/temp.fasta

while IFS= read -r line
do
    echo "$line" | tee $outputFolder/temp.txt
    id=$(echo $line | tr -d '\n')
    #extract sequences from viral reference database
    #seqtk subseq $viral_db temp.txt > $outputFolder/$id.fasta
    
    #extract sequences from local blast nt database
    if [[ ! -f "$outputFolder/$id.fasta" || ! -s "$outputFolder/$id.fasta" ]]; then
        echo "I'm extracting $id from local NCBI blastn database"
        blastdbcmd -db $nt_db -entry $id > $outputFolder/$id.fasta
    fi

    #extract sequences from NCBI website
    if [[ ! -f "$outputFolder/$id.fasta" || ! -s "$outputFolder/$id.fasta" ]]; then
        echo "I'm extracting $id from NCBI website"
        python $scripts_dir/getSeqFromNCBI.py -n $outputFolder/temp.txt -o $outputFolder
    fi   
done < "$seqIdFile"
rm -rf $outputFolder/temp.txt

#fi
