#!/usr/bin/env bash
######################################################################
#Alex Hu <xiaojun.hu@usda.gov>
#Updated: 02/04/2025
#This program run a sample contigs blast against NCBI nt
#Usage: runBlastnt.sh [sample name] [work directory] [NCBI nt directory] [script directory] [microbial Taxon Ids]
######################################################################

sample="$1"
work_dir="$2"
nt_db="$3" #"/ppq/data2/resources/ncbi_nt/nt"
scripts_dir="$4" #"/ppq/data0/software/PhytoPipe/scripts"
microbialTaxonIds="$5" #"/ppq/data2/resources/ncbi/microbial.tids"

blastEvalue=1e-40
contigs="$work_dir/annotation/${sample}.contigs.l200.fasta"
blastntOut0="$work_dir/annotation/${sample}.blastnt0.txt"
krakenReport="$work_dir/classification/${sample}.kraken2.report.txt"
blastntOut="$work_dir/annotation/${sample}.blastnt.txt",  
blastntKrona="$work_dir/annotation/${sample}.blastnt.krona.html"
blastntSummary="$work_dir/annotation/${sample}.blastnt.summary.txt"
acc_list="$work_dir/annotation/${sample}.blastnt.acc_list.txt"
acc_ti="$work_dir/annotation/${sample}.blastnt.acc_ti.txt"
blastntOut1="$work_dir/annotation/${sample}.blastnt1.txt"

# Check if any parameters were provided
if [ $# -eq 0 ]; then
  echo "Usage: $0 [sample name] [work directory] [NCBI nt directory] [script directory] [microbial Taxon Ids]"
  echo "Missing some parameters! Exit"
  exit 1
fi


#https://www.ncbi.nlm.nih.gov/books/NBK569839/
#to use taxid parameter, blastn version must >2.15.0
#taxon id: 2 is bacteria, 10239 is viruses, 4751 is fungi, 4762 is oomycetes
echo "${sample} contigs blastn NCBI nt"
/ppq/data0/software/ncbi-blast-2.13.0+/bin/blastn -db $nt_db -num_threads 8 -query $contigs -out $blastntOut0 -taxidlist $microbialTaxonIds -max_target_seqs 1 -evalue  $blastEvalue -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore stitle "  
#blastn -db $nt_db -num_threads 8 -query $contigs -out $blastntOut0 -taxids 2,10239,4751,4762 -max_target_seqs 1 -evalue  $blastEvalue -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore stitle "  

echo "Merge ${sample} blastnt and Kraken2 results "
#get taxon id (ti) for subjects in blastnt result
cut -f2 $blastntOut0 | uniq > $acc_list

blastdbcmd -db $nt_db -entry_batch $acc_list -outfmt "%i %T"  > $acc_ti 
#add kraken2 report to the blastnt result
python $scripts_dir/mergeBlastnKraken.py -b $blastntOut0 -k $krakenReport -a $acc_ti -c $contigs -o $blastntOut
python $scripts_dir/summarizeBlast.py $blastntOut $blastntSummary 
#classify blast
awk -F "\t" 'OFS="\t" { if($11 < 1e-20) print }' $blastntOut0 | sort -k 11,11g  | head -10000 > $blastntOut1
ktImportBLAST -o $blastntKrona $blastntOut1
 
rm -f $acc_list
rm -f $acc_ti
rm -f $blastntOut1

echo "${sample} blastnt is done. Please check the result in ${blastntOut} and ${blastntSummary}."
