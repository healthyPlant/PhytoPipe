#!/usr/bin/env bash
#This script is for building the custom plant virus blastn nucleotide database
#The virus taxon id is the column one in taxon.txt file
#bash /ppq/data0/software/PhytoPipe/scripts/buildBlastnDb.sh /ppq/data0/software/PhytoPipe/db/plantvirus_taxon.txt /ppq/data0/software/PhytoPipe/db/unwantedSeqId.txt /ppq/data2/resources/ncbi_nt/nt /ppq/data2/resources/ncbi/taxonomy ./plantvirus.fasta ./plantvirus.gb_taxon.txt

taxonFile="$1"  #/PhytoPipe/db/plantvirus_taxon.txt
excludeSeqId="$2" #exclude seqid in this file
nt_db="$3"   #/my/database/ncbi_nt/nt
taxonomy_db="$4"   #/my/database/taxonomy
outSeq="$5" #/my/database/ncbi/plantvirus.fasta
outSeqTaxon="$6" #/my/database/ncbi/plantvirus.gb_taxon.txt

if [[ ! -f $taxonomy_db/nucl_gb.accession2taxid || ! -f $taxonomy_db/nodes.dmp || ! -f $taxonomy_db/names.dmp ]]; then
	echo "plase make sure three files: nucl_gb.accession2taxid, nodes.dmp, names.dmp are in the folder $taxonomy_db"
	exit
fi


#get virus taxon id
cut -f1 $taxonFile | sed '/TaxonId/d'> $taxonomy_db/taxIDs

#extract sequences from NCBI nt database
blastdbcmd -db $nt_db -dbtype nucl -taxidlist $taxonomy_db/taxIDs -out $outSeq.0.fasta
#remove unwanted sequences
filterbyname.sh ow=t include=f in=$outSeq.0.fasta out=$outSeq names=$excludeSeqId

makeblastdb -in $outSeq -dbtype nucl

#get sequence taxonomy information
#get viral reference taxonomy using kaiju-addTaxonNames
grep ">" $outSeq | cut -d " " -f 1 | sed 's/>//g' > $taxonomy_db/seqid
parallel -k --pipepart -a $taxonomy_db/nucl_gb.accession2taxid --block 100M fgrep -F -f $taxonomy_db/seqid | cut -f 2,3 > $taxonomy_db/taxonId.txt
awk -v FS='\t' -v OFS='\t'  '{print "C",$1,$2}'  $taxonomy_db/taxonId.txt > $taxonomy_db/kaiju.out 
kaiju-addTaxonNames -t $taxonomy_db/nodes.dmp -n $taxonomy_db/names.dmp -i $taxonomy_db/kaiju.out -o $taxonomy_db/gb_taxon1.txt -p
cut -f 2- $taxonomy_db/gb_taxon1.txt > $outSeqTaxon

rm -rf $outSeq.0.fasta
rm -rf $taxonomy_db/taxIDs
rm -f $taxonomy_db/seqid
rm -f $taxonomy_db/taxonId.txt
rm -f $taxonomy_db/gb_taxon1.txt
rm -f $taxonomy_db/kaiju.out