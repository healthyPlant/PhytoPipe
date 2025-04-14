#! /usr/bin/env python
######################################################################
#Alex Hu <xiaojun.hu@usda.gov>
#Updated: 07/15/2021
#This program adds Kaiju report to blastx table
#
#Required Parameters:
#   -k, --kaiju X..........................Kaiju output table 
#   -b, --blast X..........................blastx table file
#   -c, --contig X.........................contig file
#   -o, --output X.........................output file 
######################################################################

import os, sys, argparse
import re
from Bio import SeqIO

#################################### Get blast table result function ####################################
#blast format: qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore stitle 
#TRINITY_DN2_c0_g1_i9	gi|1735504366|gb|CP035708.1|	89.108	2332	2317	3983768	13	2317	3973872	3971557	0	2859	Sphaerotilus natans subsp. sulfidivorans strain D-507 chromosome, complete genome
#TRINITY_DN55_c0_g1_i1	RXH91078.1	100	53	461	218	2	160	166	218	1.30E-20	110.2	RXH91078.1 hypothetical protein DVH24_020101 [Malus domestica]
def getBlast(blastFile):
    """
    Get blast result, keep the top one (highest bitscore) for a query
    """
    blastDict = {}
    with open(blastFile) as f:
        for line in f:
            cells = line.rstrip().split('\t')
            if(id not in blastDict or float(blastDict[cells[0]][11]) < float(cells[11])): #only keep top one (the highest bitscore) for a query
                blastDict[cells[0]] = cells
    return blastDict

#################################### Get kaiju table result function ####################################
#file	percent	reads	taxon_id	taxon_name
#OCT19-Q49548.kaiju.out	0.245875	8023	351495	Viruses;NA;NA;NA;Partitiviridae;NA;Raphanus sativus cryptic virus 2;
def getKaiju(kaijuFile):
    """
    Get kaiju result
    """
    kaijuDict = {} #use species name as a key
    #minimum_count = 5  #read cut off for a taxon
    with open(kaijuFile) as f:
        next(f)  #skip header
        for line in f:
            cells = line.rstrip().split('\t')
            #if(float(cells[2]) >= minimum_count and cells[3] != 'NA'): #reads >= 5 and has taxon, keep it 
            #    kaijuDict[cells[3]] = [cells[2], cells[4]] #taxon_id = (reads, name)
            if(cells[3] != 'NA'): #taxon exists
                species = cells[4].split(";")[-2]
                kaijuDict[species] = cells[2:5] # species = (reads,taxon_id, taxon path)
    return kaijuDict


#################################### Merge Kaiju and Blastx by species name ####################################
def addKaiju2Blast(blastDict, kaijuDict, contig_file):
    """
    add kaiju result (taxon id, read#, taxon name) to blastx table by common species name
    """
    #put contig sequences in a dict
    contig_dict = SeqIO.index(contig_file, "fasta")
    #print(kaijuDict)
    #add kaiju result to blast (left merge)
    for qid in blastDict:
        seq = contig_dict[qid]
        #print(seq)
        #print(seq.seq)  
        #print(blastDict[qid][12])
        try:
            species = re.search(r"\[(.*?)\]$", blastDict[qid][12]).group(1)  #get species name from blastx reference title [..]
            #print(species)
        except AttributeError:
            species = 'nan'
        #print(qid, species)
        if len(blastDict[qid]) > 13:
            blastDict[qid] = blastDict[qid][:13]
        blastDict[qid].append(species)
        if species in kaijuDict:
            blastDict[qid].extend(kaijuDict[species])
        else:
            blastDict[qid].extend(['nan','nan', 'nan'])
        #print(contig_dict[qid])
        blastDict[qid].append(str(seq.seq))
        
#################################### Get arguments function ####################################
def parseArguments():
    """
    Get arguments
    """
    parser = argparse.ArgumentParser(description = "Merge Kaiju result to Blastx table")
    parser.add_argument('-b','--blast',dest='blast_file', required=True, default="", help='Blastx result table')
    parser.add_argument('-k','--kaiju',dest='kaiju_file', required=True, default="", help='Kaiju table output from kaiju2table')
    parser.add_argument('-c','--contig',dest='contig_file', required=True, default="", help='Contig file from assembler output, ex. contigs.txt ')
    parser.add_argument('-o','--output',dest='output', required=True, default="", help='Output file' )
    parser.set_defaults(append=False)
    return parser.parse_args()

def main():
    ### Input arguments
    options = parseArguments()
    blast_file = options.blast_file
    kaiju_file = options.kaiju_file
    contig_file = options.contig_file
    output_file = options.output
    
    blastDict = getBlast(blast_file)
    kaijuDict = getKaiju(kaiju_file)
    addKaiju2Blast(blastDict, kaijuDict, contig_file)
    
    #output result
    fout = open(output_file, 'w')
    fout.write("ContigName\tReference\tPercentIdentity\tAlignLength\tContigLength\tRefLen\tContigAlignStart\tContigAlignEnd\tRefAlignStart\tRefAlignEnd\tEvalue\tBitscore\tRefTitle\tRefTaxonId\tKaijuRead\tSpecies\tTaxonPath\tContig\n") 
    #print(blastDict)
    for id in blastDict:
        #print(blastDict[id])
        #fout.write("\t".join(blastDict[id]) + "\n")
        fout.write("\t".join(blastDict[id][0:13]) + "\t")
        fout.write(blastDict[id][15] + "\t" + blastDict[id][14] + "\t" + blastDict[id][13] + "\t" + blastDict[id][16] + "\t" + blastDict[id][17] +"\n")
    fout.close()

if __name__ == '__main__':
	main()
