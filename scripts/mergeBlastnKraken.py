#! /usr/bin/env python
######################################################################
#Alex Hu <xiaojun.hu@usda.gov>
#Updated: 07/15/2021
#This program adds Kraken2 report to blastn table
#
#Required Parameters:
#   -k, --kraken X.........................Kraken2 report 
#   -b, --blast X..........................blastn table file
#   -a, --accesion X.......................Accession and its taxon id table
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

#################################### Get kraken2 result function ####################################
#  3.79  170939  170939  U       0       unclassified
# 96.21  4344464 41368   R       1       root
# 95.01  4290160 62774   R1      131567    cellular organisms
# 92.26  4165782 192822  D       2759        Eukaryota
def getKraken(kFile):
    """
    read kraken2 report, keep Number of reads assigned directly to the taxon, taxon id, species name
    """
    #minimum_count = 5  #read cut off for a taxon
    kDict = {}
    bFlag = vFlag = fFlag = oFlag = 0
    with open(kFile) as f:
        for line in f:
            line = line.rstrip()
            #cells = line.split("\t")
            cells = list(map(str.strip, line.split("\t")))
            
            if cells[3]=='D' or cells[3]=='K':
                bFlag = fFlag = oFlag = 0

            if cells[4] == '4751': #Fungi
                bFlag = vFlag = oFlag = 0
                fFlag = 1
                continue
            elif cells[4] == '2': #Bacteria
                vFlag = fFlag = oFlag = 0
                bFlag = 1
                continue
            elif cells[4] == '10239': #Viruses
                bFlag = fFlag = oFlag = 0
                vFlag = 1
                continue
            elif cells[4] == '4762': #Oomycota
                bFlag = vFlag = fFlag = 0
                oFlag = 1
                continue

            #if float(cells[2]) >= minimum_count and cells[3] == 'S':
            if cells[3].startswith('S'):  #species
                if fFlag ==1 :  #species
                    kDict[cells[4]] = [cells[2],cells[5], 'Fungi'] #taxon id = (read#, species name)
                    continue
                elif bFlag ==1 :
                    kDict[cells[4]] = [cells[2],cells[5], 'Bacteria'] #taxon id = (read#, species name)
                    continue
                elif vFlag ==1:  #species
                    kDict[cells[4]] = [cells[2],cells[5], 'Viruses'] #taxon id = (read#, species name)
                    continue
                elif oFlag ==1:  #species
                    kDict[cells[4]] = [cells[2],cells[5], 'Oomycetes'] #taxon id = (read#, species name)
                    continue

    return kDict

#################################### Get accession and its taxon id ####################################
#gb|ABM53544.1| 417305
#emb|CBA31922.1| 667019
#ref|XM_030624382.1| 3483
#dbj|BAG46934.1| 395019
def getAccTi(accTiFile):
    """
    get taxon id for accession, the accFile from " blastdbcmd -db nt -entry_batch acc_list.txt -outfmt "%i %T" > acc_list_ti.txt" 
    """
    accTiDict = {}
    with open(accTiFile) as f:
        for line in f:
            if line.startswith("Error: "):
                continue
            cells = line.rstrip().split()
            accTiDict[cells[0]] = cells[1]

    return accTiDict

#################################### Merge Kaiju and Blast using the relation in acc2ti ####################################
def addKraken2Blast(blastDict, krakenDict, accTiDict, contig_file):
    """
    add kaiju result (taxon id, read#, taxon name) to blastn table
    """
    #put contig sequences in a dict
    contig_dict = SeqIO.index(contig_file, "fasta")

    #add taxon id to blast table
    for qid in blastDict:
        if len(blastDict[qid]) > 13:
            blastDict[qid] = blastDict[qid][:13]
            
        sid = re.sub("gi\|\d+\|", "", blastDict[qid][1])  #remove gi|1731737253|
        if sid in accTiDict:
            ti = accTiDict[sid]
            blastDict[qid].append(ti)
        else:
            blastDict[qid].append('nan')

    #add kraken result to blast
    for qid in blastDict:
        ti = blastDict[qid][-1]
        if ti in krakenDict:
            blastDict[qid].extend(krakenDict[ti])
        else:
            blastDict[qid].extend(['nan','nan', 'nan'])
        blastDict[qid].append(str(contig_dict[qid].seq))


#################################### Get arguments function ####################################
def parseArguments():
    """
    Get arguments
    """
    parser = argparse.ArgumentParser(description = "Merge Kraken result to Blastn table")
    parser.add_argument('-b','--blast',dest='blast_file', required=True, default="", help='Blast result table')
    parser.add_argument('-k','--kraken',dest='kraken2_file', required=True, default="", help='Kraken2 report')
    parser.add_argument('-a','--accesion',dest='acc2ti_file', required=True, default="", help='Accession and its taxon id table')
    parser.add_argument('-c','--contig',dest='contig_file', required=True, default="", help='Contig file from assembler output, ex. contigs.txt ')
    parser.add_argument('-o','--output',dest='output', required=True, default="", help='Output file' )
    parser.set_defaults(append=False)
    return parser.parse_args()

def main():
    ### Input arguments
    options = parseArguments()
    blast_file = options.blast_file
    kraken_file = options.kraken2_file
    acc2ti_file = options.acc2ti_file
    contig_file = options.contig_file
    output_file = options.output
    
    blastDict = getBlast(blast_file)
    krakenDict = getKraken(kraken_file)
    accTiDict = getAccTi(acc2ti_file)
    #print(krakenDict)

    addKraken2Blast(blastDict, krakenDict, accTiDict, contig_file)
    
    #output selected blast result
    fout = open(output_file, 'w')
    fout.write("ContigName\tReference\tPercentIdentity\tAlignLength\tContigLength\tRefLen\tContigAlignStart\tContigAlignEnd\tRefAlignStart\tRefAlignEnd\tEvalue\tBitscore\tRefTitle\tRefTaxonId\tKrakenRead\tKrakenSpecies\tKrakenMicrobialGroup\tContig\n") 
    for id in blastDict:
        #print(blastDict[id])
        fout.write("\t".join(blastDict[id]) + "\n")
    fout.close()

if __name__ == '__main__':
	main()
