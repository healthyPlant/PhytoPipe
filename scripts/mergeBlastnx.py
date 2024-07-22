#merge blastn and blastx results and contig sequences
#! /usr/bin/env python

import os, sys, argparse
from Bio import SeqIO

#################################### Get blast table result function ####################################
#blastn format: qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore stitle taxonid taxonpath
#TRINITY_DN31_c0_g1_i1	gi|746718929|gb|KM025045.1|	74.281	7652	9923	9624	984	8560	750	8327	0.0	3090	Sorghum mosaic virus isolate YN-YT168, complete genome	32619	Sorghum mosaic virus	Sorghum mosaic virus	Viruses
def getBlastn(blastnFile):
    """
    Get blastn result, if the records have the same query id [0] and taxonomy id [13], keep the higher bit score one
    """
    blastnDict = {}
    idDict = {} #keep unique id
    with open(blastnFile) as f:
        for line in f:
            cells = line.rstrip().split('\t')
            #gb = cells[1]
            #if("gb|" in cells[1]):
                #get gi and gb
            #    gids = cells[1].split('|') #gi|746718929|gb|KM025045.1|
            #    gi = gids[1]  #746718929
            #    gb = gid[3]
                #ngb = gids[3].split('.')[0]  #get genbank id, KM025045
            if len(cells) < 14: #for missing taxonomy
                cells.extend(["","0"])
            id = cells[0]+"|"+cells[13] #query id + tax id 

            if((id not in idDict) or float(idDict[id]) < float(cells[11])): #keep the higher bitscore one
                idDict[id] = cells[11]
                blastnDict[cells[0]] = cells

    return blastnDict

#################################### Get blast table result function ####################################
#blast format: qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore stitle taxonid taxonpath
#NODE_33_length_4928_cov_612.276850	NC_003491.1	78.378	259	4928	5722	1013	1268	4030	3775	8.50e-38	163	NC_003491.1 Beet mild yellowing virus, complete genome	156690	Viruses; Riboviria; Luteoviridae; Polerovirus; Beet mild yellowing virus;
#NODE_289_length_2910_cov_8.945994	acc|GENBANK|QBK89270.1|GENBANK|MK500411|metallopeptidase	28.7	352	2910	369	205	1134	4	347	5.2e-26	124.4	acc|GENBANK|QBK89270.1|GENBANK|MK500411|metallopeptidase family M24 [Mimivirus LCMiAC02]	2506609	Viruses; Mimiviridae; unclassified Mimiviridae; Mimivirus LCMiAC02;
def getBlastx(blastFile):
    """
    Get blastx result, keep the top one (highest bitscore) for a query
    """
    blastDict = {}
    bitscoreDict = {} #keep unique id: query id + taxon id
    with open(blastFile) as f:
        for line in f:
            cells = line.rstrip().split('\t')
            #cells[10] = float(cells[10]) #convert string to sci number
            id = cells[0] + "|" + cells[13] #query id + tax id
            if len(cells) < 14: #for missing taxonomy
                cells.extend(["","0"])            
            if(id not in blastDict or float(bitscoreDict[id]) < float(cells[11])): #only keep top one (the highest bitscore) for a query
                #blastDict.setdefault(cells[0], []).append([cells])
                blastDict[cells[0]] = cells
                bitscoreDict[id] = cells[11]  #bitscore
    return blastDict


#################################### merge blastx and blastn  ####################################
def mergeBlast(blastnDict, blastxDict, contig_file):
    """
    Merge blastn and blastx
    """
    outDict = {}
    idDict = {}

    #put contig sequences in a dict
    contig_dict = SeqIO.index(contig_file, "fasta")
    #print(contig_dict)
    #get all sequence id
    for id in blastnDict:
        idDict[id] = 1  #only in blastn
    for id in blastxDict:
        if id in idDict:
            idDict[id] = 2  #in both blastn and blastx
        else:
            idDict[id] = 0  #only in blastx
    
    #merge blastn and blastx
    empty = "\t"*16
    for id in idDict:
        seq = contig_dict[id]
        if idDict[id] == 2:  #id in both blastx and blastn
            outDict[id] = "\t".join(blastxDict[id]) + "\t" + "\t".join(blastnDict[id]) + "\t" + seq.seq
        elif idDict[id] == 1: #id in only blastn
            outDict[id] = empty + "\t".join(blastnDict[id]) + "\t" + seq.seq
        else:  #id in only blastx
            outDict[id] = "\t".join(blastxDict[id]) + "\t" + empty + seq.seq

    contig_dict.close()
    
    return outDict

#################################### Get arguments function ####################################
def parseArguments():
    """
    Get arguments
    """
    parser = argparse.ArgumentParser(description = "Select top reference from blastn and blastx result.")
    parser.add_argument('-n','--blastn',dest='blastn_file', required=True, default="", help='Blastn against nt database result')
    parser.add_argument('-x','--blastx',dest='blastx_file', required=True, default="", help='Blastx against protein database (nr or virus RVDB) result')
    parser.add_argument('-c','--contig',dest='contig_file', required=True, default="", help='Contig file from assembler output, ex. contigs.txt ')
    parser.add_argument('-o','--output',dest='output_file', required=True, default="", help='output file for merged blastn and blastx' )
    parser.set_defaults(append=False)
    return parser.parse_args()

def main():
    ### Input arguments
    options = parseArguments()
    blastn_file = options.blastn_file
    blastx_file = options.blastx_file
    contig_file = options.contig_file
    output_file = options.output_file

    blastnDict = getBlastn(blastn_file)
    #print(blastnDict)
    blastxDict = getBlastx(blastx_file)

    outDict = mergeBlast(blastnDict, blastxDict, contig_file)
    #output selected blast result
    fout = open(output_file, 'w')
    fout.write("Contig_blastx\tReference_blastx\tPercent_identity_blastx\tAlignment_length_blastx\tContig_length_blastx\tReference_length_blastx\tAlign_contig_start_blastx\tAlign_contig_end_blastx\tAlign_ref_start_blastx\tAlign_ref_end_blastx\tE-value_blastx\tBitscore_blastx\tReference_name_blastx\tTaxonomy_id_blastx\tSpecies_blastx\tTaxonomy_path_blastx\t")
    fout.write("Contig_blastn\tReference_blastn\tPercent_identity_blastn\tAlignment_length_blastn\tContig_length_blastn\tReference_length_blastn\tAlign_contig_start_blastn\tAlign_contig_end_blastn\tAlign_ref_start_blastn\tAlign_ref_end_blastn\tE-value_blastn\tBitscore_blastn\tReference_name_blastn\tTaxonomy_id_blastn\tSpecies_blastn\tTaxonomy_path_blastn\tContig\n") 
    for id in outDict:
        fout.write("%s\n" % (outDict[id]))
    fout.close()

    
if __name__ == '__main__':
	main()
