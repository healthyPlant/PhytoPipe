#! /usr/bin/env python
######################################################################
#Alex Hu <xiaojun.hu@usda.gov>
#Updated: 01/06/2023
#This program compare a consensus sequence and a contig, both are close to the same referenece.
#If the conitg length is more than the consensus sequence length (without N) 90%, then use the contig not the consensus for the final report  
#
#Required Parameters:
#   -r, --ref X..............................reference name
#   -s, --sample X...........................sample name
#   -c, --consensus X........................conesensus sequence file, ex. sample.refNamex.consensus.N.fasta
#   -b, --blast X............................Blast result for selected references
#   -g, --contig X...........................de novo assembled contigs, ex. sample/contigs.fasta
#   -o, --output X...........................output folder
######################################################################

import os,sys,argparse 
from Bio import SeqIO

#################################### Get reference and contig association function ####################################
#Contig	Reference	Percent_identity	Alignment_length	Contig_length	Reference_length	Align_contig_start	Align_contig_end	Align_ref_start	Align_ref_end	E-value	Bitscore	Reference_name	Taxonomy_id	Species	Taxonomy_path	blastType
#NODE_6_length_663_cov_38.759386	NC_001642.1	98.793	663	663	6366	1	663	431	1093	0.0	1181	NC_001642.1 Bamboo mosaic virus, complete genome	35286	 Bamboo mosaic virus	Viruses; Riboviria; Tymovirales; Alphaflexiviridae; Potexvirus; Bamboo mosaic virus;	Blastn
#TRINITY_DN201_c0_g2_i1  NC_015324.1     96.115  489     489     2824    1       489     2700    2212    0.0     798       NC_015324.1 Sweet potato golden vein associated virus, complete genome  688540   Sweet potato golden vein associated virus        Viruses; Geminiviridae; Begomovirus; unclassified Begomovirus; Sweet potato golden vein associated virus; Blastn
def getRelation(refFile):
    """
    Get contig and reference relationship from $sample.selectedRef.txt, each reference is realted to the longest contig. Identity and length of the contig are kept
    """
    refDict = {}
    lengthDict = {} #for contig length
    with open(refFile) as f:
        next(f) #skip the header
        for line in f:
            cells = line.rstrip().split('\t')
            length = int(cells[4]) #contig length
            #identity = cells[2]
            id = cells[1] #reference name
            if "acc|" in cells[1]:
                id = cells[1].split("|")[-2]
            #elif "." in cells[1]:
            #    id = cells[1].split(".")[0]
            if(id not in refDict or lengthDict[id] < length): #only keep the longest contig for a ref
                refDict[id] = [cells[0], length]
                lengthDict[id] = length  
    return refDict

#################################### Get consensus length without N function ####################################
def getConsensusLength(filename):
    """
    Get consensus length without N
    """
    record = SeqIO.read(filename, "fasta") #Bio.SeqIO.parse(...) used to loop over all the records in a file as SeqRecord objects. SeqIO.read() used on files which contain just one record
    countN = record.seq.count('N')
    length = len(record.seq) - countN

    return length

#################################### compare a contig and a consensus sequence ####################################
def compareSeq(consensusLenth, refDict, refName):
    """
    Compare a contig and a consensus sequence, if the conitg length is longer than the consensus sequence length (without N) 90%, then use the contig not the consensus for the final report   
    """    
    #print(consensusLenth)
    #print(refDict)

    contigLength = refDict[refName][1]

    #compare length
    if contigLength > consensusLenth * 0.9:
        return refDict[refName][0] #return contig name
    else:
        return 'No'


#################################### Output contigs to finalAssembly ####################################
def outputContig(consensusLenth, refDict, refName, contigsFile, outputFile, sampleName):
    """
    Output a contig
    """
    #put spade contigs in dict
    record_dict = SeqIO.index(contigsFile, "fasta")

    contigName = compareSeq(consensusLenth, refDict, refName)

    if contigName != 'No':
        #extract a contig
        contig = record_dict[contigName]
        #print(contig)
        #outName = outputFolder + "/" + sampleName + "." + refName + ".contig.fasta"
        contig.id = sampleName + "." + refName + ".contig"
        output_handle = open(outputFile, "w")
        SeqIO.write(contig, output_handle, "fasta")
        output_handle.close()

    record_dict.close()

#################################### Get arguments function ####################################
def parseArguments():
    """
    Get arguments
    """
    parser = argparse.ArgumentParser(description = "Compare a contig and a consensus sequence.")
    parser.add_argument('-r','--ref',dest='ref_name', required=True, default="", help='Genbank reference which close to both the contig and conesensus')
    parser.add_argument('-s','--sample',dest='sample_name', required=True, default="", help='Sample name')
    parser.add_argument('-c','--consensus',dest='consensus_file', required=True, default="", help='conesensus sequence file, ex. sample.refNamex.consensus.N.fasta')
    parser.add_argument('-b','--blast',dest='blast_file', required=True, default="", help='Blast result for selected references')
    parser.add_argument('-g','--contig',dest='contig_file', required=True, default="", help='De novo assembled contigs, i.e. sample/contigs.fasta' )
    parser.add_argument('-o','--output',dest='output_file', required=True, default="", help='Output a conitg')
    parser.set_defaults(append=False)
    return parser.parse_args()

def main():
    ### Input arguments
    options = parseArguments()
    refName = options.ref_name
    sampleName = options.sample_name
    consensusFile = options.consensus_file
    refBlastFile = options.blast_file
    contigsFile = options.contig_file
    outputFile = options.output_file

    #read files
    refDict = getRelation(refBlastFile)

    #get a consensus length
    consensusLenth = getConsensusLength(consensusFile)

    #copy contigs to final assebly folder
    outputContig(consensusLenth, refDict, refName, contigsFile, outputFile, sampleName)

if __name__ == '__main__':
    main()
