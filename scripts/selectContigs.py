#! /usr/bin/env python
######################################################################
#Alex Hu <xiaojun.hu@usda.gov>
#Updated: 02/18/2021
#This program selects a better contig after comparing de novo assembled contigs and Viral-NGS refined contigs.
# The novel virus contig may be missed by Viral_NGS order_and_orient based on the reference genome
#
#Required Parameters:
#   -r, --ref X..............................reference name file, one genbank id per line, ex. sample.refNamex.txt'
#   -b, --blast X............................Blast result for selected references
#   -c, --contig X...........................de novo assembled contigs, ex. sample/contigs.fasta
#   -p, --pseudo X...........................Viral-NGS refined pseudo-contigs folder, i.e. pseudoContig/sample
#   -o, --output X...........................Final assembly folder
######################################################################

import os,sys,argparse 
from Bio import SeqIO
import shutil

#################################### Get referennce name function ####################################
def getName(nameFile):
    """
    Get reference name from $sample.refNames.txt
    """    
    nameDict = {} #keep unique id
    with open(nameFile) as f:
        for line in f:
            line = line.rstrip()
            nameDict[line] = 1
    return nameDict

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
            length = int(cells[4])
            identity = cells[2]
            id = cells[1]
            if "acc|" in cells[1]:
                id = cells[1].split("|")[-2]
            #elif "." in cells[1]:
            #    id = cells[1].split(".")[0]
            if(id not in refDict or lengthDict[id] < length): #only keep the longest contig for a ref
                refDict[id] = [cells[0], length, identity]
                lengthDict[id] = length  
    return refDict

#################################### Get refined contig length without N function ####################################
def getContigLength(sample, nameDict):
    """
    Get refined contig length
    """
    lenDict = {} #for contig length
    for refName in nameDict:
        filename = sample + "/" + refName + ".pseudoContig.fasta"  #full path
        if not os.path.isfile(filename):
             lenDict[refName] = [filename, 1]
        else:
            for record in SeqIO.parse(filename, "fasta"):
                countN = record.seq.count('N')
                length = len(record.seq) - countN
                lenDict[refName] = [filename, length]
    return lenDict

#################################### compare contigs and refined contigs ####################################
def compareContigs(lenDict, refDict):
    """
    Compare Spades contigs and Viral_NGS refined contigs, if identity with reference < 80%, treat is as a novel virus and keep the original contig 
    """    
    contigNameDict = {}
    #print(lenDict)
    #print(refDict)
    for refName in lenDict:
        lenDif = lenDict[refName][1] - int(refDict[refName][1])  #compare contigs length, viral-ngs - spades
        if float(refDict[refName][2]) <= 80.0: #identity <= 80, it is a novel virus, keep the Spades contig
            contigNameDict[refName] = refDict[refName][0] #contig name in contigs.fasta 
        else:
            if lenDif >= 20: #if refined contig improve 20bp, use it
                contigNameDict[refName] = lenDict[refName][0]
            else:
                contigNameDict[refName] = refDict[refName][0]
    #print(contigNameDict)
    return contigNameDict

#################################### Output contigs to finalAssembly ####################################
def outputContigs(contigNameDict, spadesContig, outputFolder):
    """
    Output contigs
    """
    #put spade contigs in dict
    record_dict = SeqIO.index(spadesContig, "fasta")

    for refName, contigName in contigNameDict.items():
        outName = outputFolder + "/" + refName + ".contig.fasta"
        if "pseudoContig" in contigName: #if a pseudoContig, copy it
            #outName = outputFolder + "/" + os.path.basename(contigName) 
            shutil.copyfile(contigName, outName) #copy 
        else: #if a Spades contig, extract it
            #extract a contig
            contig = record_dict[contigName]
            #print(contig)
            output_handle = open(outName, "w")
            SeqIO.write(contig, output_handle, "fasta")
            output_handle.close()
    record_dict.close()

#################################### Get arguments function ####################################
def parseArguments():
    """
    Get arguments
    """
    parser = argparse.ArgumentParser(description = "Get final assemblies from Spades contigs and Viral-NGS refined contigs.")
    parser.add_argument('-r','--ref',dest='ref_name', required=True, default="", help='Genbank reference name file, one per line, i.e. sample.refNamex.txt')
    parser.add_argument('-b','--blast',dest='blast_file', required=True, default="", help='Blast result for selected references')
    parser.add_argument('-c','--contig',dest='contig_file', required=True, default="", help='Spades assembled contigs, i.e. sample/contigs.fasta' )
    parser.add_argument('-p','--pseudo',dest='pseudo_folder', required=True, default="", help='Viral-NGS refined pseudo-contigs folder, i.e. pseudoContig/sample' )
    parser.add_argument('-o','--output',dest='output_folder', required=True, default="", help='Final assembly folder')
    parser.set_defaults(append=False)
    return parser.parse_args()

def main():
    ### Input arguments
    options = parseArguments()
    refNameFile = options.ref_name
    refBlastFile = options.blast_file
    spadesContigFile = options.contig_file
    refineContigFile = options.pseudo_folder
    outFolder = options.output_folder

    #read files
    nameDict = getName(refNameFile)
    refDict = getRelation(refBlastFile)

    #get refined contig length
    lenDict = getContigLength(refineContigFile, nameDict)

    #get a better contig between sapdes contigs and viral-ngs contigs
    contigNameDict = compareContigs(lenDict, refDict)

    #copy contigs to final assebly folder
    outputContigs(contigNameDict, spadesContigFile, outFolder)

if __name__ == '__main__':
	main()
