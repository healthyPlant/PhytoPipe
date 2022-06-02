#! /usr/bin/env python
######################################################################
#Alex Hu <xiaojun.hu@usda.gov>
#Updated: 02/18/2021
#This program splits a local contig/consensus blastn result to different files
#
#Required Parameters:
#   -s, --sequence X.......................a fasta file having consensus/contig sequences
#   -b, --blast X..........................Blast table result from local blast nt
#   -o, --output X.........................output prefix 
######################################################################

import os,sys,argparse 
from Bio import SeqIO

##################################### Get blastn result ####################################
#qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore stitle qcovs
def getBlast(blastFile):
    """
    Get blast result, keep the top one (highest bitscore) for a query, calculate coverage 
    s_cov = 100 * (s_end - s_start) / s_len
    """
    blastDict = {}
    covDict = {}
    with open(blastFile) as f:
        for line in f:
            cells = line.rstrip().split('\t')
            if(cells[0] not in blastDict or float(blastDict[cells[0]][11]) < float(cells[11])): #only keep top one (the highest bitscore) for a query
                blastDict[cells[0]] = cells
                #calculate coverage
                #s_start, s_end, s_len = map(float, (cells[8], cells[9], cells[5]))
                #q_cov = 100.0 * abs(s_end - s_start) / s_len
                #q_cov = format(q_cov, '.2f')
                q_cov = cells[-1] #using qcovs: Query Coverage Per Subject (for all HSPs)

                #Coverage,Evalue,Identity,Description
                desc = cells[1] + " " + cells[12]  #concatenate id and description
                covDict[cells[0]] = [q_cov,cells[10],cells[2],desc]
    
    return covDict


#################################### Add sequences to blast result ####################################
def AddSeq(seqFile, covDict, failedFile):
    """
    Add sequence to blast result
    """
    #put fasta sequences in dict
    record_dict = SeqIO.index(seqFile, "fasta")
    outHandle = None
    if failedFile:
        outHandle = open(failedFile, "w")

    for seqName in record_dict:
        #extract a contig from assembly contigs
        seq = record_dict[seqName]
        #print(seqName)
        if seqName in covDict:
            covDict[seqName].append(seq.seq)
            print(seqName, " blast success")
        else:
            #covDict[seqName] = ['nan','nan','nan','nan', seq.seq]
            if failedFile:
                SeqIO.write(seq, outHandle, "fasta")
            print(seqName, " blast failed")
    if failedFile:
        outHandle.close()

    record_dict.close()

#################################### Add sequences to blast result ####################################
def output2File(covDict, output, blastType):
    """
    Add sequence to blast result
    """
    for seqName in covDict:
        outFile = output + "/" + seqName + "." + blastType + ".txt"
        fout = open(outFile, 'w')
        fout.write("Coverage\tEvalue\tIdentity\tDescription\tSequence\n")
        fout.write("\t".join(map(str,covDict[seqName])) + "\n")
        fout.close()          

#################################### Get arguments function ####################################
def parseArguments():
    """
    Get arguments
    """
    parser = argparse.ArgumentParser(description = "Split blastn results to different files")
    parser.add_argument('-s','--sequence',dest='sequence_file', required=True, default="", help='Consensus and contig sequence file')
    parser.add_argument('-b','--blast',dest='blast_file', required=True, default="", help='Blast result from local blast nt')
    parser.add_argument('-f','--fail',dest='blast_failed_seq', required=False, default="", help='Blast failed sequences')  
    parser.add_argument('-t','--type',dest='blast_type', required=True, default="blastn", help='Blast type:blastn or blastx')   
    parser.add_argument('-o','--output',dest='output', required=True, default="", help='Output prefix')
    parser.set_defaults(append=False)
    return parser.parse_args()

def main():
    ### Input arguments
    options = parseArguments()
    seqFile = options.sequence_file
    blastFile = options.blast_file
    failedFile = options.blast_failed_seq
    blastType = options.blast_type
    output = options.output

    #get blastn results
    blastDict = getBlast(blastFile)

    #add sequences to blastDict
    print("Running ", blastType)
    if failedFile:
        AddSeq(seqFile, blastDict, failedFile )
    else:
        AddSeq(seqFile, blastDict, "" )

    #split blastn result to different files, each sequence a file
    output2File(blastDict, output, blastType)

if __name__ == '__main__':
	main()
