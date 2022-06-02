#! /usr/bin/env python
######################################################################
#Alex Hu <xiaojun.hu@usda.gov>
#Updated: 02/18/2021
#This program remote Blast NCBI
#
#Required Parameters:
#   -i, --input X..........................fasta file
#   -b, --blast X..........................Blast type ex. blastn, blastp
#   -d, --blastdb X........................NCBI blast database ex. nt, nr
#   -o, --output X.........................Blast report from NCBI 
#https://python.hotexamples.com/examples/Bio.Blast/NCBIWWW/qblast/python-ncbiwww-qblast-method-examples.html
######################################################################

import os
import argparse
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
from urllib.error import HTTPError
import time
import urllib

##Do not contact the server more often than once every 10 seconds.
#################################### Remote blast NCBI function ####################################
def blastNCBI(seq, blastType, blastDb, blastFile):
    """
    Remote blast NCBI and save xml result
    """
    try:
        blast_handle = NCBIWWW.qblast(blastType, blastDb, seq,  descriptions=10,alignments=10, expect=0.001) #("blastn", "nt", seq) ("blastx", "nr", seq)
        #res_handle = NCBIWWW.qblast(prog, db, query, entrez_query = myEntrezQuery)
        # save the result:
        savefile = open(blastFile, "w")
        savefile.write(blast_handle.read())
        savefile.close()
        blast_handle.close()
    except HTTPError:
        print('Error in remote blast!')
    except ValueError:
        print('Error in remote blast!')
    except (ConnectionResetError, urllib.error.URLError):
        print('Remote blast connection error!')

    
#################################### Remote blast NCBI function ####################################
def parseBlast(blastFile):
    """
    Parse Blast xml, output coverage, e-value, identity, description
    """
    try:
        resultHandle = open(blastFile, "r")
        #Return a single Blast record (assumes just one query).
        blastRecord = NCBIXML.read(resultHandle)
        result = []        
        desc = blastRecord.alignments[0].title #hit description, accession, hit_id, hit_def
        #print(blastRecord.alignments[0].length) #hit length
        evalue = blastRecord.alignments[0].hsps[0].expect #e-value, score, query, match, sbjct
        #print(blastRecord.alignments[0].hsps[0].identities) #Number of identities
        #print(blastRecord.alignments[0].hsps[0].align_length) #Length of the alignment
        identity = blastRecord.alignments[0].hsps[0].identities / blastRecord.alignments[0].hsps[0].align_length
        coverage = blastRecord.alignments[0].hsps[0].align_length / blastRecord.alignments[0].length
        if coverage > 1:
            coverage = 1.0
        #print(round(coverage*100,2), '{:.2e}'.format(evalue), round(identity*100,2), desc, record.seq)
        result = [round(coverage*100,2), '{:.2e}'.format(evalue), round(identity*100,2), desc]
    except (IndexError, ValueError,FileNotFoundError):
        print("NCBI blast fail")
        result = ['nan', 'nan', 'nan','BlastFail'] 
    
    return result
    
#################################### Get arguments function ####################################
def parseArguments():
    """
    Get arguments
    """
    parser = argparse.ArgumentParser(description='Generate text report')
    parser.add_argument("-i", "--input", dest='fasta_file', required=True, help="A fasta file")
    parser.add_argument("-o", "--output", dest='output_file', required=True, help="Blast report from NCBI")
    parser.add_argument("-b", "--blast", dest='blast_type', required=True, help="Blast type, like blastn, blastp")
    parser.add_argument("-d", "--blastdb", dest='blast_db', required=True, help="NCBI blast database, like nt, nr")
    parser.set_defaults(append=False)
    return parser.parse_args()


def main():
    ### Input arguments
    options = parseArguments()
    fastaFile = options.fasta_file
    outFile = options.output_file
    blastType = options.blast_type
    blastDb = options.blast_db
    
    blastFile = outFile + ".xml"
    record = SeqIO.read(fastaFile, format="fasta")
    
    blastNCBI(record.seq, blastType, blastDb, blastFile)
    blastOut = parseBlast(blastFile)
    #print(len(blastOut))
    
    #in case blast NCBI fail
    #retry = 0
    #while (len(blastOut) == 0 and retry < 3):
    time.sleep(10)
    blastNCBI(record.seq, blastType, blastDb, blastFile)
    blastOut = parseBlast(blastFile)
    #retry += 1
    #print("Retry: ", retry)
                                      
    #output result
    fout = open(outFile, 'w')
    fout.write("Coverage\tEvalue\tIdentity\tDescription\tConsensus\n")
    fout.write("\t".join(map(str,blastOut)))
    fout.write("\t" + str(record.seq) + "\n")
    #delete blastFile
    if os.path.exists(blastFile):
      os.remove(blastFile)
    
if __name__ == "__main__":
    main()
    