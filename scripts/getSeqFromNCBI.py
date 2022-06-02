#! /usr/bin/env python
######################################################################
#Alex Hu <xiaojun.hu@usda.gov>
#Updated: 02/18/2021
#This program extracts a sequence from NCBI web by gb
#
#Required Parameters:
#   -n, --name X...........................Genbank id file, one per line 
#   -o, --output X.........................output folder, a sequence per file 
######################################################################
#https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch

import os,sys,argparse 
from Bio import Entrez 
from Bio import SeqIO
#from urllib2 import HTTPError  #python2
from urllib.error import HTTPError #python3
import time

#################################### Get arguments function ####################################
def parseArguments():
    """
    Get arguments
    """
    parser = argparse.ArgumentParser(description = "Download sequences from NCBI by Genbank id.")
    parser.add_argument('-n','--name',dest='name_file', required=True, default="", help='Genbank id file, one per line')
    parser.add_argument('-o','--output',dest='output_folder', required=True, default="", help='NCBI fasta sequences output folder')
    parser.set_defaults(append=False)
    return parser.parse_args()

#################################### Download sequences from NCBI function ####################################
def downloadSeq(nameFile, outFolder):
    """
    Download sequences from NCBI using Entrez.efetch
    """    
    Entrez.email = "somebody@example.com" # Always tell NCBI who you are.
    fin = open(nameFile, 'r')
    #fout = open(outFile, "w") 
    for line in iter(fin):
        accession = line.strip()
        try:
            net_handle = Entrez.efetch(db="nucleotide",id=accession,rettype="fasta", retmode="text")  #rettype="gb",
            time.sleep(10)
        except HTTPError:
            time.sleep(60)
            net_handle = Entrez.efetch(db="nucleotide",id=accession,rettype="fasta", retmode="text")  #rettype="gb",
        #ouput reference squences
        #outfilename = "%s.%s.refseq.fasta" % (outName,accession)
        outfilename = "%s/%s.fasta" % (outFolder, accession)
        if not os.path.isdir(outFolder): #Check if a directory exists
            os.mkdir(outFolder) 
        fout = open(outfilename, "w") 
        #records = Entrez.read(handle)
        fout.write(net_handle.read())
        fout.close() 
        time.sleep(1) # to make sure not many requests go per second to ncbi
        net_handle.close()
    fin.close()
    #fout.close()

def main():
    ### Input arguments
    options = parseArguments()
    nameFile = options.name_file
    outFolder = options.output_folder

    #call download function
    downloadSeq(nameFile, outFolder)

if __name__ == '__main__':
	main()
