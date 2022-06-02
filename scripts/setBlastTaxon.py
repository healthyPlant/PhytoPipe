#! /usr/bin/env python
######################################################################
#Alex Hu <xiaojun.hu@usda.gov>
#Updated: 02/18/2021
#This program adds taxonomy infromation to blastn/x result
#
#Required Parameters:
#   -t, --taxon X.......................taxon file which have gb id, taxon id, full taxon path 
#   -b, --blast X.......................blastn/x table file
#   -o, --output X......................output file 
######################################################################

import os, sys, argparse

#################################### Get arguments function ####################################
def parseArguments():
    """
    Get arguments
    """
    parser = argparse.ArgumentParser(description = "Add taxonomy to blast result.")
    parser.add_argument('-t','--taxon',dest='taxon_file', required=True, default="", help='Taxonomy file having gb id, taxon id, full taxon path (kindom,phylum,class,order,family,genus,species)')
    parser.add_argument('-b','--blast',dest='blast_file', required=True, default="", help='Blastx or blastn table fromat result (-outfmt 6)')
    parser.add_argument('-o','--output',dest='output_file', required=True, default="", help='output file for adding taxon id, taxon path to the blast result' )
    parser.set_defaults(append=False)
    return parser.parse_args()

#################################### Get taxonomy ####################################
#AC_000001.1	114424	Viruses; Adenoviridae; Mastadenovirus; Ovine mastadenovirus A;
#AAA02484.1	35330	Viruses; Riboviria; Reoviridae; Sedoreovirinae; Orbivirus; Bluetongue virus; Bluetongue virus 13;
def getTaxon(taxonFile):
    """
    Get taxonomy information for a genbank id
    """
    taxonDict = {}
    with open(taxonFile) as f:
        for line in f:
            cells = line.rstrip().split('\t')
            species = "NA"
            if len(cells) ==3:
                if ";" in cells[2]:
                    species = cells[-1].split(";")[-2]
                    taxonDict[cells[0]] = [cells[1], species, cells[2]]
            else:
                taxonDict[cells[0]] = [cells[1], species, "NA"]
    return taxonDict

#################################### Set taxonomy ####################################
#acc|GENBANK|ACT34024.1|GENBANK|GQ268816|PsaC
def setTaxon(blastFile, taxonFile):
    """
    Set taxonomy information for a genbank id
    """
    taxonDict = getTaxon(taxonFile)
    blastDict = {}
    #blastAnn = []
    with open(blastFile) as f:
        for line in f:
            cells = line.rstrip().split('\t')
            acc = cells[1]  #reference genbank id
            if  "acc|" in cells[1]:
                acc = cells[1].split("|")[2]
            if acc in taxonDict:
                if(cells[0] not in blastDict or float(blastDict[cells[0]][11]) < float(cells[11])): #only keep top one (the highest bitscore) for a query
                    blastDict[cells[0]] = cells + taxonDict[acc] #cells.extend(taxonDict[acc[2]])
            else:
                if(cells[0] not in blastDict or float(blastDict[cells[0]][11]) < float(cells[11])): #only keep top one (the highest bitscore) for a query
                    blastDict[cells[0]] = cells + ['NA', 'NA', 'NA']
    return blastDict


def main():
    ### Input arguments
    options = parseArguments()
    taxon_file = options.taxon_file
    blast_file = options.blast_file
    output_file = options.output_file

    blastDict = setTaxon(blast_file, taxon_file)
    fout = open(output_file, 'w')
    #print(blastxDict)
    for id in blastDict:
        fout.write("%s\n" % "\t".join(blastDict[id]))

    fout.close()


if __name__ == '__main__':
	main()

