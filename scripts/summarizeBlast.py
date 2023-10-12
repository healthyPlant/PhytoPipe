#! /usr/bin/env python
######################################################################
#Alex Hu <xiaojun.hu@usda.gov>
#Updated: 02/20/2023
# The script summarize blasts result from contigs against NCBI nt and nr, also include classification results from Kraken and Kaiju
# The input file is sample.blastnt.txt or sample.blastnr.txt
# The output file is sample.blastnt.summary.txt or sample.blastnr.summary.txt
#Usage: python summarizeBlast.py [input file] [output file]
######################################################################

import sys

if len(sys.argv) != 3:
    print("Usage: python summarizeBlast.py [input file] [output file] ")
    print("Argument missed, exit!")
    sys.exit()

blastFile = sys.argv[1] #input balst table
outFile = sys.argv[2] #output summary 

sampleName = blastFile.split("/")[-1].split(".")[0] #sample name

#read blast file
#ContigName	Reference	PercentIdentity	AlignLength	ContigLength	RefLen	ContigAlignStart	ContigAlignEnd	RefAlignStart	RefAlignEnd	Evalue	Bitscore	RefTitle	RefTaxonId	KrakenRead	KrakenSpecies	KrakenMicrobialGroup	Contig
#ContigName	Reference	PercentIdentity	AlignLength	ContigLength	RefLen	ContigAlignStart	ContigAlignEnd	RefAlignStart	RefAlignEnd	Evalue	Bitscore	RefTitle	RefTaxonId	KaijuRead	Species	TaxonPath	Contig
speciesDict = {}
with open(blastFile) as f:
    next(f) #skip header
    for line in f:
        cells = line.split("\t")
        contigLength = int(cells[4])
        classifiedRead = cells[14]
        species = cells[15]
        microbialGroup = cells[16]
        if species != "nan":
            if species not in speciesDict:
                speciesDict.setdefault(species,[])
            speciesDict[species].append([contigLength, classifiedRead, microbialGroup])

#print(speciesDict)
speciesSumDict = {}
for species in speciesDict:
    contigNum = len(speciesDict[species])
    if speciesDict[species][0][1] == 'nan':
        classifiedRead = 0
    else:
        classifiedRead = int(speciesDict[species][0][1])
    microbialGroup = speciesDict[species][0][2]
    if microbialGroup.startswith("cellular organisms;Eukaryota;Opisthokonta;Fungi;"):
        microbialGroup = "Fungi"
    elif microbialGroup.startswith("cellular organisms;Bacteria;"):
        microbialGroup = "Bacteria"
    elif microbialGroup.startswith("Viruses;"):
        microbialGroup = "Viruses"
    elif microbialGroup.startswith("cellular organisms;Eukaryota;Sar;Stramenopiles;Oomycota;"):
        microbialGroup = "Oomycetes"
    elif microbialGroup.startswith("cellular organisms;Archaea;"):
        microbialGroup = "Archaea"
    elif microbialGroup.startswith("cellular organisms;Eukaryota;"):
        microbialGroup = "Eukaryota"

    #get longest contig size 
    maxContig = 0
    for item in speciesDict[species]:
        if item[0] > maxContig:
            maxContig = item[0]
    
    speciesSumDict[species] = [sampleName, microbialGroup, species, classifiedRead, contigNum, maxContig] 

#output
header = "Sample\tTypes of pathogens\tPathogen name\tRead number by classification\tContig number\tLongest Contig size"
fout = open(outFile,'w')
fout.write(header + "\n")

for key, value in sorted(speciesSumDict.items(), key=lambda e: e[1][3], reverse=True):  #output sorted by classifiedRead
    if value[1] != 'nan' or value[4] > 1 :  #don't output unknown microbialGroup and just one contig 
        fout.write("\t".join(map(str, value)) + "\n")
fout.close()
