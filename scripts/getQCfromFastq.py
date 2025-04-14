#! /usr/bin/env python
#Get a raw fastq QC information function #
import gzip
import sys, os

def getQualFromFastq(fastq_file):
    """
    Get raw quality information from a fastq.gz file 
    """

    #@ is the read start 
    # replace ASCII quality with integer
    #ord()	Converts a character to an integer
    #ord(str(quality)) - 33
    totalBases = 0
    totalQual = 0
    q30_count = 0
    lineno = 0
    for line in gzip.open(fastq_file, "rt"):
        lineno += 1
        if lineno%4 == 0:
            line = line.strip() 
            qstr = list(line)
            for q in qstr:
                qual = ord(q) - 33
                if qual >= 30:
                    q30_count += 1
                totalQual += qual
            totalBases += len(line)
    totalRead = int(lineno/4)
    avgQ30 = 100 * float(q30_count)/float(totalBases)
    avgQual = float(totalQual)/float(totalBases)
    sample=os.path.basename(fastq_file).split(".")[0]
    
    avgQualSet = (sample, str(totalRead), "{:.0f}".format(totalBases/1000000), "{:.2f}".format(avgQ30), "{:.2f}".format(avgQual))
    return(avgQualSet)

fastq_file = sys.argv[1]
avgQualSet = getQualFromFastq(fastq_file)
print(f'Sample: {avgQualSet[0]} \n Total reads: {avgQualSet[1]} \n  Yield(Mbases): {avgQualSet[2]} \n  Percent>=Q30Bases: {avgQualSet[3]} \n MeanQualityScore: {avgQualSet[4]}')