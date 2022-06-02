#! /usr/bin/env python
######################################################################
#Alex Hu <xiaojun.hu@usda.gov>
#Updated: 02/18/2021
#This program gets novel viruses read mapping information for a contig
#
#Required Parameters:
#Required Parameters:
#   -w, --workdir X........................working directory
#   -o, --output X.........................output file 
#######################################################################!/usr/bin/env python

import glob
import os
import sys
import argparse

#################################### Get mapping coverage information function ####################################
def getCoverage(covFile):
    """
    Get coverage percentage and average depth
    """
    count = count0 = total = 0
    with open(covFile) as f:
        for line in f:
            line = line.strip()
            cells = line.split()
            cov = float(cells[-1])
            count += 1
            total += cov
            if cov == 0:
                count0 += 1
        coverage = (count - count0)/float(count) * 100
        depth = total/count
        coverage = float("{:.2f}".format(coverage))
        depth = float("{:.1f}".format(depth))
    return (coverage, depth, count)


#################################### Get mapping to reference information function ####################################
def getMapInfo(bamFile, coverageFile):
    """
    Summary mapping information:reference length, reference title, mapped reads, coverage, rpkm, 
    """

    #get reference and its mapping information
    #bamFile = 'mapping/map2Ref/NC_Apples_3_S3/NC_028132.sorted.bam'
    #bamFile = mapDir + '/' + sample + '/' + ref + '.sorted.bam'
    cmd = 'samtools view -c -F 260 ' + bamFile
    mappedReads = os.popen(cmd).read().strip() #get sys comand return value, os.system only returns the error value.
    
    #get coverage
    #coverageFile = 'mapping/map2Ref/NC_Apples_3_S3/NC_028132.coverage.txt'
    #coverageFile = mapDir + '/' + sample + '/' + ref + '.coverage.txt'
    (coverage, depth, length) = getCoverage(coverageFile)
    
    mapArr = map(str, [length, mappedReads, coverage, depth])
    
    return mapArr

#################################### Get arguments function ####################################
def parseArguments():
    """
    Get arguments
    """
    parser = argparse.ArgumentParser(description='Generate text report')
    parser.add_argument("-w", "--workdir", dest='workdir', required=True, help="Working directory, ex. /ppq/data2/pgqp_pipeline/Run45")
    parser.add_argument("-o", "--output", dest='output_file', required=True, help="output file name, ex. novelVirusMapping.txt")
    parser.set_defaults(append=False)
    return parser.parse_args()


def main():
    ### Input arguments
    options = parseArguments()
    workdir = options.workdir  #'/ppq/data2/pgqp_pipeline/Run45'
    outFile = options.output_file #workdir + '/report/novelVirusMapping.txt'

    #get samples and reference name
    #novelVirus/map2Contig/OCT19-Q50416_S4
    
    mapDir = workdir + '/novelVirus/map2Contig/'
    fasFiles = glob.glob(mapDir + '*.fasta')
    #print(fasFiles)
    smpRefDict = {}
    for fname in fasFiles:
        names = os.path.basename(fname).split(".")
        smpRefDict.setdefault(names[0],[]).append(names[1])
    #print(smpRefDict)
    
    #get mapping coverage
    fout = open(outFile, 'w')
    header = 'Sample\tRefId\tRefLength\tMappedReads\tPercentContigCovered\tMeanCoverage\tRefTitle'
    fout.write(header + '\n')
    for sample in smpRefDict:
        for ref in smpRefDict[sample]:
            coverageFile = mapDir + sample + "/" + ref + '.contig.coverage.txt' 
            bamFile = mapDir + sample + "/" + ref + '.contig.sorted.bam'
            refFile = workdir + "/mapping/ref/" + sample + "/" + ref + ".fasta"
            if os.path.exists(coverageFile) and os.path.exists(bamFile):
                mapArr = getMapInfo(bamFile, coverageFile)
                #print(sample, ref+".contig", " ".join(mapArr))
                refTitle = 'NA'
                if os.path.exists(refFile):
                    refTitle = open(refFile).readline().rstrip().replace(">","")
                fout.write(sample + '\t' + ref+'.contig\t' + '\t'.join(mapArr) + '\t' + refTitle + '\n')
    fout.close()        

if __name__ == "__main__":
    main()
