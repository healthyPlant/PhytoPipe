#! /usr/bin/env python
######################################################################
#Alex Hu <xiaojun.hu@usda.gov>
#Updated: 02/18/2021
#This program generates a txt report, including read mapping statistics, consensus/contig blast results and sequences
#The report has the following contents:
#1.	Total number of reads obtained and approximate read length
#2.	Number of reads mapping to the virus genome, the % of genome covered and the mean read coverage over the genome
#3.	The virus sequence that was used as reference
#4.	The closest nucleotide match from GenBank after BLASTn

#Required Parameters:
#   -w, --workdir X........................working directory
#   -a, --acronym X........................Acronym file from ICTV, ex. ICTV_virus_acronym2019.txt
#   -o, --output X.........................Report file name 
#Optional Parameters:
#   -f, --forward..........................Forward strand surffix in the fastq file, ex. R1
#   -f, --reverse..........................Reverse strand surffix in the fastq file, ex. R1
#######################################################################

import os, sys
from os import path
import subprocess
import shutil
import argparse   #take arguments
import pandas as pd
import re
import glob
import string
from datetime import datetime

coverageCutoff = 10.0  #genome coverage cut off 10%
#################################### Get sample information function ####################################
def getReadStat(inFile, strand1, strand2):
#def getReadStat(inFile):    
    """
    Read multiqc stats file qc/multiqc/raw_multiqc_data/multiqc_general_stats.txt and save to a hash
    """
    statDict0 = {}
    statDict = {}
    #get sample name, read length and total read number from file
    with open(inFile) as f:
        f.readline() #skip the title
        lines = f.readlines()
        for line in lines: 
            line = line.strip()
            cells = line.split("\t")
            smpName = cells[0]  #sample"_S1_R1_001"
            #statDict0[smpName] = {}
            if strand1 and strand1 in cells[0]:  #strand1 != ''
                smpName = re.sub(r"(.*)_%s.*" % strand1, "\\1", cells[0])  #cells[0].replace("_"+strand1, "")
                #smpName = re.sub(r"_R1_001", "", cells[0])
                statDict0[smpName] = {}
            elif strand2 in cells[0]:
                smpName = re.sub(r"(.*)_%s.*" % strand2, "\\1", cells[0])  #cells[0].replace("_"+strand2, "")
                #smpName = re.sub(r"_R2_001", "", cells[0])
                if smpName not in statDict0:
                    statDict0[smpName] = {}
            else:
                smpName = cells[0]
                statDict0[smpName] = {}
            #print(smpName)
            statDict0[smpName].setdefault("readLen",[]).append(cells[3])
            statDict0[smpName].setdefault("totalRead",[]).append(cells[5])
    #print(statDict0)
    #calculate average read length
    for smpName in statDict0:
        statDict[smpName] = []
        if len(statDict0[smpName]["readLen"]) == 1:
            statDict[smpName].append(float(statDict0[smpName]["readLen"][0]))
            statDict[smpName].append(float(statDict0[smpName]["totalRead"][0]))
        if len(statDict0[smpName]["readLen"]) == 2:  #for paired-end
            readLen = (float(statDict0[smpName]["readLen"][0]) + float(statDict0[smpName]["readLen"][1]))/2
            statDict[smpName].append(readLen)
            statDict[smpName].append(float(statDict0[smpName]["totalRead"][0]) + float(statDict0[smpName]["totalRead"][1]))        
    #print(statDict)
    return statDict

#################################### Get reference information function ####################################
def getRefInfo(refFile):
    """
    Get reference title and its length
    """
    with open(refFile) as f:
        title = f.readline().strip()
        lines = f.readlines()
        seq = ''
        for line in lines: 
            line = line.strip() 
            seq = seq + line
    seqLen = len(seq)
    title = title.replace(">", "")
    return (title, seqLen)

#################################### Get mapping coverage information function ####################################
def getCoverage(covFile):
    """
    Get mapping coverage percentage and average depth
    """
    count = count0 = total = 0
    depth = coverage = 0
    with open(covFile) as f:
        for line in f:
            line = line.strip()
            cells = line.split()
            cov = float(cells[-1])
            count += 1
            total += cov
            if cov == 0:
                count0 += 1
        if count >= 1:
            coverage = (count - count0)/float(count) * 100
            depth = total/count
            coverage = float("{:.2f}".format(coverage))
            depth = float("{:.1f}".format(depth))
    return (coverage, depth)

#################################### Get RPKM function ####################################
def getRPKM(numReads,geneLength,totalNumReads):
    """
    Calculate RPKM - Reads per kilo base per million mapped reads
    RPKM =   numReads / ( geneLength/1000 * totalNumReads/1,000,000 )
    """
    rpkm = int(numReads) / ( int(geneLength)/1000 * int(totalNumReads)/1000000 )
    mapPercent = int(numReads)/totalNumReads * 100.0
    mapPercent = float("{:.2f}".format(mapPercent))
    rpkm = float("{:.2f}".format(rpkm))
    return (rpkm,mapPercent)


#################################### Get mapping to reference information function ####################################
def getMapInfo(refFile, bamFile, coverageFile, totalNumReads):
    """
    Summary mapping information:reference length, reference title, mapped reads, coverage, rpkm, 
    """
    #refFile = 'mapping/ref/NC_Apples_3_S3/NC_028132.fasta' 
    #refFile = refDir + '/' + sample + '/' + ref + '.fasta' 
    #get reference title and lenth
    (refTitle, refLen) = getRefInfo(refFile)
    
    #get reference and its mapping information
    #bamFile = 'mapping/map2Ref/NC_Apples_3_S3/NC_028132.sorted.bam'
    #bamFile = mapDir + '/' + sample + '/' + ref + '.sorted.bam'
    cmd = 'samtools view -c -F 260 ' + bamFile
    mappedReads = os.popen(cmd).read().strip() #get sys comand return value, os.system only returns the error value.
    
    #get coverage
    #coverageFile = 'mapping/map2Ref/NC_Apples_3_S3/NC_028132.coverage.txt'
    #coverageFile = mapDir + '/' + sample + '/' + ref + '.coverage.txt'
    (coverage, depth) = getCoverage(coverageFile)
    #get RPKM
    #totalNumReads = readStat[sample][1]
    #print(mappedReads,refLen,totalNumReads)
    (rpkm,mapPercent) = getRPKM(mappedReads,refLen,totalNumReads)
    
    #mapArr = map(str, [refTitle, refLen, mappedReads,rpkm, mapPercent,coverage, depth])
    mapArr = [refTitle, refLen, mappedReads,rpkm, mapPercent,coverage, depth]
    
    return mapArr
    
#################################### get NCBI blast result function ####################################
def getBlast(blastFile):
    """
    Get NCBI blast result, it is the second line of the blastFile
    """
    #T45588B2e_S8.NC_035759.consensus.blastn.txt
    #blastFile = ncbiBlastDir + '/' + sample + '.' + ref + '.consensus.blastn.txt'/'.contig.blastn.txt'/'contig.blastx.txt'
    cmd = 'head -2 ' + blastFile + ' | tail -1'
    blast = os.popen(cmd).read().strip()

    #print(blastFile)
    #print(blast)
    
    return blast

#################################### get taxonomomy function ####################################
def getTax(annDir):
    """
    Get a reference taxonomomy from sample.selectedRef.txt in annotation dir
    """
    taxonDict = {}
    #get all sample.selectedRef.txt
    refFiles = [f for f in glob.glob(annDir + "/*.selectedRef.txt", recursive=True)]
    for file in refFiles:
        #print(file)
        with open(file) as fin:
            for line in fin.readlines()[1:]:
                cells = list(map(str.strip, line.split("\t"))) #remove space at the start and end
                id = cells[1]  #reference name, NC_003844.1
                if "acc|" in id:  #acc|GENBANK|QGY72600.1|GENBANK|MN532657|RNA 
                    id = id.split("|")[-2]
                if "." in id:
                    id = id.split(".")[0]
                taxon = cells[13:16] #[taxon_id, Species, taxon path]

                taxonDict[id] = taxon
    
    return taxonDict

#################################### get Acronym function ####################################
def getAcronym(acronymFile):
    """
    Get a acronym for a taxon id
    """
    acronymDict = {}
    with open(acronymFile) as f:
        for line in f:
            line = line.strip()
            cells = line.split("\t")
            if len(cells) > 1:
              acronymDict[cells[0]] = cells[1] 

    return acronymDict    

    
#################################### Put all sample and mapping information together function ####################################
def setReport(reportDict, refNameFile, refDir, mapDir, ncbiBlastDir, sample, trimReadStat, taxonDict, acronymDict, blastSuffix):
    """
    Put all (raw read stat, trimmed read stat, mapping infor, blast) together
    """

    refName = [line.strip() for line in open(refNameFile, 'r')]
    #print(refName)
    #print(sample)
    #smp = re.sub(r'_S\d+$', '',   sample)
    #print(readStat[sample])
    #totalNumReads = readStat[sample][1]
    cleanNumReads = trimReadStat[sample][1]
    #readStat0 = reportDict.pop(sample + "|nan")
    readStat0 = reportDict[sample + "|nan"]
    #print(readStat0)
    for ref in refName:
        refFile = refDir + '/' + sample + '/' + ref + '.fasta'
        bamFile = mapDir + '/' + sample + '/' + ref + '.sorted.bam'
        coverageFile = mapDir + '/' + sample + '/' + ref + '.coverage.txt'
        
        #get mapping to reference information
        #mapArr = getMapInfo(refFile, bamFile, coverageFile, totalNumReads)
        mapArr = getMapInfo(refFile, bamFile, coverageFile, cleanNumReads)

        #filter low coverage (<10% genome) virus
        coverage = mapArr[5]
        #print(list(mapArr))
        if "consensus.blastn" in blastSuffix and coverage < coverageCutoff: 
        #   #print(refFile, coverage)
            print("Low genome coverage virus filtered: ", ref, " in ", sample)
            continue
        
        #get concensus ncbi Blast information
        blastnFile = ncbiBlastDir + '/' + sample + '.' + ref + blastSuffix #'.consensus.blastn.txt', '.contig.blastx.txt'
        blastxFile = ncbiBlastDir + '/' + sample + '.' + ref + '.contig.blastx.txt' #for config blastx result, it's possible a novel virus
        if os.path.isfile(blastnFile):
            blast = getBlast(blastnFile)
        elif os.path.isfile(blastxFile):
            blastSuffix = '.contig.blastx.txt'
            blast = getBlast(blastxFile)
        else:
            blast = "nan\tnan\tnan\tnan\tnan"

        id = sample + "|" + ref
        #in case reference id in both consensus and contig
        if id in reportDict:
            print(id, " has both consensus and contig in the report")
            id = sample + "|" + ref + ".1"

        reportDict.setdefault(id, [])
        #reportDict[id].extend(map(str,[round(readStat[sample][0],2), int(readStat[sample][1]), round(trimReadStat[sample][0],2), int(trimReadStat[sample][1])]))
        reportDict[id].extend(readStat0)
        #print(id, ref)
        #print(taxonDict)
        #print(taxonDict[ref])
        reportDict[id].append(ref)
        try: 
            #print(taxonDict[ref][1], acronymDic[taxonDict[ref][1]])
            reportDict[id].extend(taxonDict[ref])
            reportDict[id].append(acronymDict[taxonDict[ref][1]]) #species vs. acronym
        except KeyError:
            #print(ref, "ref is missed")
            reportDict[id].append("nan")
        #reportDict[id].extend(mapArr)
        reportDict[id].extend(map(str, mapArr))

        blastType = blastSuffix.split(".")[-3] + "." + blastSuffix.split(".")[-2]
        reportDict[id].append(blastType)
        reportDict[id].extend(blast.split("\t")[1:]) #skip coverage
    
#################################### Get arguments function ####################################
def parseArguments():
    """
    Get arguments
    """
    parser = argparse.ArgumentParser(description='Generate text report')
    parser.add_argument("-w", "--workdir", dest='workdir', required=True, help="Working directory, ex. /ppq/data2/pgqp_pipeline/Run45")
    parser.add_argument('-f','--forward',dest='strand1', required=False, default="R1", help='Forward strand, ex. R1')
    parser.add_argument('-r','--reverse',dest='strand2', required=False, default="R2", help='Reverse strand, ex. R2')
    parser.add_argument("-a", "--acronym", dest='acronym_file', required=True, help="Acronym file, ex. ICTV_virus_acronym2019.txt")
    parser.add_argument("-o", "--output", dest='output_file', required=True, help="Report file name, ex. Run45.report.txt")
    parser.set_defaults(append=False)
    return parser.parse_args()


def main():
    ### Input arguments
    options = parseArguments()
    workdir = options.workdir
    strand1 = options.strand1
    strand2 = options.strand2
    acronymFile = options.acronym_file #'/ppq/data2/resources/ncbi/ICTV_virus_acronym2019.txt'  
    outFile = options.output_file
    annDir = workdir + '/annotation'
    ncbiBlastDir = workdir + '/report/ncbiBlast'
    #blastFile = workdir + '/report/consensus.blast.txt'
    readStatFile = workdir + '/qc/multiqc/raw_multiqc_data/multiqc_general_stats.txt'
    trimReadStatFile = workdir + '/qc/multiqc/trimmed_multiqc_data/multiqc_general_stats.txt'
    refDir = workdir + '/mapping/ref'
    mapDir = workdir + '/mapping/map2Ref'

    #get total reads
    readStat = getReadStat(readStatFile, strand1, strand2)
    #get clean clean reads, afer removing duplicate reads and PhiX174 contaminant, triming low quality bases
    trimReadStat = getReadStat(trimReadStatFile, strand1, strand2)

    #print(readStat)
    #print(trimReadStat)
    #sync the sample names in the dicts: readStat, trimReadStat
    for id1 in trimReadStat: #Q51550AB_S1
        #print("trim id: ", id1)
        for id2 in list(readStat): #Q51550AB_S1_R1_001, '_S1_R1_001' is added by bcl2fastq
            #print("raw id: ", id2)
            if id1 in id2:
                readStat[id1] = readStat.pop(id2) #change dict key

    #get reference taxonomomy
    taxonDict = getTax(annDir) #ref: [taxonid, species, taxonomy]
    #print(taxonDict)
    
    #get acronym 
    acronymDict = getAcronym(acronymFile) # taxonid: acronym
        
    #get sample names
    #samples = [os.path.basename(x) for x in glob.glob(annDir + "/*.blastn.txt", recursive=True)]
    #samples = [x.replace(".blastn.txt","") for x in samples]
    samples = list(trimReadStat.keys())
    #print(samples)

    #put all samples report tegother
    reportDict = {}
    #add readStat and trimReadStat
    for sample in samples:
        id = sample+"|nan"
        reportDict.setdefault(id, []).append(sample) #remove _S1 from sample name
        reportDict[id].extend(map(str,[round(readStat[sample][0],2), int(readStat[sample][1]), round(trimReadStat[sample][0],2), int(trimReadStat[sample][1])]))
    #print(reportDict)

    #add mapping and blast information
    for sample in samples:
        #print(sample)
        refNamenFile = annDir + '/' + sample + '.refNamen.txt'  #reference from blastn ex. 'annotation/NC_Apples_3_S3.refNamen.txt'
        refNamexFile = annDir + '/' + sample + '.refNamex.txt' #reference from blastx ex. 'annotation/NC_Apples_3_S3.refNamex.txt'
        refnSize = refxSize = 0
        if path.exists(refNamenFile):
            refnSize =  os.stat(refNamenFile).st_size #check whether file is empty, os.stat("file").st_size == 0
        if path.exists(refNamexFile):
            refxSize = os.stat(refNamexFile).st_size
            #print(type(refxSize))
        
        #blastn results
        if refnSize > 0: 
            blastSuffix = '.consensus.blastn.txt'
            setReport(reportDict, refNamenFile, refDir, mapDir,ncbiBlastDir, sample, trimReadStat, taxonDict, acronymDict, blastSuffix)
            print("Blastn: ", sample)
        #blastx results
        if refxSize > 0: 
            #get reference names
            #refDir = 'mapping/finalAssembly'
            #mapDir = 'mapping/map2Contig'
            print("Blastx: ", sample)
            blastSuffix = '.contig.blastn.txt'
            setReport(reportDict, refNamexFile, refDir, mapDir,ncbiBlastDir, sample, trimReadStat, taxonDict, acronymDict, blastSuffix)
        if refnSize == 0 and refxSize == 0:
            id = sample+"|nan"
            #reportDict.setdefault(id, []).append(sample) #remove _S1 from sample name
            #reportDict[id].extend(map(str,[round(readStat[sample][0],2), int(readStat[sample][1]), round(trimReadStat[sample][0],2), int(trimReadStat[sample][1])]))
            #reportDict[id].extend(["nan","nan","nan","nan","nan","nan","nan","nan","nan","nan","nan","nan"])
            #reportDict[id].extend(["nan"]*18)
            reportDict[id].extend(["nan"]*17)
            print("Sample no virus: ", sample)
    
    #Filter noise, just keep virus, viroid, satellite
    for id in list(reportDict):
        if len(reportDict[id][-2]) < 22:
            del reportDict[id]  #remove dict item only having readstat
        elif len(reportDict[id][-2]) < 6 or not ("virus" in reportDict[id][-2] or "viroid" in reportDict[id][-2] or "satellite" in reportDict[id][-2] or "viridae" in reportDict[id][-2] or "virales" in reportDict[id][-2] or "nan" in reportDict[id][-3]):
            print("Filter non viral annotation: ", id, ",", reportDict[id][-2])
            del reportDict[id]

    #check missing samples 
    for sample in samples:
        smpExist = False
        for id in reportDict:
            if sample in id:
                smpExist = True
                break
        if not smpExist: #if sample is filtered, add it back
            id = sample+"|nan"
            reportDict.setdefault(id, []).append(sample) #remove _S1 from sample name
            reportDict[id].extend(map(str,[round(readStat[sample][0],2), int(readStat[sample][1]), round(trimReadStat[sample][0],2), int(trimReadStat[sample][1])]))
            #reportDict[id].extend(["nan"]*18)
            reportDict[id].extend(["nan"]*17)

    #output
    #print(reportDict.values())
    #output file
    fout = open(outFile, 'w')
    header = "Sample\tReadLength\tTotalReads\tCleanReadLength\tCleanReads\t"
    header = header + "RefId\tTaxonId\tSpecies\tTaxonomy\tAcronym\tRefTitle\tRefLen\tMappedReads\tRPKM\tPercentMappedreads\tPercentGenomeCovered\tMeanCoverage\tBlastType\tBlastEvalue\tBlastIdentity\tBlastDescription\tSequence"
    fout.write(header + "\n")
    for id in reportDict:
        fout.write("\t".join(reportDict[id]) + "\n")
    
    fout.close()

if __name__ == "__main__":
    main()

    