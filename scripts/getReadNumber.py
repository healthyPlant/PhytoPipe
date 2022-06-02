#! /usr/bin/env python
######################################################################
#Alex Hu <xiaojun.hu@usda.gov>
#Updated: 02/18/2021
#This program get read numbers at each QC step: 
#1. raw reads QC: read#, yield, >=Q30Bases%, mean quality
#2. host ribosomal RNA %
#3. read# after remove duplicates
#4. control read#
#5. read# after trimming
#6. read# after removing host

#Required Parameters:
#   -w, --workdir X........................working directory
#   -m, --format X.........................read file format, ex. fastq.gz
#   -o, --output X.........................Report file name 
#Optional Parameters:
#   -f, --forward..........................Forward strand surffix in the fastq file, ex. R1
#   -f, --reverse..........................Reverse strand surffix in the fastq file, ex. R1
#######################################################################

import sys, argparse
import glob
import subprocess
import os
import json
import re
import multiprocessing 
import gzip

#################################### Get raw read information function ####################################  
def parseJson(jsonFile):
    """
    Get raw quality information from a jason file which generate by Illumina bcl2fastq
    """
    with open(jsonFile) as f:
        data = json.load(f) #json.load() output a dict
    #print(data["ConversionResults"]) # a list
    qualDict = {} #a dict of lists to store sample read quality
    for cr in data["ConversionResults"]: # a list of dict
        #print(cr.keys()) #'LaneNumber', 'TotalClustersRaw', 'TotalClustersPF', 'Yield', 'DemuxResults', 'Undetermined'
        for sample in cr['DemuxResults']: # a list of dict
            #print(sample.keys()) #'SampleId', 'SampleName', 'IndexMetrics', 'NumberReads', 'Yield', 'ReadMetrics'
            #print(sample['SampleId'], sample['NumberReads'], sample['Yield'], end=' ')
            
            yieldQ30 = qualityScoreSum = 0
            for read in sample['ReadMetrics']: # a list of dict
                #print(read["Yield"], read['YieldQ30'], read['QualityScoreSum'])
                yieldQ30 += read['YieldQ30']
                qualityScoreSum += read['QualityScoreSum']
            avgQual = float(qualityScoreSum/sample['Yield'])
            q30Percent = float(yieldQ30/sample['Yield'] * 100)
            #print(yieldQ30,qualityScoreSum, "{:.2f}".format(q30Percent), "{:.2f}".format(avgQual))    
            #qualDict.setdefault(sample['SampleId'],[]).append([sample['NumberReads'], sample['Yield'],q30Percent, avgQual])
            qualDict.setdefault(sample['SampleName'],[]).append([sample['NumberReads'], sample['Yield'],q30Percent, avgQual])
                
    #there are several lanes , get total yield and average qualily for a sample
    avgQualDict = {}
    for sample in qualDict:
        #print(sample, end=' ')
        #print(qualDict[sample])
        laneNum = len(qualDict[sample])
        totalRead = totalYield = totalQ30 = totalQual = 0
        for read in qualDict[sample]:
            totalRead += read[0]
            totalYield += read[1]
            totalQ30 += read[2]
            totalQual += read[3]
        avgQ30 = totalQ30/laneNum
        avgQual = totalQual/laneNum
        #print(totalYield, avgQ30, avgQual)
        avgQualDict[sample] = [totalRead, "{:.0f}".format(totalYield/1000000), "{:.2f}".format(avgQ30), "{:.2f}".format(avgQual)]
    return(avgQualDict)

#################################### Get a raw fastq information function ####################################  
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
        #if lineno%4 == 1: flag = (line in ids)
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

#################################### Get all raw fastq information function ####################################  
def sumQual(fastqFolder, input_format, strand1, strand2):
    """
    Run parallel to get each fastq file quality
    """
    #fastqFiles = [os.path.basename(x) for x in glob.glob(fastqFolder + '/*.' +  input_format)]
    fastqFiles = [x for x in glob.glob(fastqFolder + '/*.' +  input_format)]
    #Parallel Processing multiple fastq.gz file 
    pool = multiprocessing.Pool() 
    outputs_async = pool.map_async(getQualFromFastq, fastqFiles) 
    outputs = outputs_async.get() #return a list of sets
    #print("Output: {}".format(outputs)) # a list

    #convert set to list
    for i in range(len(outputs)):
        outputs[i] = list(outputs[i])
    #remove _001 in sample name    
    for output in outputs:
        if output[0].endswith("_001"):
            output[0] = output[0].replace("_001", "") 
               
    #check paired-end or single-end read
    paired = False
    for output in outputs:
        if output[0].endswith(strand2):
            paired = True
            break

    #change a list to a dict
    avgQualDict = {}    
    #for paired-end reads
    avgQualDictR1 = {}
    avgQualDictR2 = {}

    for output in outputs:
        if output[0].endswith(strand1):
            sname = output[0].replace("_"+strand1, "")
            avgQualDictR1[sname] = output[1:len(output)]
        elif output[0].endswith(strand2): #paired-end read
            sname = output[0].replace("_"+strand2, "")
            avgQualDictR2[sname] = output[1:len(output)]
        else: #single-end read
            avgQualDict[output[0]] = output[1:len(output)]

    if paired: #combine R1 and R2
        for sname in avgQualDictR2:
            avgQualDict[sname] = (avgQualDictR1[sname][0] + "(F)|" + avgQualDictR2[sname][0] + "(R)",  avgQualDictR1[sname][1] + "(F)|" + avgQualDictR2[sname][1] + "(R)", avgQualDictR1[sname][2] + "(F)|" + avgQualDictR2[sname][2] + "(R)", avgQualDictR1[sname][3] + "(F)|" + avgQualDictR2[sname][3] + "(R)" )
    elif not avgQualDict:
        avgQualDict = avgQualDictR1

    return(avgQualDict)

#################################### Get arguments function ####################################
def parseArguments():
    """
    Get arguments
    """
    parser = argparse.ArgumentParser(description='Generate text report')
    parser.add_argument("-w", "--workdir", dest='workdir', required=True, help="Working directory, ex. /ppq/data2/pgqp_pipeline/Run45")
    parser.add_argument('-f','--forward',dest='strand1', required=False, default="R1", help='Forward strand, ex. R1')
    parser.add_argument('-r','--reverse',dest='strand2', required=False, default="R2", help='Reverse strand, ex. R2')
    parser.add_argument("-m", "--format", dest='format', required=True, help="read file format, ex. fastq.gz")
    parser.add_argument("-o", "--output", dest='output_file', required=True, help="output file name, ex. qcReadNumber.txt")
    parser.set_defaults(append=False)
    return parser.parse_args()


def main():
    ### Input arguments
    options = parseArguments()
    workdir = options.workdir  #'/ppq/data2/pgqp_pipeline/Run45'
    strand1 = options.strand1  #R1
    strand2 = options.strand2  #R2
    input_format = options.format #"fastq.gz"
    outFile = options.output_file #workdir + '/report/qcReadNumber.txt'

    #QC folders
    removeLogDir = workdir + '/logs/removeDuplicate'
    contaminantLogDir = workdir + '/logs/removeControl'
    trimLogDir = workdir + '/logs/trimmomatic'
    cleanDir = workdir + '/cleaned'
    rRNADir = workdir + '/logs/rRNA_qc'

    statJsonFile = workdir + '/raw/Stats/Stats.json'
    rawDir = workdir + '/raw'

    samples = [os.path.basename(x) for x in glob.glob(trimLogDir + '/*.log')]
    samples = [x.replace(".log","") for x in samples]
    #print(samples)


    orderSmp = [" " for x in samples]
    if os.path.exists(statJsonFile):
        #sort sample by S\d number
        for sample in samples:
            #print(sample)
            m = re.match("\S+_S(\d+)$",sample)
            if m:
                #print(m.groups())
                sNum = int(m.groups()[0]) - 1
                orderSmp[sNum] = sample
        #print(orderSmp)
    else:
        orderSmp = samples

    header="Sample\tRawReads\tRawYield(Mbases)\tPercent>=Q30Bases\tRawMeanQualityScore\tRibosomalRNA\tReadsAfterRemoveDuplicates\tControlReads\tReadsAfterTrim\tPossiblePathogenReads"

    fout = open(outFile, 'w')
    fout.write(header + "\n")

    qualDict = {}
    if os.path.exists(statJsonFile):
        qualDict = parseJson(statJsonFile)
    #elif(os.path.isdir(rawDir)):
    else:
        qualDict = sumQual(rawDir, input_format, strand1, strand2)
        
    #print(qualDict)
    for sample in orderSmp:
        #print(sample)
        rmFile = removeLogDir + "/" + sample + ".log"
        ctFile = contaminantLogDir + "/" + sample + ".log"
        tmFile = trimLogDir + "/" + sample + ".log"
        clFile = cleanDir + "/" + sample + ".pathogen.fastq.gz"
        rnaFile = rRNADir + "/" + sample + ".log"
        #print(rnaFile)
        paired = False
        if not os.path.exists(clFile):
            clFile = cleanDir + "/" + sample + "_" + strand2 + ".pathogen.fastq.gz"
            paired = True
        
        #get rRNA %
        #bbduk reports total reads (2 * read pairs for paired-end reads)
        cmd = "grep 'Total Removed:' " + rnaFile + " | awk '{gsub(/[()]/,\"\");print $5}' "
        rrna = subprocess.check_output(cmd, shell=True)
        rrna = rrna.decode("utf-8").strip()
        
        cmd = "grep -m 1 -w 'Input:' " + rnaFile + " | awk '{print $2}'" #-m1 for first match
        raw = subprocess.check_output(cmd, shell=True)
        #print(sample, raw)
        raw = raw.decode("utf-8").strip() # using decode() function to convert byte string to string

        #clumpy.sh reports read number (not read pairs for paired-end)
        cmd = "tail -10 " + rmFile + " | grep -w 'Reads Out:' | awk '{print $3}'"
        #print(cmd)
        #clumpify.sh report all reads, paired is double pairs
        dup = subprocess.check_output(cmd, shell=True)
        #print(sample, dup)
        dup = dup.decode("utf-8").strip()
        #in case no dup value, assign 0
        if not dup:
          dup = str(0)
        
        #bbsplit reports read1 and read2 for paired-end        
        cmd = "grep -w 'mapped:' " + ctFile + " | awk '{print $3}'"
        ctm = subprocess.check_output(cmd, shell=True) 
        ctm = ctm.decode("utf-8").strip()
        if paired:
            ctm = ctm.replace("\n", "(F)|")
            ctm = ctm + "(R)"
        #print(ctm)
        
        #Trimmomatic report Read Pairs for paired-end
        cmd = "grep 'Input Read' " + tmFile + " | awk -F: '{print $3}' | cut -d' ' -f2"
        trim = subprocess.check_output(cmd, shell=True)
        trim = trim.decode("utf-8").strip()
        #prsint(trim)
        if paired:
            trim = str(int(trim) * 2)

        #count read number from fastq.gz file
        #for paired-end, 2*read#
        cmd = "zcat " + clFile + " | wc -l"
        clean = subprocess.check_output(cmd, shell=True)
        clean = int(clean.decode("utf-8").strip())/4
        if paired:
            clean = str(int(clean) * 2)

        #raw = qualDict[smp][0]
        #print("Dup:", dup, "Raw:", raw)
        percentDup = "{:.2f}".format(float(dup)/float(raw)*100)
        percentTrim = "{:.2f}".format(float(trim)/float(raw)*100)
        percentClean = "{:.2f}".format(float(clean)/float(raw)*100)
        #print(raw,dup,percentDup,ctm,trim,percentTrim, clean, percentClean)
        #     

        if os.path.exists(statJsonFile):
            #remove _S\d at the end of sample
            smp = re.sub(r'_S\d+$', '', sample)
            fout.write(sample + "\t")
            try:
                if paired and type(qualDict[smp][0]) == int:
                    qualDict[smp][0] = int(qualDict[smp][0])*2
                rawQual = '\t'.join(map(str, qualDict[smp]))
                #print(rawQual)

                fout.write(rawQual) #sample + "\t" + rawQual + 
            except KeyError:
                print(smp, " missed!") 
                fout.write("NA\tNA\tNA\tNA")
        else:
            fout.write(sample + "\t")
            fout.write('\t'.join(qualDict[sample]))
        fout.write("\t" + rrna + "\t" + dup + " (" + percentDup + "%)" + "\t" + ctm + "\t" + str(trim) + " (" + percentTrim + "%)" + "\t" + str(int(clean)) + " (" + percentClean + "%)" + "\n") #  
    fout.close()  

if __name__ == "__main__":
    main()

