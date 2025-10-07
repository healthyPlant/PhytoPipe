#! /usr/bin/env python
######################################################################
#Alex Hu <xiaojun.hu@usda.gov>
#Updated: 04/30/2022
#This program filter report.txt
#Usage: python filterReport.py [report file] [filterKey file] [output file]
######################################################################

import os
import sys
import re

if len(sys.argv) < 4:
    print("Usage: python filterReport.py [report file] [filterKey file] [output file]")
    print("Argument missed, exit!")
    sys.exit()
rptFile = sys.argv[1] #workdir + '/report/report.txt'
keyFile = sys.argv[2] #filterKeys.txt
outFile = sys.argv[3] #workdir + '/report/report.filtered.txt''

#count records for each sample
countDict = {}
recordDict = {}
sampleInfor = {}
deleted = {}
filtered = {}

print("Filter viral results")
#read report 
#Sample	RawReadLength	TotalReads	CleanReadLength	CleanReads	RefId	TaxonId	Species	Taxonomy	Acronym	RefTitle	RefLen	MappedReads	RPKM	PercentMappedreads	PercentGenomeCovered	MeanCoverage	BlastType	BlastQueryCover	BlastEvalue	BlastIdentity	BlastDescription	Sequence
with open(rptFile) as f:
    header = f.readline()
    for line in f:
        line = line.rstrip()
        cells = line.split("\t")
        if cells[0] in countDict:
            countDict[cells[0]] += 1
        else:
            countDict[cells[0]] = 1
            
        id = cells[0]+"|"+cells[6]+"|"+cells[5]  #Sample + TaxonID + RefId 
        recordDict[id] = cells
        #sampleInfor[cells[0]] = cells[0:15]
print("start records: ", len(recordDict))
#print(recordDict)

#get filter keys
filterKeys = []
with open(keyFile) as f:
    for line in f:
        if line.startswith('#'):
            continue
        line = line.rstrip()
        filterKeys.append(line)
#print(filterKeys)

#filter 1: filter by key words: fungi, insect, plant host, and other unwanted
print("***************************************************************")
print("Filter enviromental viruses or host sequences")
for record in list(recordDict): #using list to prevent dict size change
    #print(recordDict[record])
    for fkey in filterKeys:
        #check RefID, Taxonomy, RefTitle, or BlastDescription
        #print(recordDict[record][5], recordDict[record][8], recordDict[record][10], recordDict[record][-2])
        try:
            if fkey in recordDict[record][5] or fkey in recordDict[record][8] or fkey in recordDict[record][10] or fkey in recordDict[record][-2]:
                print("Filter: ", recordDict[record][0], ",", recordDict[record][-2])
                countDict[recordDict[record][0]] -= 1 #keep counting
                #If all records in a sample are noise, change records to nan
                if countDict[recordDict[record][0]] == 0:
                    id = recordDict[record][0] + "|nan|nan"
                    recordDict[id] = recordDict[record][0:5] + ["nan"]*18 #extend a list
                    #print(recordDict[id])
                del recordDict[record]  #deleted unwanted records
        except KeyError:
            #print(record, " is deleted.")
            #continue
            deleted[record] = 1
            
#print(recordDict)
#print("After keyword filter records: ", len(recordDict))
print(f"Not plant viruses records: {len(deleted)}")
print("***************************************************************")
print("Filter repeated viruses, which may be viral quasispecies")

#filter 2: filter repeated taxons in a sample. They may be quansispecies
count = 0
for record in recordDict:
    #print(record)
    #skip multiple segment viruses
    if 'segment' in recordDict[record][10]: #segment in RefTitle
        filtered[record] = recordDict[record]
    elif 'nan' in record:
        filtered[record] = recordDict[record]
    else:
        #if one sample has many same taxon Id viruses, keep the hihger MappedReads and PercentGenomeCovered one 
        id = recordDict[record][0] + "|" + recordDict[record][6]  #Sample + TaxonId
        #print(id, recordDict[record])
        if id in filtered:
            print("Filter: ", filtered[id][0], ",", filtered[id][-2])
            if int(filtered[id][12]) < int(recordDict[record][12]): #compare MappedReads
                filtered[id] = recordDict[record]
            elif filtered[id][12] == recordDict[record][12]:
                if float(filtered[id][15]) < float(recordDict[record][15]): #compare PercentGenomeCovered
                    filtered[id] = recordDict[record]
            count += 1
        else:
            filtered[id] = recordDict[record]
                
#print(filtered)
print("***************************************************************")
print("Repeat taxon records: ", count)
print("Final records: ", len(filtered))


#get order of the samples: S1, S2, ....
sampleOrder = {}
sampleNums = set()
orderSample = False
for id in filtered:
    seqId = id.split("|")[0]
    if re.search(r"_S(\d+)$", seqId):
        orderSample = True
        sampleNum = int(re.search(r"_S(\d+)$", seqId).group(1))
        #print(sampleNum.group(1), seqId)
        sampleNums.add(sampleNum)
        sampleOrder.setdefault(sampleNum, []).append(id) 
        #if sampleNum >= sampleSize:  #get the max sample number (sample range)
        #    sampleSize = sampleNum

#sampleSize = len(sampleNums)

#print(sampleOrder)
#print(sampleNums)
#output to file
fout = open(outFile, 'w')
fout.write(header)
if orderSample:
    for n in sorted(sampleNums): #range(1, sampleSize+1): 
        for id in sampleOrder[n]:
            #print(sampleOrder[n])
            #print(filtered[id][0:20])
            fout.write("\t".join(filtered[id]) + "\n")
else:
    for id in filtered:
        fout.write("\t".join(filtered[id]) + "\n")

fout.close()  

