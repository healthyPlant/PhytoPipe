#! /usr/bin/env python
######################################################################
#Alex Hu <xiaojun.hu@usda.gov>
#Updated: 02/18/2021
#This program generate a html report, which includes FastQC, multiQC, Kraken2, Kaiju, blastn, blastx
#Usage: python getHtmlReport.py [working directory] [monitor pathogen file] [report] [output file]
######################################################################

import pandas as pd
import os
import glob
import shutil
import sys


if len(sys.argv) < 6:
    print("Usage: python getHtmlReport.py [working directory] [monitor pathogen file] [report folder] [report file] [output file]")
    print("Argument missed, exit!")
    sys.exit()
workdir = sys.argv[1]
mpFile = sys.argv[2] #monitor pathogens
rptdir = sys.argv[3] #workdir + '/report'
rptFile = sys.argv[4] #workdir + '/report/report.txt'
outFile = sys.argv[5] #workdir + '/report/report.html'

#rptFile = rptdir + '/report.txt'

#blast threshold to display
blastnThreshold = 1e-100 #evalue
blastxThreshold = 1e-200 #evalue

pd.set_option('display.max_colwidth',100)

def getMonitorPathogen(mpFile):
    mpDict = {}
    with open(mpFile) as f:
        for line in f:
            line = line.rstrip()
            cells = line.split("\t")
            mpDict[cells[0]] = cells[1]
    return mpDict

def extractKraken(kFile, mpFile):
    """
    Extract unclassified (0), Eukaryota (2759), Bacteria (2), Viruses(10239), Fungi (4751), Oomycetes(4762)
    Output species if Eukaryota (2759), Bacteria (2), Fungi (4751), Oomycetes(4762) Percentage of fragments > 1%, Viruses(10239) Number of fragments > 100
    keep percentage, mapped reads, Number of reads assigned directly to this taxon, species name
    """
    #20.06  3234437 3234437 U       0       unclassified
    viral_cutoff = 100  #mapped viral reads cut off
    kDict = {}
    mpDict = {}
    if os.path.exists(mpFile):
        mpDict = getMonitorPathogen(mpFile) #get monitor pathogens

    bcount = fcount = ccount = 5 #limit display top 5 species for Bacteria, Fungi and Oomycetes

    #initial kDict
    for key in ['Eukaryota', 'Bacteria', 'Viruses', 'Fungi', 'Oomycetes']:
        kDict[key] = []

    kraken = []
    eFlag = bFlag = vFlag = fFlag = oFlag = 0
    ecount = 0
    with open(kFile) as f:
        for line in f:
            line = line.rstrip()
            #cells = line.split("\t")
            cells = list(map(str.strip, line.split("\t")))
            if cells[5] == 'uncultured organism':
                continue

            if eFlag != 1 and (cells[3]=='D' or cells[3]=='K'):
                eFlag = bFlag = fFlag = oFlag = 0 #reset marker

            if cells[4] == '2759': #Eukaryota
                eFlag = 1
                bFlag = vFlag = fFlag = oFlag = 0
                continue
            if eFlag == 1 and cells[3].startswith('S'): 
                ecount += 1
                if ecount == 1: #add first Eukaryota species, it's host
                    kDict['Eukaryota'].append([cells[0],cells[2],cells[4],cells[5],'Eukaryota'])
                elif float(cells[0]) > 5.0: #keep >5% Eukaryota species (host)
                    kDict['Eukaryota'].append([cells[0],cells[2],cells[4],cells[5],'Eukaryota'])

            if cells[4] == '4751':
                fFlag = 1
                eFlag = bFlag = vFlag = oFlag = 0
                continue
            if fFlag == 1 and cells[3].startswith('S') and 'uncultured' not in cells[5]: 
                if (cells[4] in mpDict and float(cells[2]) > 100) or (float(cells[2]) > 1000 and fcount > 0): #keep read >1000 species
                    kDict['Fungi'].append([cells[0],cells[2],cells[4],cells[5],'Fungi'])
                    fcount -= 1
                
            if cells[4] == '2':
                bFlag = 1
                eFlag = vFlag = fFlag = oFlag = 0
                continue
            if bFlag == 1 and cells[3].startswith('S') and 'uncultured' not in cells[5]:
                if float(cells[2]) > 1000 and bcount > 0: #keep read >1000 species
                    kDict['Bacteria'].append([cells[0],cells[2],cells[4],cells[5],'Bacteria'])
                    bcount -= 1
                elif 'Candidatus Phytoplasma' in cells[5] and float(cells[2]) > 20: #keep read >20 for Phytoplasma
                    #report Candidatus Phytoplasma
                    kDict['Bacteria'].append([cells[0],cells[2],cells[4],cells[5],'Bacteria'])
                elif cells[4] in mpDict and float(cells[2]) > 100:
                    kDict['Bacteria'].append([cells[0],cells[2],cells[4],cells[5],'Bacteria'])
                
            if cells[4] == '4762':
                oFlag = 1
                eFlag = bFlag = vFlag = fFlag = 0
                continue    
            if oFlag == 1 and cells[3].startswith('S'): 
                if (cells[4] in mpDict and float(cells[2]) > 100) or (float(cells[2]) > 1000 and ccount > 0): #keep read >1000 species
                    kDict['Oomycetes'].append([cells[0],cells[2],cells[4],cells[5],'Oomycetes'])
                    ccount -= 1
                
            if cells[4] == '10239':
                vFlag = 1
                eFlag = bFlag = fFlag = oFlag = 0
                continue    
            if vFlag == 1 and cells[3].startswith('S') and 'synthetic construct' not in cells[5]: #keep read >100 species
                if (cells[4] in mpDict and float(cells[2]) > 10) or float(cells[2]) > viral_cutoff: 
                    kDict['Viruses'].append([cells[0],cells[2],cells[4], cells[5],'Viruses'])
            
                
    for key in ['Eukaryota', 'Bacteria', 'Viruses', 'Fungi', 'Oomycetes']:
        if kDict[key]:
            kraken.extend(kDict[key])
    return kraken

def extractKaiju(kjFile, mpFile):
    """
    Extract Bacteria (2), Viruses(10239), Fungi (4751), Oomycetes(4762) from Kaiju table result
    Output species if Bacteria (2), Fungi (4751), Oomycetes(4762) read number > 1000, Viruses(10239) Number of fragments > 100
    keep percentage, mapped reads, Number of reads assigned directly to this taxon, species name
    """
    #file    percent reads   taxon_id        taxon_name
    #classification/SRX6762507.kaiju.out     2.099517        515840  470     cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Moraxellaceae;Acinetobacter;Acinetobacter calcoaceticus/baumannii complex;Acinetobacter baumannii;
    viral_cutoff = 100  #mapped viral reads cut off
    kjReport = []
    mpDict = {}
    bcount = fcount = ccount = 5 #limit display top 5 species for Bacteria, Fungi and Oomycetes

    if os.path.exists(mpFile):
        mpDict = getMonitorPathogen(mpFile) #get monitor pathogens
        
    with open(kjFile) as f:
        f.readline() #skip header
        for line in f:  
            line = line.rstrip()
            if "environmental samples" in line:
                continue

            cells = line.split("\t")
            domain = 'None'
            species = 'None'
            taxon = cells[4].split(";")

            if cells[4].startswith("cellular organisms;Eukaryota;Sar;Stramenopiles;Oomycota;"):
                if (float(cells[2]) >= 1000 and ccount > 0) or cells[3] in mpDict:
                    if len(taxon) >= 8:
                        domain = 'Oomycetes'
                        species = taxon[-2]
                        ccount -= 1
            elif cells[4].startswith("cellular organisms;Eukaryota;Opisthokonta;Fungi;"): 
                if (float(cells[2]) >= 1000 and fcount > 0) or cells[3] in mpDict:
                    if len(taxon) >= 8:
                        domain = 'Fungi'
                        species = taxon[-2]
                        fcount -= 1
            elif cells[4].startswith("cellular organisms;Bacteria;"): 
                if "Candidatus Phytoplasma" in cells[4] and float(cells[2]) >= 10:
                    domain = 'Bacteria'
                    species = taxon[-2]
                elif (float(cells[2]) >= 1000 and bcount > 0) or cells[3] in mpDict:
                    if len(taxon) >= 5:
                        domain = 'Bacteria'
                        species = taxon[-2]
                        bcount -= 1
            elif cells[4].startswith("Viruses;") and (float(cells[2]) >= viral_cutoff or (cells[3] in mpDict and float(cells[2]) >= 10) ):
                domain = 'Viruses'
                if len(taxon) >= 5 and " " in taxon[-2] :
                    species = taxon[-2]
            else:
                continue
            if species != 'None':
                #kjReport.append([cells[1],cells[2],cells[3],species, domain, cells[4]])
                kjReport.append([cells[1],cells[2],cells[3],species, domain])
    return kjReport

# First table: Use Pandas <a href="http://pandas.pydata.org/pandas-docs/version/0.15.1/generated/pandas.DataFrame.to_html.html">to_html()</a> routine to convert Pandas data frame to HTML table. Replace default table styling with <a href="http://getbootstrap.com/css">Bootstrap</a> table styling.

#2. raw read QC from multiQC
rawqcFile = workdir + '/qc/multiqc/raw_multiqc_data/multiqc_fastqc.txt'
df = pd.read_csv(rawqcFile, sep='\t')
df1 = df.iloc[:,[0,4,5,6,7,8,9]].copy()
#df1.loc[:,'Sample'] = df1.iloc[:,0].str.replace('_R1_001', '', regex=False) #remove _R1_001 in the sample name
rawQC = df1.to_html(index=False).replace('<table border="1" class="dataframe">','<table class="table table-striped" data-name="mytable">') # use bootstrap styling
#get rawQC sample name
rawQCsample = df['Sample'].copy()

#3. read numbers in each QC step
numFile = rptdir + '/qcReadNumber.txt'
df = pd.read_csv(numFile, sep='\t', header = 0)
numQC = df.to_html(index=False).replace('<table border="1" class="dataframe">','<table class="table table-striped" data-name="mytable">') # use bootstrap styling
samples = df['Sample'].copy()

#4. clean read QC from multiQC
cleanqcFile = workdir + '/qc/multiqc/trimmed_multiqc_data/multiqc_fastqc.txt'
df = pd.read_csv(cleanqcFile, sep='\t')
df2 = df.iloc[:,[0,4,5,6,7,8,9]]
#df2.iloc[:,0].str.replace('_R1_001', '', regex=False)
cleanQC = df2.to_html(index=False).replace('<table border="1" class="dataframe">','<table class="table table-striped" data-name="mytable">') # use bootstrap styling
#get cleanQC sample name
cleanQCsample = df['Sample'].copy()

#5. assembly QC
#samples = df1['Sample'].copy()
#samples = list(map(lambda st: str.replace(st, '_R1_001', ''), samples)) 
#print(samples)
#put single quast togther

quastFile = workdir + '/qc/quast/' + samples[0] + '.quast/transposed_report.tsv'
df = pd.read_csv(quastFile, sep='\t', header = 0)
df.loc[0,'Assembly'] = samples[0]
for sample in samples[1:]:
    quastFile = workdir + '/qc/quast/' + sample + '.quast/transposed_report.tsv'
    df1 = pd.read_csv(quastFile, sep='\t', header=0)
    df1.loc[0,'Assembly'] = sample
    df = df.append(df1, sort = False)
df3 = df.iloc[:,[0,1,2,13,14,15,16,17,18,19,20,21]] 
#df1['Assembly'] = samples   
assemblyQC = df3.to_html(index=False).replace('<table border="1" class="dataframe">','<table class="table table-striped" data-name="mytable">') # use bootstrap styling

#6. Kraken2 report
kReport = []
for sample in samples:
    krakenFile = workdir + '/classification/' + sample + '.kraken2.report.txt'
    report = extractKraken(krakenFile, mpFile)
    for rpt in report:
        rpt.insert(0,sample)  #add it at the first position
        #print(rpt)
        kReport.append(rpt)
    #print(kReport)
#By default the keys of the dict become the DataFrame columns, Specify orient='index' to create the DataFrame using dictionary keys as rows:    
#df4 = pd.DataFrame.from_dict(kDict, orient='index')
df4 = pd.DataFrame(kReport, columns =['Sample', 'PercentageOfMappedReads', 'NumberOfReads', "TaxonId", 'Species', 'Domain']) 
#print(df4)
kraken2 = df4.to_html(index=False).replace('<table border="1" class="dataframe">','<table class="table table-striped" data-name="mytable">') # use bootstrap styling

#7. Kaiju report
kjReport = []
for sample in samples:
    kjFile = workdir + '/classification/' + sample + '.kaiju.table.txt'
    #with open(kjFile) as f:
    #    f.readline() #skip header
    #    for line in f:  
    #        line = line.rstrip()
    #        cells = line.split("\t")
    #        if 'Viruses' in cells[4] and float(cells[2]) >= threshold and cells[4] != 'unclassified' and '(non-viral)' not in cells[4] :
    #            kjReport.append([sample, cells[1],cells[2],cells[3],cells[4]])
    #        elif('Viruses' not in cells[4] and 'NA' not in cells[3] and float(cells[2]) >= 1000):
    #            kjReport.append([sample, cells[1],cells[2],cells[3],cells[4]])

    report = extractKaiju(kjFile, mpFile)
    for rpt in report:
        rpt.insert(0,sample)  #add it at the first position
        #print(rpt)
        kjReport.append(rpt)
    #print(kjReport)

df5 = pd.DataFrame(kjReport, columns =['Sample', 'PercentageOfMappedReads', 'NumberOfReads', 'TaxonId','Species', 'Domain'])  #'TaxonName'
#print(df5)
kaiju = df5.to_html(index=False).replace('<table border="1" class="dataframe">','<table class="table table-striped" data-name="mytable">') # use bootstrap styling

#8. Blastn

#remove empty files
blastnFiles = []
blastn = ""
bsamples = []
for sample in samples:
    blastnFile = workdir + '/annotation/' + sample + '.blastn.txt'
    # check if size of file is 0
    if os.stat(blastnFile).st_size > 0:
        blastnFiles.append(blastnFile)
        bsamples.append(sample)
#add blastn result        
if len(blastnFiles) > 0:
    df = pd.read_csv(blastnFiles[0], sep='\t',header=None)
    df.columns = ['Contig', 'Reference', 'Percent_identity', 'Alignment_length', 'Contig_length', 'Reference_length', 'Align_contig_start', 'Align_contig_end', 'Align_ref_start', 'Align_ref_end', 'E-value', 'Bitscore', 'Reference_name', 'Taxonomy_id', 'Species', 'Taxonomy_path']
    #print(df.dtypes)
    df.insert(loc=0, column='Sample', value=bsamples[0])  #add sample name at column one
    for i in range(1, len(blastnFiles)):
        df0 = pd.read_csv(blastnFiles[i], sep='\t',header=None)
        df0.columns = ['Contig', 'Reference', 'Percent_identity', 'Alignment_length', 'Contig_length', 'Reference_length', 'Align_contig_start', 'Align_contig_end', 'Align_ref_start', 'Align_ref_end', 'E-value', 'Bitscore', 'Reference_name', 'Taxonomy_id', 'Species', 'Taxonomy_path']
        df0.insert(loc=0, column='Sample', value=bsamples[i])
        #print(blastnFiles[i])
        df = df.append(df0, ignore_index = True)
        #print(df.shape[0])
    #display viruses and virod using different criteria
    df = df[(df['E-value'] <= blastnThreshold) | ((df['E-value'] <= 1e-30) & (df['Species'].str.contains("viroid")) )]
    df6 = df.iloc[:,[0,1,3,11,13,15]].copy()
    df6['E-value'] = df6['E-value'].map('{:,.2e}'.format)  #format scitific     
    #print(df6.head())
    blastn = df6.to_html(index=False, col_space=200).replace('<table border="1" class="dataframe">','<table class="table table-striped" data-name="mytable">') # use bootstrap styling

#9. Blastx
#remove empty files
blastxFiles = []
blastx = ""
bsamples = []
for sample in samples:
    blastxFile = workdir + '/annotation/' + sample + '.blastx.txt'
    # check if size of file is 0
    if os.stat(blastxFile).st_size > 0:
        blastxFiles.append(blastxFile)
        bsamples.append(sample)
#add blastn result        
if len(blastxFiles) > 0:
    df = pd.read_csv(blastxFiles[0], sep='\t', header=None)
    df.columns = ['Contig', 'Reference', 'Percent_identity', 'Alignment_length', 'Contig_length', 'Reference_length', 'Align_contig_start', 'Align_contig_end', 'Align_ref_start', 'Align_ref_end', 'E-value', 'Bitscore', 'Reference_name', 'Taxonomy_id', 'Species', 'Taxonomy_path']
    df = df[df.iloc[:,15].str.contains('virus',na=False)] #select only virus
    #df1 = df[df.iloc[:,10] <= threshold]
    df.insert(loc=0, column='Sample', value=bsamples[0])  #add sample name at column one
    #df7 = df1.iloc[:,[0,1,3,11,13,15]]
    for i in range(1,len(blastxFiles)):
        df0 = pd.read_csv(blastxFiles[i], sep='\t', header=None)
        df0.columns = ['Contig', 'Reference', 'Percent_identity', 'Alignment_length', 'Contig_length', 'Reference_length', 'Align_contig_start', 'Align_contig_end', 'Align_ref_start', 'Align_ref_end', 'E-value', 'Bitscore', 'Reference_name', 'Taxonomy_id', 'Species', 'Taxonomy_path']
        df0 = df0[df0.iloc[:,15].str.contains('virus',na=False)] #select only virus
        #df1 = df[df.iloc[:,10] <= threshold]
        df0.insert(loc=0, column='Sample', value=bsamples[i])
        #df1 = df1.iloc[:,[0,1,3,11,13,15]]
        #print(df1.head())
        df= df.append(df0, sort = False)
    df = df[df['E-value'] <= blastxThreshold]
    df7 = df.iloc[:,[0,1,3,11,13,15]].copy()
    df7['E-value'] = df7['E-value'].map('{:,.2e}'.format)  #format scitific     
    #print(df6.head())
    blastx = df7.to_html(index=False, col_space=200).replace('<table border="1" class="dataframe">','<table class="table table-striped" data-name="mytable">') # use bootstrap styling

#10. Selected contigs
selectFile = workdir + '/annotation/' + samples[0] + '.selectedRef.txt'
df = pd.read_csv(selectFile, sep='\t', header = 0)
df.insert(loc=0, column='Sample', value=samples[0])  #add sample name at column one
for sample in samples[1:]:
    selectFile = workdir + '/annotation/' + sample + '.selectedRef.txt'
    #print(selectFile)
    df1 = pd.read_csv(selectFile, sep='\t', header = 0)
    df1.insert(loc=0, column='Sample', value=sample)  #add sample name at column one
    df = pd.concat([df, df1], axis=0) # Stack the DataFrames on top of each other
#calculate genome/reference coverage
coverage = df.iloc[:,4]/df.iloc[:,6] 
df8 = df.iloc[:,[0,1,3,5,6,11,13,15,17]].copy()
df8['E-value'] = df8['E-value'].map('{:,.2e}'.format)  #format scitific 
df8.insert(loc=5, column='RefCoverage', value=coverage)  #add sample name at column one
selectContig = df8.to_html(index=False).replace('<table border="1" class="dataframe">','<table class="table table-striped" data-name="mytable">') # use bootstrap styling

#11. mapping reads to reference information
#rptFile = workdir + '/report/report.txt'
df = pd.read_csv(rptFile, sep='\t', header = 0, index_col=False, na_filter = False)
#print(df.head())
#df9 = df.iloc[:,[0,15,17,19,21,22,23,24,25,26,18,20]].copy()
df9 = df[["Sample","RefId","Species","Acronym","RefLen","MappedReads","RPKM","PercentMappedreads","PercentGenomeCovered","MeanCoverage","Taxonomy","RefTitle"]]
#col1 = df9.iloc[:,0].str.split('_', n=1, expand = True) #remove date from seqId
#df9.iloc[:,0] = col1[1]
#df9.rename(columns={'SeqID':'Sample'}, inplace=True) #change a column name
#print(df9.head())
#Add image html links to the table
df9c = df9.copy()  #ro prevent SettingWithCopyWarning
df9c['MappingGraph'] = "<a href=\"image/" + df9['Sample'] + "." + df9['RefId'] + ".coverage.png\" alt=\"" + df9['Sample'] + " " + df9['RefId'] + " coverage graph\" target=\"_blank\">" + df9['Sample'] + " " + df9['RefId'] +"</a>"
#if no maping graph file exists, set value to None
for i in range(len(df9c.index)):
    filename = rptdir + "/image/" + df9c.loc[i].Sample + "." + df9c.loc[i].RefId + ".coverage.png"
    #print(filename)
#    if not os.path.exists(filename):
    if not os.path.isfile(filename): 
        df9c.at[i,'MappingGraph'] = None
       
map2Ref = df9c.to_html(index=False, escape=False, col_space=200).replace('<table border="1" class="dataframe">','<table class="table table-striped" data-name="mytable">') # use bootstrap styling

#df10 = df.iloc[:,[1,8,15,17,27,28,29,30,31]].copy()
#df10 = df[["Sample","RefId","Species","BlastType","BlastQueryCover","BlastEvalue","BlastIdentity","BlastDescription"]]
df10 = df[["Sample","RefId","Species","BlastType","BlastEvalue","BlastIdentity","BlastDescription"]]
#df10['BlastEvalue'] = df10['BlastEvalue'].map('{:,.2e}'.format)  #format scitific
ncbiBlast = df10.to_html(index=False).replace('<table border="1" class="dataframe">','<table class="table table-striped" data-name="mytable">') # use bootstrap styling

#generate a fasta file for consensus (ref consensus or contig consensus)
# iterating over multiple columns
def makeFasta(sample, ref, abbr, seq):
    fasta = ""
    if not pd.isna(ref):   #ref != 'nan' is always true, use pd to check nan
        #print(sample, abbr , ref, seq )
        fasta = ">" + sample + "." + str(abbr) + "." + str(ref) + "\n" + str(seq)
    return fasta

result = [makeFasta(col[0], col[1], col[2], col[3]) for col in df[['Sample', 'RefId', 'Acronym','Sequence']].values]
#print(result)
#output final assemblies
outfile = rptdir + '/finalConsensus.txt'
fout = open(outfile, 'w')
for rst in result:
    if rst:
        fout.write(rst + "\n")
fout.close()

#12. mapping reads to novel virus contig/consensus
mapFile = rptdir + '/novelVirusMapping.txt'
map2Contig = ""
df11 = pd.DataFrame()
if os.path.exists(mapFile) and os.stat(mapFile).st_size != 0:
    df11 = pd.read_csv(mapFile, sep='\t')
    #print(df11.head())
    #add ref annotation
    df_temp = df10.iloc[:,[2,3]].copy().drop_duplicates() #RefId, Species; drop_duplicates() function removes all the duplicate rows and returns only unique rows.
    #df_temp['RefId'] = df_temp['RefId'] + ".contig"
    #df11 = pd.merge(df11, df_temp, on='RefId', how='left')
    df11['MappingGraph'] = "<a href=\"image/" + df11['Sample'] + "." + df11['RefId'] + ".coverage.png\" alt=\"" + df11['Sample'] + " " + df11['RefId'] + " coverage graph\" target=\"_blank\">" + df11['Sample'] + " " + df11['RefId'] +"</a>"
    map2Contig = df11.to_html(index=False, escape=False).replace('<table border="1" class="dataframe">','<table class="table table-striped" data-name="mytable">') # use bootstrap styling

#make html links
#Handle different raw sequence names: sample, sample_R1|R2, sample_R1|R2_001
rawFastQC = ""
count=0
rawNameFormat = 1
for sample in rawQCsample:
    fname1 = workdir + "/qc/raw_fastqc/" + sample + "_fastqc.html"
    fname2 = workdir + "/qc/raw_fastqc/" + sample.replace("_001","") + "_fastqc.html"  #for _001 name format
    if os.path.exists(fname1):
        rawFastQC += "<a href=\"html/" + sample + "_fastqc.html\">" + sample + "   </a> &emsp;&emsp;" 
    elif os.path.exists(fname2):
        rawFastQC += "<a href=\"html/" + sample.replace("_001","") + "_fastqc.html\">" + sample.replace("_001","") + "   </a> &emsp;&emsp;" 
    else:
        rawNameFormat = 0
        break
    count += 1
    if count%2 == 0:
        rawFastQC += "<br />" #add break
#print(rawFastQC)

if not rawNameFormat:
    rawFastQC = ""
    count=0
    for sample in samples:
        #check a file exists
        if os.path.exists(workdir + "/qc/raw_fastqc/" + sample + "_fastqc.html"):
            rawFastQC += "<a href=\"html/" + sample + "_fastqc.html\">" + sample + "   </a> &emsp;&emsp;" 
        count += 1
        if count%2 == 0:
            rawFastQC += "<br />" #add break

#for trimmed sequence fastQC    
trimFastQC = ""
count=0
for sample in cleanQCsample:
    trimFastQC += "<a href=\"html/" + sample + ".trimmed_fastqc.html\">" + sample + "   </a> &emsp;&emsp;"
    count += 1
    if count%2 == 0:
        trimFastQC += "<br />" #add break

quastQC = ""
krakenReport = ""
kaijuReport = ""
blastnReport = ""
blastxReport = ""
blastntReport = ""
blastnrReport = ""
count=0
for sample in samples:
    #rawFastQC += "<a href=\"html/" + sample + "_fastqc.html\">" + sample + "   </a> " 
    #trimFastQC += "<a href=\"html/" + sample + ".trimmed_fastqc.html\">" + sample + "   </a> " 
    quastQC += "<a href=\"html/" + sample + ".quast.html\">" + sample + "   </a> &emsp;&emsp;" 
    krakenReport += "<a href=\"html/" + sample + ".kraken2.report.html\">" + sample + "   </a> &emsp;&emsp;"
    kaijuReport += "<a href=\"html/" + sample + ".kaiju_krona.html\">" + sample + "   </a> &emsp;&emsp;"
    blastnReport += "<a href=\"html/" + sample + ".blastn.krona.html\">" + sample + "   </a> &emsp;&emsp;"
    blastxReport += "<a href=\"html/" + sample + ".blastx.krona.html\">" + sample + "   </a> &emsp;&emsp;"

    blastntFile = rptdir + '/html/' + sample + '.blastnt.krona.html'
    blastnrFile = rptdir + '/html/' + sample + '.blastnr.krona.html'

    if os.path.exists(blastntFile):
        blastntReport += "<a href=\"html/" + sample + ".blastnt.krona.html\">" + sample + "   </a> &emsp;&emsp;"
    
    if os.path.exists(blastnrFile):
        blastnrReport += "<a href=\"html/" + sample + ".blastnr.krona.html\">" + sample + "   </a> &emsp;&emsp;"

    count += 1
    if count%2 == 0:
        #rawFastQC += "<br />" #add break
        #trimFastQC += "<br />"
        quastQC += "<br />"
        krakenReport += "<br />"
        kaijuReport += "<br />"
        blastnReport += "<br />"
        blastxReport += "<br />"
        if blastntReport != "":
            #print("blastntreport is not empty", blastntReport)
            blastntReport += "<br />"
        if blastnrReport != "":
            blastnrReport += "<br />"

        
#print(rawFastQC)
#print(trimFastQC)
#print(quastQC)
#print(krakenReport)
#print(kaijuReport)

#make image links
#To open in the same tab use: target = "_self"
def makeLink(sample, ref):
    link = ""
    #if not pd.isna(ref):  #np.isnan(float('nan')); math.isnan(float('nan'))
    if ref != "nan" and ref != '':
        link = "<a href=\"image/" + sample + "." + ref + ".coverage.png\" alt=\"" + sample + " " + ref + " coverage graph\" target=\"_blank\">" + sample + " " + ref +"</a>"
    return link

refLink0 = [makeLink(col[0], col[1]) for col in df9[['Sample', 'RefId']].values]
refLink = [x for x in refLink0 if x] #remove empty strings
refLinkStr = "<br />".join(refLink)

contigLinkStr = ""
if len(df11)>1:
    contigLink0 = [makeLink(col[0], col[1]) for col in df11[['Sample', 'RefId']].values]
    contigLink = [x for x in contigLink0 if x] #remove empty strings
    contigLinkStr = "<br />".join(contigLink)
#print(refLinkStr)
#print(contigLinkStr)

#copy files to report folder
#copy html files
htmlDir = rptdir + '/html'
if not os.path.isdir(htmlDir): #check a folder exists
    os.mkdir(htmlDir) #make a folder
imgDir = rptdir + '/image'
if not os.path.isdir(imgDir): #check a folder exists
    os.mkdir(imgDir) #make a folder

#copy files
source_dir = workdir + '/classification'
for filename in glob.glob(os.path.join(source_dir, '*.html')):
    try:
        shutil.copy(filename, htmlDir)
    except PermissionError:
        print("Destination file exits!")
source_dir = workdir + '/qc/raw_fastqc'
for filename in glob.glob(os.path.join(source_dir, '*.html')):
    shutil.copy(filename, htmlDir)
source_dir = workdir + '/qc/trimmed_fastqc'
for filename in glob.glob(os.path.join(source_dir, '*.html')):
    shutil.copy(filename, htmlDir)   
source_dir = workdir + '/qc/multiqc'
for filename in glob.glob(os.path.join(source_dir, '*.html')):
    shutil.copy(filename, htmlDir)

#copy blast classification krona files
source_dir = workdir + '/annotation'
for filename in glob.glob(os.path.join(source_dir, '*.html')):
    try:
        shutil.copy(filename, htmlDir)
    except PermissionError:
        print("Destination file exits!")

for sample in samples:
    fromName = workdir + '/qc/quast/' + sample + '.quast/report.html'
    toName = htmlDir + '/' + sample + '.quast.html'
    shutil.copy(fromName, toName)

source_dir = workdir + '/mapping/map2Ref'
for filename in glob.glob(os.path.join(source_dir, '*.png')):
    shutil.copy(filename, imgDir) 
source_dir = workdir + '/novelVirus/map2Contig'
if os.path.exists(source_dir):
    for filename in glob.glob(os.path.join(source_dir, '*.png')):
        shutil.copy(filename, imgDir)


#get kraken2 and kaiju overlapped viruses
#df4 = pd.DataFrame(kReport, columns =['Sample', 'PercentageOfMappedReads', 'NumberOfReads', "TaxonId", 'Species', 'Domain']) 
#df5 = pd.DataFrame(kjReport, columns =['Sample', 'PercentageOfMappedReads', 'NumberOfReads', 'TaxonId','TaxonName']) 
df4['Sample_taxon'] = df4.Sample.str.cat(df4.TaxonId.astype(str), sep='_')
df5['Sample_taxon'] = df5.Sample.str.cat(df5.TaxonId.astype(str), sep='_')
kk_overlap = pd.merge(df4, df5, how='inner', on=['Sample_taxon'])
#print(kk_overlap.head())
kk_overlap1 = kk_overlap[['Sample_x','NumberOfReads_x', 'NumberOfReads_y','Species_x','Domain_x']] #'TaxonName'
# Change the column names 
kk_overlap1.columns =['Sample','NumberOfReads_Kraken', 'NumberOfReads_Kaiju','Species','Domain'] #,'TaxonName'
kk_overlapped = kk_overlap1.to_html(index=False).replace('<table border="1" class="dataframe">','<table class="table table-striped" data-name="mytable">') # use bootstrap styling
    
### Step 2: Generate HTML report as a string and write to file
html_string = '''
<!DOCTYPE html>
<html lang="en" style="height: 100%;">
<head>
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.1/css/bootstrap.min.css">
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <style>
    * {
        box-sizing: border-box;
       }
    body {
        font-family: Arial, Helvetica, sans-serif;
    }

    /* Style the header */
    header {
        background-color: #666;
        padding: 10px;
        text-align: center;
        font-size: 35px;
        color: white;
    }
    .sidenav {
      height: 100%;
      width: 200px;
      position: fixed;
      z-index: 1;
      top: 0;
      left: 0;
      background-color: #eee;
      overflow-x: hidden;
      padding-top: 20px;
    }

    .sidenav a {
      padding: 6px 8px 6px;
      text-decoration: none;
      font-size: 16px;
      color: #2196F3;
      display: block;
    }

    .sidenav a:hover {
      color: #064579;
    }

    .main {
      margin-left: 200px; /* Same as the width of the sidenav for sticky navigation*/
      font-size: 16px; /* Increased text to enable scrolling */
      padding: 0px 10px;
    }
    
    @media screen and (max-height: 450px) {
      .sidenav {padding-top: 15px;}
      .sidenav a {font-size: 12px;}
    }
    
    article {
        float: left;
        padding: 20px;
        background-color: #f1f1f1;
        width:100%;
    }


    </style>

    <script>
    function myFunction() { // for filter
        var input, filter, table, tr, td, i,alltables;
        alltables = document.querySelectorAll("table[data-name=mytable]");
        input = document.getElementById("myInput");
        filter = input.value.toUpperCase();
        alltables.forEach(function(table){
            tr = table.getElementsByTagName("tr");
            for (i = 0; i < tr.length; i++) {
                td = tr[i].getElementsByTagName("td")[0];
                if (td) {
                    if (td.innerHTML.toUpperCase().indexOf(filter) > -1) {
                        tr[i].style.display = "";
                    } else {
                        tr[i].style.display = "none";
                    }
                }
            }
        });
    }
    </script>
    
    </head>
    
    <body>
        <header>
            <h1>HTS Report</h1>
        </header>

        <section>
            <!-- Sidebar -->
            <div class="sidenav">
                    <a href="#summary">Summary</a>
                    <h4>Sample QC</h4>
                    <a href="#rawQC">Raw reads QC</a>
                    <a href="#cleanQC">Clean reads QC</a>
                    <a href="#numQC">Reads number track</a>
                    <a href="#assemblyQC">Assembly QC</a>
                    <h4>Read classification</h4>
                    <a href="#kraken2">Kraken2 (NCBI nt)</a>
                    <a href="#kaiju">Kaiju (Pathogen nr)</a>
                    <a href="#kk_overlap">Kraken2 and Kaiju overlap</a>
                    <h4>Contig annotation</h4>
                    <a href="#blastn">Blastn (viral refSeq)</a>
                    <a href="#blastx">Blastx (RVDB)</a>
                    <a href="#selectContig">Selected contigs</a>
                    <h4>Mapping</h4>
                    <a href="#mapping">Map reads to reference</a>
                    <h4>Novel virus</h4>
                    <a href="#novel">Map reads to contig</a>                    
            </div>
            
            <div class="main">
            <br>
            <div class="card-body row no-gutters align-items-center pull-right">
                  <label>Filter sample </label>
                  <input type="text" id="myInput" onkeyup="myFunction()" placeholder="A sample name " title="Type in a sample name">
            </div>

            <article id=summary>
            <h1>Virus NCBI blast</h1>''' + ncbiBlast + '''
            <h3><a href="finalConsensus.txt" target="_blank">View consensus</a></h3><br />
            </article>
            
            <article id='rawQC'>
            <h1>Raw reads QC</h1>''' + rawQC + '''
            
            <h3>View raw reads fastQC report:</h3> ''' + rawFastQC + ''' <br><br>
            <a href="html/raw_multiqc.html">All samples raw reads multiQC</a>
            </article>
            
            <article id='cleanQC'>
            <h1>Clean reads QC</h1>''' + cleanQC + '''
            <h3>View clean reads fastQC report:</h3> ''' + trimFastQC + ''' <br><br>
            <a href="html/trimmed_multiqc.html">All samples clean reads multiQC</a>
            </article> 

            <article id='numQC'>
            <h1>Read numbers in each QC step</h1>''' + numQC + '''
            </article> 
            
            <article id='assemblyQC'>
            <h1>Assembly QC</h1>''' + assemblyQC + '''
            <h3>View assembly QC quast report:</h3> ''' + quastQC + ''' <br><br>
            <a href="html/quast_multiqc.html">All samples quast multiQC</a>
            </article> 


            <article id='kraken2'>
            <h1>Kraken2 report</h1>''' + kraken2 + '''
            <h3>View Kraken2 krona report:</h3> ''' + krakenReport + ''' 
            </article> 

            <article id='kaiju'>
            <h1>Kaiju report</h1>''' + kaiju + '''
            <h3>View Kaiju krona report:</h3> ''' + kaijuReport + ''' 
            </article>
            
            <article id='kk_overlap'>
            <h1>Pathogens detected by both Kraken2 and Kaiju:</h1> ''' + kk_overlapped + ''' 
            </article>
            
            <article id='blastn'>
            <!-- <h1>Blastn: contigs against NCBI viral references</h1>''' + blastn + ''' -->
            <h1>View viral blastn krona report:</h1> ''' + blastnReport + ''' 
            </article>
                
            <article id='blastx'>
            <!-- <h1>Blastx: contigs against Reference Viral Database (RVDB)</h1>''' + blastx + ''' -->
            <h1>View viral blastx krona report:</h1> ''' + blastxReport + ''' 
            </article>
            '''
if blastntReport:
    html_string += '''<article id='blastnt'>
                <h1>View contigs blastn against NCBI nt krona report:</h1> ''' + blastntReport + ''' 
                </article>'''
if blastnrReport:
    html_string += '''<article id='blastnr'>
            <h1>View contigs blastx against NCBI nr krona report:</h1> ''' + blastnrReport + ''' 
            </article>'''
            
html_string += '''<article id='selectContig'>               
            <h1>Selected contigs and their annotation</h1>''' + selectContig + '''
            </article>
            
            <article id='mapping'>
            <h1>Mapping reads to viral references</h1>''' + map2Ref + '''
            <!--
            <h3>View reads mapping to reference graph:</h3> ''' + refLinkStr + ''' 
            -->
            </article>

            
            <article id='novel'>
            <h1>Mapping reads to novel contigs</h1>
            '''

if contigLinkStr: #skip no novel contigs
    #html_string += map2Contig + ''' <h3>View reads mapping to assembly graph:</h3> ''' + contigLinkStr
    html_string += map2Contig
else:
    html_string += '''<h3>No novel contigs</h3> '''
    
html_string += '''</article>
            </div>
    </body>
</html>'''


# Finally, write the html string to a local file.
fout = open(outFile,'w')
fout.write(html_string)
fout.close()