#! /usr/bin/env python
######################################################################
#Alex Hu <xiaojun.hu@usda.gov>
#Updated: 07/15/2021
#This program selects top references from kraken2, kaiju, blastn and blastx resultx
#
#Required Parameters:
#   -n, --blastn X..........................Blastn table file
#   -x, --blastx X..........................Blastx table file
#   -k, --kraken X..........................Kraken2 report
#   -j, --kaiju X...........................Kaiju table output
#   -c, --contig X..........................contig file
#   -s, --select X..........................output file for selected references
#   -a, --all X.............................output file for merged blastn and blastx
#   -l, --listn X...........................output file for selected reference names from blastn
#   -r, --listx X...........................output file for selected reference names from blastx
######################################################################

import os, sys, argparse
from Bio import SeqIO

kraken_viral_cutoff = 5 #mapped viral reads cut off
kaiju_cutoff = 100 # use reads number as cutoff instead of percentage 0.05
evalue_cutoff_blastn = 1e-30 #blastn evalue cut off 
evalue_cutoff_blastx = 1e-40 # blastx evalue for reference filter
identity_cutoff = 80  #blastn identity cutoff
virod_cutoff = 1e-10  #E-value cutoff for viroid blastn
 
#################################### Get arguments function ####################################
def parseArguments():
    """
    Get arguments
    """
    parser = argparse.ArgumentParser(description = "Select top reference from blastn and blastx result.")
    parser = argparse.ArgumentParser(description = "Select top reference from blastn and blastx result.")
    parser.add_argument('-n','--blastn',dest='blastn_file', required=True, default="", help='Blastn against nt database result')
    parser.add_argument('-x','--blastx',dest='blastx_file', required=True, default="", help='Blastx against protein database (nr or virus RVDB) result')
    parser.add_argument('-k','--kraken',dest='kraken_file', required=True, default="", help='Kraken2 report')
    parser.add_argument('-j','--kaiju',dest='kaiju_file', required=True, default="", help='Kaiju table output from kaiju2table')
    parser.add_argument('-c','--contig',dest='contig_file', required=True, default="", help='Contig file from assembler output, ex. contigs.txt ')
    parser.add_argument('-s','--select',dest='out_selectedContig', required=True, default="", help='output file for selected references' )
    parser.add_argument('-a','--all',dest='output_mergedBlast', required=True, default="", help='output file for merged blastn and blastx' )
    parser.add_argument('-l','--listn',dest='output_blastnRef', required=True, default="", help='selected reference name list from blastn' )
    parser.add_argument('-r','--listx',dest='output_blastxRef', required=True, default="", help='selected reference name list from blastx' )
    parser.set_defaults(append=False)
    return parser.parse_args()

#################################### Get blast table result function ####################################
#blastn format: qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore stitle taxonid taxonpath
#TRINITY_DN31_c0_g1_i1	gi|746718929|gb|KM025045.1|	74.281	7652	9923	9624	984	8560	750	8327	0.0	3090	Sorghum mosaic virus isolate YN-YT168, complete genome	32619	Sorghum mosaic virus	Sorghum mosaic virus	Viruses
def getBlastn(blastnFile):
    """
    Get blastn result, if the records have the same query id [0] and taxonomy id [13], keep the higher bit score one
    """
    blastnDict = {}
    idDict = {} #keep unique id
    with open(blastnFile) as f:
        for line in f:
            cells = line.rstrip().split('\t')
            #gb = cells[1]
            #if("gb|" in cells[1]):
                #get gi and gb
            #    gids = cells[1].split('|') #gi|746718929|gb|KM025045.1|
            #    gi = gids[1]  #746718929
            #    gb = gid[3]
                #ngb = gids[3].split('.')[0]  #get genbank id, KM025045
            if len(cells) < 14: #for missing taxonomy
                cells.extend(["","0"])
            id = cells[0]+"|"+cells[13] #query id + tax id 

            if((id not in idDict) or float(idDict[id]) < float(cells[11])): #keep the higher bitscore one
                idDict[id] = cells[11]
                blastnDict[cells[0]] = cells

    return blastnDict

#################################### Get blast table result function ####################################
#blast format: qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore stitle taxonid taxonpath
#NODE_33_length_4928_cov_612.276850	NC_003491.1	78.378	259	4928	5722	1013	1268	4030	3775	8.50e-38	163	NC_003491.1 Beet mild yellowing virus, complete genome	156690	Viruses; Riboviria; Luteoviridae; Polerovirus; Beet mild yellowing virus;
#NODE_289_length_2910_cov_8.945994	acc|GENBANK|QBK89270.1|GENBANK|MK500411|metallopeptidase	28.7	352	2910	369	205	1134	4	347	5.2e-26	124.4	acc|GENBANK|QBK89270.1|GENBANK|MK500411|metallopeptidase family M24 [Mimivirus LCMiAC02]	2506609	Viruses; Mimiviridae; unclassified Mimiviridae; Mimivirus LCMiAC02;
def getBlastx(blastFile):
    """
    Get blastx result, keep the top one (highest bitscore) for a query
    """
    blastDict = {}
    bitscoreDict = {} #keep unique id: query id + taxon id
    with open(blastFile) as f:
        for line in f:
            cells = line.rstrip().split('\t')
            #cells[10] = float(cells[10]) #convert string to sci number
            id = cells[0] + "|" + cells[13] #query id + tax id
            if len(cells) < 14: #for missing taxonomy
                cells.extend(["","0"])            
            if(id not in blastDict or float(bitscoreDict[id]) < float(cells[11])): #only keep top one (the highest bitscore) for a query
                #blastDict.setdefault(cells[0], []).append([cells])
                blastDict[cells[0]] = cells
                bitscoreDict[id] = cells[11]  #bitscore
    return blastDict

#################################### Get kaiju table result function ####################################
#file	percent	reads	taxon_id	taxon_name
#OCT19-Q49548.kaiju.out	0.245875	8023	351495	Viruses;NA;NA;NA;Partitiviridae;NA;Raphanus sativus cryptic virus 2;
def getKaiju(kaijuFile):
    """
    Get kaiju result, use reads > 100  #use >5% as 'percent' cut-off 
    """
    kaijuDict = {}
    with open(kaijuFile) as f:
        next(f)  #skip header
        for line in f:
            cells = line.rstrip().split('\t')
            if(float(cells[2]) >= 10 and cells[3] != 'NA'): #reads >= 10 and has taxon, keep it  #cells[1]>0.05
                kaijuDict[cells[3]] = cells[2]  #taxon_id, reads
    return kaijuDict

#################################### Get Kraken2 report function ####################################
#percent	reads_clade_rotted  reads_taxon rank	taxon_id	scientific_name
# 15.55  4683889 4683889 U       0       unclassified
# 84.45  25434922        102908  R       1       root
#  0.01  1520    0       D       10239     Viruses
def getKraken(kFile):
    """
    Extract Viruses(10239) from Kraken2 report
    keep taxon_id  and assigned reads
    """
    #20.06  3234437 3234437 U       0       unclassified
    #kraken_viral_cutoff = 50  
    krakenDict = {}
    vFlag = 0
    with open(kFile) as f:
        for line in f:
            line = line.rstrip()
            #cells = line.split("\t")
            cells = list(map(str.strip, line.split("\t")))
            if cells[5] == 'uncultured organism':
                continue
            if cells[4] == '10239':  #Viruses clade starts
                vFlag = 1
                continue    
            if vFlag == 1 and cells[3].startswith('S') and float(cells[2]) >= kraken_viral_cutoff and 'synthetic construct' not in cells[5]: 
                #[cells[0],cells[2],cells[4], cells[5],'Viruses'] = ['PercentageOfMappedReads', 'NumberOfReads', "TaxonId", 'Species', 'Domain']
                krakenDict[cells[4]] = cells[2] #taxon_id => reads

    return krakenDict

#################################### Select top reference from both blastx and blastn  ####################################
def selectTop(blastnDict, blastxDict, krakenDict, kaijuDict, evalue_cutoff_blastn, evalue_cutoff_blastx):
    """
    Select top candidate references based on blastn E-value (pathogen:0.0, viroid <1e-30) or blastx (e-value == 0.0 or found in both blastx and blastn) or both kaiju and blastx/blastn found it
    """
    outDict = {} #save selected contigs ant their blast result
    taxonDict = {}   #save selected taxons
    #1. select viroid and top virus from blastn
    for id in blastnDict:  # id = contig name
        #1)select a virus from blastn, if it's evalue=0 and identity >= 80
        if( float(blastnDict[id][10]) == 0.0 and float(blastnDict[id][2]) >= identity_cutoff):  #and (blastnDict[id][-1] == "Viruses" or blastnDict[id][-1] == "Bacteria" or blastnDict[id][-1] == "Fungi" or blastnDict[id][-1] == "Archaea") ): #only pathogens
            outDict[id] = "\t".join(blastnDict[id]) + "\tBlastn"
            taxonDict[blastnDict[id][13]] = 1
            #print("Blastn ", id, " ", blastnDict[id][13])
            #del blastnDict[id]  #delete the key, value
        #2)select a viroid if its evalue<=1e-10 and have viroid in the annotation
        if("viroid" in blastnDict[id][12] and float(blastnDict[id][10]) <= virod_cutoff):  
            outDict[id] = "\t".join(blastnDict[id]) + "\tBlastn"
            taxonDict[blastnDict[id][13]] = 1
            #print("Blastn viroid", id, " ", blastnDict[id][13])

    #2. select a virus if a viral protein is identified (evalue=0)
    for id in blastxDict:
        #if matched protein is found, keep it
        if float(blastxDict[id][10]) == 0.0 and id not in outDict:
            #print("Blastx ", id, " ", blastxDict[id][13])
            taxonDict[blastxDict[id][13]] = 1            
            outDict[id] = "\t".join(blastxDict[id]) + "\tBlastx"  #+ "\t"*5, blastxDict[id][16]

    #3. if virus is found in both kraken2 and blastn
    for tid in krakenDict:  #taxon id in Kaiju
        for id in blastnDict:
            if id not in outDict and blastnDict[id][13] not in taxonDict:  #if a contig or a taxon is not previously selected 
                #if taxon_id in both kraken2 and blastn, and blastn evalue < cutoff, keep it
                #print(id, tid, blastnDict[id][13], blastnDict[id][10])
                if tid == blastnDict[id][13] and float(blastnDict[id][10]) <= evalue_cutoff_blastn: #if the same taxon id
                    outDict[id] = "\t".join(blastnDict[id]) + "\tBlastn"
                    taxonDict[blastnDict[id][13]] = 1
                    #print("Kraken2 + blastn ", id, " ", blastnDict[id][13])

    #4. for the virus having several segments, if one fragment passes the threshold, the other fragments should be included, even the contig e-value is higher. The contig with highest bit score one is kept
    refTaxonDict = {}
    for bid in blastnDict:  
        for oid in outDict:
            #if taxon id was selected, and reference was not selected, keep the reference with the higher bitscore
            selected = outDict[oid].split("\t")
            if blastnDict[bid][13] == selected[13]: #taxonId equal
                if blastnDict[bid][1] != selected[1]: #refId different
                    rtId = selected[1] + "_" + selected[13] #=refId + taxonId
                    #refTaxonDict[rtId] = selected  
                    #keep the contig with the higher bitScore [11] 
                    if rtId not in refTaxonDict or float(refTaxonDict[rtId][11]) < float(blastnDict[bid][11]):
                        refTaxonDict[rtId] = blastnDict[bid]
                        #print(blastnDict[bid])
    #print(outDict.keys())                    
    for rtId in refTaxonDict:
        newid = refTaxonDict[rtId][0] #contig name
        if newid not in outDict:
            outDict[newid] = "\t".join(refTaxonDict[rtId]) + "\tBlastn"
            #print("Virus segments ", newid)
        
    #5. if virus is found in both kaiju and blastx
    for id in blastxDict:
        #if matched protein is found, keep it
       for tid in kaijuDict:
            if id not in outDict and blastxDict[id][13] not in taxonDict and float(kaijuDict[tid]) >= kaiju_cutoff: #if Kaiju reads# >=100:
                #if taxon_id in both kaiju and blastx and blastx evalue < cutoff, keep it
                if tid == blastxDict[id][13] and float(blastxDict[id][10]) <= evalue_cutoff_blastx: #if the same taxon id
                    #print("Kaiju + Blastx ", id, " ", blastxDict[id][13])
                    taxonDict[blastxDict[id][13]] = 1    
                    outDict[id] = "\t".join(blastxDict[id]) + "\tBlastx"  #+ "\t"*5, blastxDict[id][16]
    
    #Note: can't use contig coverage filter, some virus titer is very low
    #del outDict[id]

    return outDict

#################################### Get reference names from selected contigs and interesting viruses  ####################################
def getRefName(selectDict):
    """
    Get reference names for a taxon id, keep the highest bitscore one 
    """
    idDictn = {} #for blastn
    idDictx = {} #for blastx
    bitscoreDictn = {}
    bitscoreDictx = {}
    taxDict = {}
    for id in selectDict:
        cells = selectDict[id].split("\t")
        #for blastn
        if cells[-1] == 'Blastn':
            rid = cells[1] #reference
            if "." in cells[1]:
                rid = cells[1].split(".")[0]
            #combined id: refid + taxid, for some viruses having several segments
            cid = rid + "_" + cells[13]
            if (cid not in idDictn or float(cells[11]) > bitscoreDictn[cid]):#keep the taxon having the higher bitscore
                #idDictn[cid] = id
                idDictn[cid] = rid
                #print(id)
                bitscoreDictn[cid] = float(cells[11]) #refid+taxonid, bitscore          
                taxDict[cells[13]] = 1

    #search blastx if a different taxid is found, possible a new virus
    for id in selectDict:
        cells = selectDict[id].split("\t")
        rid = cells[1]
        if "acc|" in cells[1]:
            rid = cells[1].split("|")[-2]    
        #for blastx, taxid not in blastn       
        if (cells[-1] == 'Blastx' and float(cells[10]) <= evalue_cutoff_blastx and (float(cells[10]) == 0.0 or cells[13] not in taxDict)):
            if (cells[13] not in idDictx or float(cells[11]) > bitscoreDictx[cells[13]]): #keep the taxon having the higher bitscore  
                cid = rid + "_" + cells[13]
                idDictx[cid] = rid #one ref id per tax id
                bitscoreDictx[cells[13]] = float(cells[11]) #taxon_id, bitscore

    return idDictn, idDictx


#################################### merge blastx and blastn  ####################################
def mergeBlast(blastnDict, blastxDict, contig_file):
    """
    Merge blastn and blastx
    """
    outDict = {}
    idDict = {}

    #put contig sequences in a dict
    contig_dict = SeqIO.index(contig_file, "fasta")
    #print(contig_dict)
    #get all sequence id
    for id in blastnDict:
        idDict[id] = 1  #only in blastn
    for id in blastxDict:
        if id in idDict:
            idDict[id] = 2  #in both blastn and blastx
        else:
            idDict[id] = 0  #only in blastx
    
    #merge blastn and blastx
    empty = "\t"*16
    for id in idDict:
        seq = contig_dict[id]
        if idDict[id] == 2:  #id in both blastx and blastn
            outDict[id] = "\t".join(blastxDict[id]) + "\t" + "\t".join(blastnDict[id]) + "\t" + seq.seq
        elif idDict[id] == 1: #id in only blastn
            outDict[id] = empty + "\t".join(blastnDict[id]) + "\t" + seq.seq
        else:  #id in only blastx
            outDict[id] = "\t".join(blastxDict[id]) + "\t" + empty + seq.seq

    contig_dict.close()
    
    return outDict

def main():
    ### Input arguments
    options = parseArguments()
    blastn_file = options.blastn_file
    blastx_file = options.blastx_file
    kraken_file = options.kraken_file
    kaiju_file = options.kaiju_file
    output_file1 = options.out_selectedContig
    output_file2 = options.output_mergedBlast
    output_file3 = options.output_blastnRef
    output_file4 = options.output_blastxRef
    contig_file = options.contig_file
    
    blastnDict = getBlastn(blastn_file)
    #print(blastnDict)
    blastxDict = getBlastx(blastx_file)
    krakenDict = getKraken(kraken_file)
    kaijuDict = getKaiju(kaiju_file)

    outDict1 = selectTop(blastnDict, blastxDict, krakenDict, kaijuDict, evalue_cutoff_blastn, evalue_cutoff_blastx)
    #output selected blast result
    fout1 = open(output_file1, 'w')
    fout1.write("Contig\tReference\tPercent_identity\tAlignment_length\tContig_length\tReference_length\tAlign_contig_start\tAlign_contig_end\tAlign_ref_start\tAlign_ref_end\tE-value\tBitscore\tReference_name\tTaxonomy_id\tSpecies\tTaxonomy_path\tblastType\n") 
    for id in outDict1:
        fout1.write("%s\n" % outDict1[id])
    fout1.close()

    outDict2 = mergeBlast(blastnDict, blastxDict, contig_file)
    #output selected blast result
    fout2 = open(output_file2, 'w')
    fout2.write("Contig_blastx\tReference_blastx\tPercent_identity_blastx\tAlignment_length_blastx\tContig_length_blastx\tReference_length_blastx\tAlign_contig_start_blastx\tAlign_contig_end_blastx\tAlign_ref_start_blastx\tAlign_ref_end_blastx\tE-value_blastx\tBitscore_blastx\tReference_name_blastx\tTaxonomy_id_blastx\tSpecies_blastx\tTaxonomy_path_blastx\t")
    fout2.write("Contig_blastn\tReference_blastn\tPercent_identity_blastn\tAlignment_length_blastn\tContig_length_blastn\tReference_length_blastn\tAlign_contig_start_blastn\tAlign_contig_end_blastn\tAlign_ref_start_blastn\tAlign_ref_end_blastn\tE-value_blastn\tBitscore_blastn\tReference_name_blastn\tTaxonomy_id_blastn\tSpecies_blastn\tTaxonomy_path_blastn\tContig\n") 
    for id in outDict2:
        fout2.write("%s\n" % (outDict2[id]))
    fout2.close()

    #output reference list, one per line
    (outDict3,outDict4) = getRefName(outDict1)
    #output selected references from blastn
    fout3 = open(output_file3, 'w')
    for id in outDict3:
        fout3.write("%s\n" % (outDict3[id]))
    fout3.close()

    #output selected references from blastx
    fout4 = open(output_file4, 'w')
    for id in outDict4:
        fout4.write("%s\n" % (outDict4[id]))
    fout4.close()
    
if __name__ == '__main__':
	main()
