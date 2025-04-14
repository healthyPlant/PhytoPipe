#! /usr/bin/env python
#########################################################################################
#Katie Ko
#Updated: 7/28/2022
#This program prepares a nucleotide and protein database using plant virus sequences maintained 
#   with the VirusDetect software package (http://virusdetect.feilab.net/cgi-bin/virusdetect/index.cgi).
#   The database can be used to identify plant viruses present in RNA-seq samples.
#Usage: python createPlantViralDb.py [-db] [-n2] [-p2] [-nd] [-nm] [-ni] [-a] [-g] [-i]
#
#Required Arguments:
#   -db X.....................................full path to desired database directory
#   -n2 X.....................................full path to taxonomy nucl_gb.accession2taxid file
#   -p2 X.....................................full path to taxonomy prot.accession2taxid file
#   -nd X.....................................full path to taxonomy nodes.dmp file
#   -nm X.....................................full path to taxonomy names.dmp file
#   -ni X.....................................full path to fasta file containg VirusDetect nucleotide database sequences
#   -a X......................................full path to fasta file containg VirusDetect protein database sequences
#   -g X......................................full path to VirusDetect Genbank file (vrl_genbank.info.gz) for annotations
#   -i X......................................full path to VirusDetect file (vrl_idmapping.gz) containing ID mappings
#
#Optional Arguments:
#   --nuc_prefix X............................specify desired prefix for nucleotide database output files
#   --prot_prefix X...........................specify desired prefix for protein database output files
#
#Output
# - This program generates an annotated nucleotide fasta file, a nucleotide gb_taxon.txt file, a dmnd file,
#   and a protein gb_taxon.txt file for use as a Plant Virus Database with PhytoPipe
#########################################################################################

from Bio import SeqIO
import argparse, subprocess, re
import time
import gzip

################################### Generate gb_taxon.txt file ###################################
#make_gb_taxon_file
#Usage: generates gb_taxon.txt file containing accession, taxonID, and taxonomy of each virus sequence in three tab-delimited columns
#Input: 
# - dictionary containing sequence accessions as keys and matching taxonIDs as values
# - dictionary generated from taxonomy nodes.dmp file with key as taxonID and value as parent node taxonID of key
# - dictionary generated from taxonomy names.dmp file with key as taxonID and value as scientific name of key
# - full path of desired gb_taxon.txt file output
#Generates:
# - gb_taxon.txt file for use with Phytopipe, format:
#  [accession]     [taxonID]        [taxonomy]
#   KY203336        2057945     Viruses; unclassified viruses; Pleurochrysis sp. endemic virus 1b;
#
def make_gb_taxon_file(accession_to_taxon_dict, nodes_dict, names_dict, gb_taxon_file_name):

    gb_taxon_file = open(gb_taxon_file_name, "w")   #open gb_taxon_file containing accession, taxonID, and taxonomy in tab delimited columns for each sequence in database

    for accession in accession_to_taxon_dict:                         

        lineage = ""

        taxonId = accession_to_taxon_dict[accession]    #obtain taxonID for each accession in database

        curr_id = taxonId   #set curr_id to accession's taxonID

        while curr_id in names_dict and curr_id != "1":     #recurses from accession's taxonID to root taxonID (1) while curr_id is in names_dict

            lineage = names_dict[curr_id] + "; " + lineage      #add rank of curr_id to beginning of lineage string
            
            if curr_id in nodes_dict:   #change curr_id to taxonID of parent node if parent node exists in nodes_dict

                curr_id = nodes_dict[curr_id]

            else:   #print error if parent node is not in nodes_dict

                print("Error: Can't find " + curr_id)

        gb_taxon_file.write(accession + "\t" + str(taxonId) + "\t" + lineage + "\n")    #write accession, taxonID of accession, and obtained lineage in three tab delimited columns

    gb_taxon_file.close()

################################### Make dictionary from taxonomy names.dmp file ###################################
#make_names_dict
#Usage: iterates through taxonomy names.dmp file and matches line for taxonID and scientific name of each node to make names_dmp dictionary
#Input: full path to taxonomy names.dmp file, format: 
# [taxon_id]|      [rank]   |           [...]           |   [category]  |
#   1       |       all     |                           |       synonym |
#   1       |       root    |                           |       scientific name |
#   2       |   Bacteria    |   Bacteria <bacteria>     |       scientific name |
#Returns: dictionary containing taxonID as key and scientific name as value
#
def make_names_dict(names_dmp):

    names_dict = {}     #initailize names_dict as key: node_id, value: scientifc name

    with open(names_dmp, "r") as names_file:                
    
        for line in names_file:

            if len(line) > 0 and "scientific name" in line:     #operates on lines labeled "scientific name"

                split_line = line.split("|")        #splits line using "|"
                taxon_id = split_line[0].strip()    #split_line[0] contains taxonID of node
                rank = split_line[1].strip()        #split_line[0] contains rank of node

                names_dict[taxon_id] = rank


    return names_dict

################################### Make dictionary from taxonomy nodes.dmp file ###################################
#make_nodes_dict
#Usage: iterates through taxonomy nodes.dmp file and matches line for taxonID of current node and taxonID of its parent node to make nodes_dmp dictionary
#Input: 
# - full path to taxonomy nodes.dmp file, format:
#  [child_id]	|	[parent_id] |     [rank]        | [...]
#	   2	    |	  131567	|	superkingdom	| [...]
#Returns: 
# - dictionary with key of taxonID of current node and value of taxonID of parent node
#
def make_nodes_dict(nodes_dmp):

    nodes_dict = {}     #intiailze nodes_dict as key: child node taxonID, value: parent node taxonID

    with open(nodes_dmp, "r") as nodes_file:                

        for line in nodes_file:

            if len(line) > 0:

                split_line = line.split("|")            #splits line using "|"
                child_node = split_line[0].strip()      #split_line[0] contains taxonID of child node
                parent_node = split_line[1].strip()     #split_line[1] contains taxonID of parent node
                    
                nodes_dict[child_node] = parent_node

    return nodes_dict

################################### Make dictionary of protein accessions to taxonIDs ###################################
#add_prot_taxIDs
#Usage: iterates through prot.accession2taxid and matches for accessions and taxonIDs for lines with accession.version in
# prot_accessions_hash
#Input: 
# - hash of desired sequence accession.version for protein database
# - full path to taxonomy prot.accession2taxid file containing accession, accession.version, taxid, and gi in four tab delimited columns, format: 
#     accession       accession.version       taxid       gi
#     A0A009IHW8         A0A009IHW8.1         1310613 1835922267
#Returns:
# - dictionary with keys as accessions of prot_accessions_hash and values of matching taxonIDs
#
def add_prot_taxIDs(prot_accessions_hash, pa2t_file):
                         
    pa2t = open (pa2t_file, "r")

    accession_version_to_taxon_dict = {}    #key: accession.version, value: taxonID

    for line in pa2t:

        split_line = line.split("\t")       #splits tab delimited lines

        if len(line) > 3:
                
            line_accession_version = split_line[1].strip()      #split_line[1] contains accession.version

            if line_accession_version in prot_accessions_hash:      #operates on pa2t_file lines with accession.version in the prot_accessions_hash

                accession_version_to_taxon_dict[line_accession_version] = split_line[2].strip()     #split_line[2] contains taxonID

    pa2t.close()

    return accession_version_to_taxon_dict

################################### Make dictionary of nucleotide accessions to taxonIDs ###################################
#add_nuc_taxIDs
#Usage: iterates through nuc_gb.accession2taxid and matches for accessions and taxonIDs for lines with accession in nuc_accessions_hash
#Input: 
# - hash of desired sequence accessions for nucleotide database
# - full path to taxonomy nuc_gb.accession2taxid file containing accession, accession.version, taxid, and gi in four tab delimited columns, format:
#	accession	accession.version	taxid	gi
#	  A00001	    A00001.1	    10641	58418
#Returns:
# - dictionary with keys of accession of nuc_accessions_hash and values of matching taxonIDs
#
def add_nuc_taxIDs(nuc_accessions_hash, na2t_file):

    na2t = open(na2t_file, "r")

    for line in na2t:   #key: accession, value: taxonID                               

        split_line = line.split("\t")       #split tab-delimited file           
        
        if len(line) > 3:
            
            line_accession = split_line[0]  #split_line[0] contains the accession number

            if line_accession in nuc_accessions_hash:

                nuc_accessions_hash[line_accession] = split_line[2]    #split_line[2] contains the taxonID 

    na2t.close()

    return nuc_accessions_hash

################################### Create a hash of plant virus accessions in fasta file ###################################
#hash_accessions
#Usage: iterates through fasta file and creates hash of accessions of sequences
#Input:
# - full path to fasta file containing sequences for database
#Returns:
# - dictionary with keys of accessions in fasta file and values all set to 1
#
def hash_accessions(fasta_path):

    fasta = open(fasta_path, "r")

    accession_dict = {}     #key: accession, value: 1

    for sequence in SeqIO.parse(fasta, "fasta"):    #iterate through fasta file and add accession to accession_dict

        seq_id = str(sequence.id)
        accession_dict[seq_id] = 1

    fasta.close()

    return accession_dict

################################### Annotate plant virus protein sequence fasta ###################################
#annotate_aa_fasta
#Usage: generates annotated file with headers of the following format:
#>acc|GENBANK|[protein_id]|GENBANK|[accession]|[annotation]
#
#Input:
# - full path to unannotated protein fasta file containing desired protein sequences for database
# - full path to tab delimited file containing accession in first column and protein id in second column, format:
#   [accession]      [protein_id]
#     AB000048        BAA19009.1
# - full path to desired location of annotated protein fasta file
# - dictionary with keys of accessions and values of annotations obtained from genbank.info file
#Generates:
# - protein fasta file containing annotations from genbank.info file
#
def annotate_aa_fasta(aa_input_file, id_file, prot_ann_path, genbank_dict):

    id_to_info = {}     #key: protein_id, value: annotation

    #annotate prot file
    viral_ids = gzip.open(id_file, "rt")

    for viral_line in viral_ids:    #iterate through viral_ids and make dict of protein_id to annotation
                
        viral_line_split = viral_line.split("\t")

        prot_id = viral_line_split[1].strip()       #viral_line_split[1] is the protein_id of the line

        accession = viral_line_split[0].strip()     #viral_line_split[0] is the accession of the line

        annotation = genbank_dict[accession]

        id_to_info[prot_id] = (accession, annotation)

    viral_ids.close()

    unannotated = open(aa_input_file, "r")
    with open(prot_ann_path, "w") as annotated:

        for sequence in SeqIO.parse(unannotated, "fasta"):  #iterate through ids in unannotated fasta and get annotation from id_to_ann

            prot_id = str(sequence.id)

            #id_to_info[prot_id][0] returns the accession of that protein_id, while id_to_info[prot_id][1] returns the annotation of that protein_id
            header = ">acc|GENBANK|" + prot_id + "|GENBANK|" + id_to_info[prot_id][0] + "|" + id_to_info[prot_id][1] + "\n"     #create annotated fasta header
            annotated.write(header)
            annotated.write(str(sequence.seq) + "\n")

    unannotated.close()

################################### Annotate plant virus nucleotide sequence fasta ###################################
#annotate_nuc_fasta
#Usage: 
#Input:
# - full path to unannotated nucleotide fasta file containing desired nucleotide sequences for database
# - full path to desired location of annotated nucleotide fasta file
# - dictionary with accessions as keys and annotations from genbank.info file as values
#Generates:
# - nucleotide fasta file with annotations from genbank.info
#
def annotate_nuc_fasta(nuc_input_file, nuc_ann_path, genbank_dict):

    new_fasta = open(nuc_ann_path, "w")                                                        
   
    with open(nuc_input_file, "r") as old_fasta:                                          
       
        for sequence in SeqIO.parse(old_fasta, "fasta"):    #match accession to genbank_dict for annotation                                   
           
            new_header = ">" + sequence.id + " " + genbank_dict[sequence.id.strip()] + "\n"        #create annotated fasta header
            new_fasta.write(new_header)                                                     
            new_fasta.write(str(sequence.seq) + "\n")                                      

    new_fasta.close()

################################### Make a dictionary of accessions to annotations using Genbank.info file ###################################
#make_genbank_dict
#Usage: iterates through Genbank.info file to generate dictionary containing accession and annotation of each line
#Input:
# - full path to Genbank.info file, format: 
#   [accession]   [isolation]     [annotation]                                                                                             [...]
#   AB000048        2007       protoparvovirus Feline panleukopenia virus gene for nonstructural protein 1, complete cds, isolate: 483.    [...]
#Returns:
# - dictionary with accessions as keys and annotations from genbank.info file as values
#
def make_genbank_dict(genbank_file):

    genbank_dict = {}   #key: accession, value: annotation                                            

    with gzip.open (genbank_file, "rt") as genbank:                        
                                              
        for line in genbank:      
            
            split_line = line.split("\t")                       
            accession = split_line[0].strip()   #split_line[0] is accession of each line
            annotation = split_line[3].strip()  #split_line[3] is the annotation of each line
            genbank_dict[accession] = annotation

    return genbank_dict

################################### Parse arguments ###################################
#parse_arguments
#Usage: parse through user provided arguments
#Returns: 
# - ArgumentParser object containing parameters passed to script
#
def parse_arguments():

    parser = argparse.ArgumentParser()

    #required arguments with full path to database location and taxonomy database files
    parser.add_argument("-db", dest="db_path", type=str, required=True, help="specify full path to database directory, like /path/to/plant_virus_db/")
    parser.add_argument("-n2", dest="na2t", type=str, required=True, help="specify full path to taxonomy nucl_gb.accession2taxid file, like /path/to/taxonomy/nucl_gb.accession2taxid")
    parser.add_argument("-p2", dest="pa2t", type=str, required=True, help="specify full path to taxonomy prot.accession2taxid file, like /path/to/taxonomy/prot.accession2taxid")
    parser.add_argument("-nd", dest="nodes", type=str, required=True, help="specify full path to taxonomy nodes.dmp file, like /path/to/taxonomy/nodes.dmp")
    parser.add_argument("-nm", dest="names", type=str, required=True, help="specify full path to taxonomy names.dmp file, like /path/to/taxonomy/names.dmp")

    #required arguments with full path to VirusDetect files  
    parser.add_argument("-ni", dest="nuc_input", type=str, required=True, help="unzip Plants_248_U95.tar.gz file and provide full path to fasta file from VirusDetect containing VirusDetect nucleotide database sequences, /path/to/vrl_Plants_248_U95")
    parser.add_argument("-a", dest="aa_input", type=str, required=True, help="unzip Plants_248_U95.tar.gz file and provide full path to fasta file from VirusDetect containing VirusDetect protein database sequences, like /path/to/vrl_Plants_248_U95_prot")
    parser.add_argument("-g", dest="genbank", type=str, required=True, help="specify full path to VirusDetect Genbank.info.gz file from VirusDetect for annotations, like /path/to/vrl_genbank.info.gz")
    parser.add_argument("-i", dest="id", type=str, required=True, help="specify full path to vrl_idmapping.gz file from VirusDetect containing ID mappings, like /path/to/vrl_idmapping.txt.gz")

    #optional arguments to name outputs
    parser.add_argument("--nuc_prefix", dest="nuc_prefix", required=False, type=str, help="specify desired prefix for nuc database files, default=\"vrl_Plants_248_U95\"", default="vrl_Plants_248_U95")
    parser.add_argument("--prot_prefix", dest="prot_prefix", required=False, type=str, help="specify desired prefix for prot database files, default=\"vrl_Plants_248_U95_prot\"", default="vrl_Plants_248_U95_prot")
   
    options = parser.parse_args()

    return options

################################### Main ###################################
def main():

    start = time.time()

    print("Parsing arguments")
    #parse arguments
    options = parse_arguments()

    #assign parsed arguments to variables
    full_path = options.db_path
    na2t_file = options.na2t
    pa2t_file = options.pa2t
    nodes = options.nodes
    names = options.names
    nuc_input_file = options.nuc_input
    aa_input_file = options.aa_input
    genbank_file = options.genbank
    id_file = options.id
    nuc_prefix = options.nuc_prefix
    prot_prefix = options.prot_prefix

    #build file names
    nuc_ann_path = full_path + nuc_prefix + "_annotated.fasta"
    nuc_gb_taxon_txt_file = full_path + nuc_prefix + "_annotated.gb_taxon.txt"
    prot_ann_path = full_path + prot_prefix + "_annotated.fasta"
    diamond_path = full_path + prot_prefix + "_annotated"
    prot_gb_taxon_txt_file = full_path + prot_prefix + "_annotated.gb_taxon.txt"

    #create annotated files
    print("Annotating fasta files")
    genbank_dict = make_genbank_dict(genbank_file)
    annotate_nuc_fasta(nuc_input_file, nuc_ann_path, genbank_dict)              #generates annotated nucleotide fasta file
    annotate_aa_fasta(aa_input_file, id_file, prot_ann_path, genbank_dict)      #generates annotated protein fasta file

    #prepare taxonomy files
    print("Preparing taxonomy files")
    nodes_dict = make_nodes_dict(nodes)     #make dictionary from nodes.dmp file
    names_dict = make_names_dict(names)     #make dictionary from names.dmp file

    #make nucleotide database output files
    print("Generating nucleotide database")
    subprocess.run("makeblastdb -in " + nuc_ann_path + " -dbtype nucl -title " + nuc_prefix , shell=True)
    nuc_accessions_hash = hash_accessions(nuc_ann_path)
    nuc_accession_to_taxon_dict = add_nuc_taxIDs(nuc_accessions_hash, na2t_file)
    make_gb_taxon_file(nuc_accession_to_taxon_dict, nodes_dict, names_dict, nuc_gb_taxon_txt_file)      #generates nuc gb_taxon file

    #make protein database output files
    print("Generating protein database")
    subprocess.run("diamond makedb --in " + prot_ann_path + " -d " + diamond_path, shell=True)          #generates prot dmnd file
    prot_accessions_hash = hash_accessions(aa_input_file)
    print("Protein accessions hash length: " + str(len(prot_accessions_hash)))
    prot_accession_to_taxon_dict = add_prot_taxIDs(prot_accessions_hash, pa2t_file)
    print("Prot_accession_to_taxon_dict length: " + str(len(prot_accession_to_taxon_dict)))
    make_gb_taxon_file(prot_accession_to_taxon_dict, nodes_dict, names_dict, prot_gb_taxon_txt_file)    #generates prot gb_taxon file

    end = time.time()
    print("Time to make database: " + str(end-start))

if __name__ == "__main__":
    main()
