#! /usr/bin/env python
######################################################################
#extract_kraken_reads.py takes in a kraken-style output and kraken report
#and a taxonomy level to extract reads matching that level
#Copyright (C) 2019 Jennifer Lu, jlu26@jhmi.edu
#
#This file is part of KrakenTools
#KrakenTools is free software; oyu can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 3 of the license, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of 
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program; if not, see <http://www.gnu.org/licenses/>.
#
######################################################################
#Jennifer Lu, jlu26@jhmi.edu
#Updated: 06/03/2019
#
#This program extracts reads classified by Kraken as a 
#specified taxonomy ID. Those reads are extracted into a new FASTA file.
#
#Required Parameters:
#   -k, --kraken, --kraken-file X.......kraken output file
#   -s, -s1, -1, -U X...................read file 
#                                       [FASTA/FASTQ - may be gzipped]
#   -s2, -2, X..........................second read file if paired 
#                                       [FASTA/FASTQ - may be gzipped]
#   -o, --output X......................output FASTA file with reads 
#   -t, --taxid, --taxids X.............list of taxonomy IDs to extract 
#                                       [separated by spaces]
#   -r, --report-file X.................kraken report file
#                                       [required only with --include-children/parents]
#Optional Parameters:
#   -h, --help..........................show help message.
#   --max X.............................only save the first X reads found
#   --include-children **...............include reads classified at lower levels 
#   --include-parents **................include reads classified at parent levels 
#                                       of taxids 
#   --append............................append extracted reads to output file if existing
#   --noappend..........................rewrite file if existing [default] 
# ** by default, only reads classified exactly at taxids provided will be extracted
# ** if either of these are specified, a report file must also be provided 
######################################################################
import os, sys, argparse
from time import gmtime
from time import strftime
import re
#################################################################################
#Tree Class 
#usage: tree node used in constructing taxonomy tree  
#   includes only taxonomy levels and genomes identified in the Kraken report
class Tree(object):
    'Tree node.'
    def __init__(self, taxid, level_num, level_id, children=None, parent=None):
        self.taxid = taxid
        self.level_num = level_num
        self.level_id = level_id
        self.children = []
        self.parent = parent
        if children is not None:
            for child in children:
                self.add_child(child)
    def add_child(self, node):
        assert isinstance(node,Tree)
        self.children.append(node)
#################################################################################
#process_kraken_output
#usage: parses single line from kraken output and returns taxonomy ID and readID
#input: kraken output file with readid and taxid in the
#   second and third tab-delimited columns
#returns: 
#   - taxonomy ID
#   - read ID
def process_kraken_output(kraken_line):
    l_vals = kraken_line.split('\t')
    searchObj = re.search( r'\(taxid (\d+)\)', l_vals[2])  #modified to fit Kraken2 --use-name output by Xiaojun
    taxid = 0
    if searchObj:
        tax_id = searchObj.group(1) #return str
    read_id = l_vals[1]
    #if (tax_id == 'A'):
    #    tax_id = 81077
    #else:
    #    tax_id = int(tax_id)
    return [tax_id, read_id]

#process_kraken_report
#usage: parses single line from report output and returns taxID, levelID
#input: kraken report file with the following tab delimited lines
#   - percent of total reads
#   - number of reads (including at lower levels)
#   - number of reads (only at this level)
#   - taxonomy classification of level
#       (U, - (root), - (cellular org), D, P, C, O, F, G, S)
#   - taxonomy ID (0 = unclassified, 1 = root, 2 = Bacteria...etc)
#   - spaces + name
#returns:
#   - taxonomy ID 
#   - level number (number of spaces before name)
#   - level_type (type of taxonomy level - U, R, D, P, C, O, F, G, S, etc) 
def process_kraken_report(report_line):
    l_vals = report_line.strip().split('\t')
    try:
        int(l_vals[1])
    except ValueError:
        return []
    #Extract relevant information
    level_type = l_vals[3]
    taxid = l_vals[4]
    #Get spaces to determine level num
    spaces = 0
    for char in l_vals[-1]:
        if char == ' ':
            spaces += 1
        else:
            break
    level_num = int(spaces/2)
    return[taxid, level_num, level_type]
#################################################################################
#Main method 
def main():
    #Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-k', dest='kraken_file', required=True,
        help='Kraken output file to parse')
    parser.add_argument('-t', dest='taxid', required=True,
        type=int, nargs='+',
        help='Taxonomy ID[s] of reads to extract (space-delimited)')
    parser.add_argument('-o', dest='output_file', required=True,
        help='Output read IDs, one id per line')
    parser.add_argument('--max', dest='max_reads', required=False, 
        default=100000000, type=int,
        help='Maximum number of reads to save [default: 100,000,000]')
    parser.add_argument('-r','--report',dest='report_file', required=False,
        default="",
        help='Kraken report file. [required only if --include-parents/children \
        is specified]')
    parser.add_argument('--include-parents',dest="parents", required=False, 
        action='store_true',default=False,
        help='Include reads classified at parent levels of the specified taxids')
    parser.add_argument('--include-children',dest='children', required=False,
        action='store_true',default=False,
        help='Include reads classified more specifically than the specified taxids')
    parser.add_argument('--exclude', dest='exclude', required=False,
        action='store_true',default=False,
        help='Instead of finding reads matching specified taxids, finds all reads NOT matching specified taxids') 
    parser.set_defaults(append=False)

    args=parser.parse_args()
    
    #Initialize taxids
    save_taxids = {}
    for tid in args.taxid:
        save_taxids[str(tid)] = 0  #force str by Xiaojun
    #print("initial taxid: ", save_taxids) 
    main_lvls = ['R','K','D','P','C','O','F','G','S']

    #STEP 0: READ IN REPORT FILE AND GET ALL TAXIDS 
    if args.parents or args.children:
        #check that report file exists
        if args.report_file == "": 
            sys.stderr.write(">> ERROR: --report not specified.")
            exit(1)
        #sys.stdout.write(">> STEP 0: PARSING REPORT FILE %s\n" % args.report_file)
        #create tree and save nodes with taxids in the list 
        base_nodes = {} 
        r_file = open(args.report_file,'r')
        prev_node = -1
        level_id = 'U' #initalize level_id by Xiaojun
        for line in r_file:
            #extract values
            report_vals = process_kraken_report(line)
            if len(report_vals) == 0:
                continue
            [taxid, level_num, level_type] = report_vals
            #print(taxid, level_num, level_type)
            #tree root
            #print("taxid ", taxid)
            if(taxid == '0'): #skip unclissified row, fixed by Xiaojun
                continue
            if taxid == '1':
                root_node = Tree(taxid, level_num, level_type)
                prev_node = root_node
                #save if needed
                if taxid in save_taxids:
                    base_nodes[taxid] = root_node
                continue
            #move to correct parent
            while level_num != (prev_node.level_num + 1):
                prev_node = prev_node.parent 
            #determine correct level ID 
            if level_id == '-' or len(level_id) > 1:
                if prev_node.level_id in main_lvls:
                    level_id = prev_node.level_id + '1'
                else:
                    num = int(prev_node.level_id[-1]) + 1
                    level_id = prev_node.level_id[:-1] + str(num)
            #make node
            curr_node = Tree(taxid, level_num, level_type, None, prev_node)
            prev_node.add_child(curr_node)
            prev_node = curr_node
            #save if taxid matches
            #print("taxid ",taxid, type(taxid))
            if taxid in save_taxids:
                base_nodes[taxid] = curr_node 
        r_file.close()
        #FOR SAVING PARENTS
        if args.parents:
            #For each node saved, traverse up the tree and save each taxid 
            for tid in base_nodes:
                curr_node = base_nodes[tid]
                while curr_node.parent != None:
                    curr_node = curr_node.parent
                    save_taxids[curr_node.taxid] = 0
        #FOR SAVING CHILDREN 
        if args.children:
            for tid in base_nodes:
                curr_nodes = base_nodes[tid].children
                while len(curr_nodes) > 0:
                    #For this node
                    curr_n = curr_nodes.pop()
                    if curr_n.taxid not in save_taxids:
                        save_taxids[curr_n.taxid] = 0
                    #Add all children
                    if curr_n.children != None:
                        for child in curr_n.children:
                            curr_nodes.append(child)

    ##############################################################################
    sys.stdout.write("process file %s\n" % args.kraken_file)
    sys.stdout.write("%i taxonomy IDs to parse\n" % len(save_taxids))
    #sys.stdout.write(">> STEP 1: PARSING KRAKEN FILE FOR READIDS %s\n" % args.kraken_file)
    #Initialize values
    count_kraken = 0
    read_line = -1
    
    #PROCESS KRAKEN FILE FOR CLASSIFIED READ IDS
    k_file = open(args.kraken_file, 'r')
    #sys.stdout.write('\t0 reads processed') #change by Xiaojun
    #sys.stdout.flush() #change by Xiaojun
    #Evaluate each sample in the kraken file
    save_readids = {}
    for line in k_file:
        count_kraken += 1
        #if (count_kraken % 10000 == 0): #stop output number by Xiaojun
        #    sys.stdout.write('\r\t%0.2f million reads processed' % float(count_kraken/1000000.))
        #    sys.stdout.flush()
        #Parse line for results
        [tax_id, read_id] = process_kraken_output(line)
        
        #Skip if reads are human/artificial/synthetic
        if (tax_id in save_taxids) and not args.exclude:
            save_taxids[tax_id] += 1
            save_readids[read_id] = 0 
        elif (taxid not in save_taxids) and args.exclude:
            save_taxids[tax_id] += 1
            save_readids[read_id] = 0 
        if len(save_readids) >= args.max_reads:
            break 
    #Update user
    k_file.close()
    sys.stdout.write('t%0.2f million reads processed\n' % float(count_kraken/1000000.))
    sys.stdout.write('%i read IDs saved\n' % len(save_readids))
    #print(save_readids.keys())  
    #output read id to a file
    #print("\n".join("{0}".format(id) for id in save_readids))  #{!r}, {0} format without quote
    with open(args.output_file, 'w') as fout:
        for id in save_readids:
            print("{0}".format(id), file=fout)

    #End of program
    #time = strftime("%m-%d-%Y %H:%M:%S", gmtime())
    #sys.stdout.write("PROGRAM END TIME: " + time + '\n')
    exit(0)

#################################################################################

if __name__ == "__main__":
    main()

