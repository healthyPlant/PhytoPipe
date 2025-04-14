#! /usr/bin/env python
######################################################################
#Alex Hu <xiaojun.hu@usda.gov>
#Updated: 02/18/2021
#This program draws a coverage graph for a virus genome
#
#Required Parameters:
#   -i, --input X..........................coverage tsv file 
#   -o, --output X.........................coverage grpah 
######################################################################

import sys, argparse
import numpy as np
import pandas as pd
from collections import OrderedDict, defaultdict
import csv
import math

#to solve error: couldn't connect to display "localhost"
import matplotlib
try:
    matplotlib.use("Agg")
except ValueError:
    pass
import matplotlib.pyplot as plt

#################################### Plot coverage function ####################################
def plot_pileup_coverage(tabFile, graphFile):
    """
    plot coverage graph based on samtools depth/pileup file
    """
    title = tabFile.replace(".txt","")
    #title = tabFile.replace(".pileup.tab","")
    table = pd.read_csv(tabFile, sep='\t', header=None)
    #print(table)
    alignment = pd.DataFrame(data=table)
    refName = alignment.iloc[0,0]
    position = alignment.iloc[:,1]
    depth = alignment.iloc[:,2]
    max_dep = max(depth)
    min_dep = min(depth)
    #mean_dep = sum(depth)/len(depth)
    #mean_dep = np.mean(depth)
    median_dep = np.median(depth)

    #plt.plot(position, depth, 'bo', markersize=1)
    plt.stackplot(position,depth,colors='cyan',edgecolor='black')
    plt.xlabel('Position in Genome ' + refName)
    plt.ylabel('Depth of Coverage')
    plt.ylim(min_dep+0.5, max_dep*1.10)
    #add mean value
    plt.axhline(median_dep, color='red', alpha=0.5, linestyle='--', linewidth=1)
    plt.legend(('median: {:5.0f}'.format(median_dep),'depth'),loc='upper right')
    plt.title(title)
    #plt.show()

    plt.gcf().set_size_inches(10,6)
    plt.savefig("%s" % graphFile)

#################################### Get arguments function ####################################
def parseArguments():
    """
    Get arguments
    """
    parser = argparse.ArgumentParser(description = "Draw coverage graph.")
    parser.add_argument('-i','--input',dest='input_file', required=True, default="", help='Coverage tsv file')
    parser.add_argument('-o','--output',dest='output_file', required=True, default="", help='Coverage grpah')
    parser.set_defaults(append=False)
    return parser.parse_args()

def main():
    ### Input arguments
    options = parseArguments()
    tsvFile = options.input_file
    graphFile = options.output_file

    #call download function
    plot_pileup_coverage(tsvFile, graphFile)

if __name__ == '__main__':
	main()
