#!/usr/bin/env python 

import argparse
import os
import sys

#argument function
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-v",
        "--working_directory",
        help="current working directory",
        required=True)
    parser.add_argument(
    "-r",
    "--input_reads",
    help="input file for downsampled forward reads for draft genomes",
    required=True)   
    parser.add_argument(
        "-o",
        "--organism_name",
        help="bug of interest",
        required=True)
    parser.add_argument(
        "-n",
        "--number of strains",
        help="number of genomes that will be compared as the reference",
        required=True)
    parser.add_argument(
        "-e",
        "--experiment_name",
        help="name of experiment",
        required=True)
    parser.add_argument(
        "-d",
        "--output",
        help="output folder",
        required=True)
    return(parser.parse_args())

#function to form command to use docker container plentyofbugs
def pob(wd, readfile, organism, strains, experiment, output):
    pobcmd =
    "docker run --rm -t -v "
    + wd +
    ":/data/ nickp60/plentyofbugs:0.87  -f"
    + readfile +
    " -o "
    + organism +
    " -n "
    + strains +
    "-e "
    + experiment +
    " -d "
    +  output 
    return()

#function to form command to use get_genomes.py. Downloads genomes, in this case, the sra chosen by pob
def get_genomes(ref, output):
    ggcmd = "get_genomes.py -q " + ref + " -o" + output 
    return()
 
#calling the above functions to obtain best reference genome and download it 
if __name__ == '__main__':
 args=get_args()
 inputreads = args.input_reads
 for read in inputreads: 
     best_ref = pob(wd=args.directory, readfile=read, organism=args.organism, strains=args.maxstrain, experiment=name, output=poboutput)
 for sra in best_ref(poboutput):
     os.mkdir(ggoutput)
     get_genomes(ref=sra, output=ggoutput)
 
 #download the sra using get_genomes():
 
 
 
    
