#!/usr/bin/env python

import argparse
import sys
import os
import subprocess




   
#srapure = file containing only [0]accession and [1]genuspecies for 16S sequences
#remove " and ' from the lines

def filter_srapure(path, organism_name):
    results = []
    with open(path, "r") as infile:
        for line in infile:
            split_line = [x.replace('"', '').replace("'", "") for x in line.strip().split("\t")]
            print(split_line)
            if split_line[1] == organism_name:
                results.append(split_line[0])
    return(results)

#create two suboutputs, raw and downampled
#create a variable that is the fastq-dump command, with relevant SRA and output

def download_SRA(SRA, destination):
    suboutput_dir_raw = os.path.join(destination, "raw")
    suboutput_dir_downsampled = os.path.join(destination, "downsampled")
    os.mkdir(suboutput_dir_raw)
    os.mkdir(suboutput_dir_downsampled)
    cmd = "fastq-dump --split-files " + SRA + " -O " + suboutput_dir_raw
    subprocess.run(cmd,
                   shell=sys.platform !="win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
    downcmd = "seqtk sample -s100 " + suboutput_dir_raw +  "1000000 > " + suboutput_dir_down + "reads" + SRA + ".fq"
    subprocess.run(downcmd,
                  shell=sys.platform !="win32",
                  stdout=subprocess.PIPE,
                  sterr=subprocess.PIPE,
                   check=True)
    return()

    
#argument function

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output_dir", help="path to output", required=True)
    parser.add_argument("-n", "--organism_name", help="genus species in quotes", required=True)
    parser.add_argument("-s", "--sra_path", help="path to sra file with accession and genus/species", required=True)
    return(parser.parse_args())

if __name__ == '__main__':
# main()
#calling the previously made functions
 args=get_args()
 makedirs(args.output)
 filtered_sras = filter_srapure(path=args.srapure, organism_name=args.organism)
 for i, accession in enumerate(filtered_sras):
     this_ouput=os.path.join(args.output, i)
     os.mkdir(this_output)
     download_SRA(SRA=accession, destination=this_output)







