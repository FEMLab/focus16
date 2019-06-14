#!/usr/bin/env python

import argparse
import sys
import os
import subprocess
import random

import get_n_genomes as gng
import fetch_sraFind_data as fsd

def get_args():
    parser = argparse.ArgumentParser(
    description = "given a genus species, downloads the raw reads for each SRA using " +
    "fastq-dump. Also downloads downsampled [1000000 reads] version of each SRA")
    parser.add_argument("-o", "--output_dir", help="path to output", required=True)
    parser.add_argument("-n", "--organism_name", help="genus species in quotes", required=True)
    parser.add_argument("-s", "--sraFind_path", help="path to sraFind file with accession and genus/species", required=False)
    parser.add_argument(
        "-g",
        "--genomes_dir",
        help="path to genomes directory",
        required=True)
    parser.add_argument(
        "-p",
        "--prokaryotes",
        action="store",
        help="path_to_prokaryotes.txt",
        default="./prokaryotes.txt",
        required=False)
    parser.add_argument("-S", "--nstrains", help="number of strains",
                        type=int, required=True)
    return(parser.parse_args())


#srapure = file containing only [0]accession and [1]genuspecies for 16S sequences
#remove " and ' from the lines

def filter_srapure(path, organism_name, strains):
    results = []
    with open(path, "r") as infile:
        for line in infile:
            split_line = [x.replace('"', '').replace("'", "") for x in line.strip().split("\t")]
            if split_line[1] == organism_name:
                results.append(split_line[0])
    random.shuffle(results)
    if strains == 0:
       return(results)
    else:
        return(results[0:strains])

#create two suboutputs, raw and downampled
#create a variable that is the fastq-dump command, with relevant SRA and output

def download_SRA(SRA, destination):
    suboutput_dir_raw = os.path.join(destination, "raw")
    suboutput_dir_downsampled = os.path.join(destination, "downsampled")
    os.makedirs(suboutput_dir_raw)
    os.makedirs(suboutput_dir_downsampled)
    cmd = "fastq-dump --split-files " + SRA + " -O " + suboutput_dir_raw
    subprocess.run(cmd,
                   shell=sys.platform !="win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
    downcmd = "seqtk sample -s100 " + suboutput_dir_raw + "/" + SRA + "_1.fastq 1000000 > " + suboutput_dir_down + "reads" + SRA + ".fq"
    subprocess.run(downcmd,
                  shell=sys.platform !="win32",
                  stdout=subprocess.PIPE,
                  sterr=subprocess.PIPE,
                   check=True)
    return()


def pob(genomes_dir, readsf):
    suboutput_dir = os.path.join(output_dir, "plentyofbugs")
    os.makedirs(suboutput_dir)
    pobcmd = "plentyofbugs -g" + genomes_dir +  "-f" + readsf
    subprocess.run(pobcmd,
                   shell=sys.platform !="win32",
                   stdout=subprocess.PIPE,
                   sterr=subprocess.PIPE,
                   check=True)


if __name__ == '__main__':
    # main()
    #calling the previously made functions
    args=get_args()
    os.makedirs(args.output_dir)
    if not os.path.exists(args.sra_path):
        args.sra_path = fsd.main(args)
    filtered_sras = filter_srapure(
        path=args.sra_path,
        organism_name=args.organism_name,
        strains=args.nstrains)
    gng.main(args)
    for i, accession in enumerate(filtered_sras):
        this_output=os.path.join(args.output_dir, str(i))
        os.makedirs(this_output)
        print("Downloading " + accession)
        download_SRA(SRA=accession, destination=this_output)
        readsf = os.path.join(this_output, "/downsampled/reads1.fq")
        pob(genomes_dir=genomes_dir, readsf=readsf)
