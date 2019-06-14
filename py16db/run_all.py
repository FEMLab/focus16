#!/usr/bin/env python

import argparse
import sys
import os
import subprocess
import random
import shutil

import get_n_genomes as gng
import fetch_sraFind_data as fsd

def get_args():
    parser = argparse.ArgumentParser(
    description = "given a genus species, downloads the raw reads for each SRA using " +
    "fastq-dump. Also downloads downsampled [1000000 reads] version of each SRA")
    parser.add_argument("-o", "--output_dir", help="path to output", required=True)
    parser.add_argument("-n", "--organism_name", help="genus species in quotes", required=True)
    parser.add_argument("-s", "--sraFind_path", dest="sra_path",
                        default="srapure",
                        help="path to sraFind file with accession and genus/species", required=False)
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
    parser.add_argument("--get_all", help="number of strains",
                        actions="store_true",
                        required=True)
    return(parser.parse_args())


#srapure = file containing only [0]accession and [1]genuspecies for 16S sequences
#remove " and ' from the lines

def filter_srapure(path, organism_name, strains, get_all):
    results = []
    with open(path, "r") as infile:
        for line in infile:
            split_line = [x.replace('"', '').replace("'", "") for x in line.strip().split("\t")]
            if split_line[11] == organism_name:
                results.append(split_line[17])
    random.shuffle(results)
    if strains != 0:
        results = results[0:strains]

    sras = []
    for result in results:
        these_sras = results.split(",")
        if get_all:
            for sra in these_sras:
                sras.append(sra)
        else:
            sras.append(these_sras[0])
    return(sras)

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


def pob(genomes_dir, readsf, output_dir):
    pobcmd = "plentyofbugs -g " + genomes_dir +  " -f " + readsf + " -o " + output_dir
    subprocess.run(pobcmd,
                   shell=sys.platform !="win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
    best_ref = os.path.join(output_dir, "best_reference")
    with open(best_ref, "r") as infile:
        for line in infile:
            print(line)
            sraacc = line.strip().split('\t')
    return(sraacc)


if __name__ == '__main__':
    # main()
    #calling the previously made functions
    args=get_args()
    os.makedirs(args.output_dir)
    if not os.path.exists(args.sra_path):
        sraFind_output_dir = os.path.join(args.output_dir, "sraFind")
        args.sra_path = fsd.main(args, output_dir=sraFind_output_dir)
    filtered_sras = filter_srapure(
        path=args.sra_path,
        organism_name=args.organism_name,
        strains=args.nstrains,
        get_all=args.get_all)
    if os.path.exists(args.genomes_dir):
        if len(os.listdir(args.genomes_dir)) == 0:
            print("Warning: genome directory exists but is " +
                  "empty: downloading genomes")
            shutil.rmtree(args.genomes_dir)
            gng.main(args)
        else:
            pass
    else:
        gng.main(args)
    for i, accession in enumerate(filtered_sras):
        this_output=os.path.join(args.output_dir, str(i))
        os.makedirs(this_output)
        print("Downloading " + accession)
        download_SRA(SRA=accession, destination=this_output)
        readsf = os.path.join(this_output, "/downsampled/reads1.fq")
        pob_dir = os.path.join(output_dir, "plentyofbugs")
        pob(genomes_dir=genomes_dir, readsf=readsf, output_dir=pob_dir)
