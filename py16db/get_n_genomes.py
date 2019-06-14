#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Nick Waters
"""

import os
import sys
import subprocess
import argparse
import shutil
from random import shuffle

def get_args():  # pragma nocover
    parser = argparse.ArgumentParser(
        description="Blast assembly against core genome to find and " +
        "eliminate truncated genes due to bad assembly, " +
        "returning a assembly  with genes that meet a set of criteria. ",
        add_help=True)
    parser.add_argument(
        "-o",
        "--organism_name",
        help="organism_name",
        required=True)
    parser.add_argument(
        "-g",
        "--genomes_dir",
        help="path to output dir",
        required=True)
    parser.add_argument(
        "-n", "--nstrains",
        help="",
        type=int,
        required=True)

    parser.add_argument(
        "-p",
        "--prokaryotes",
        action="store",
        help="path_to_prokaryotes.txt",
        required=True)
    return(parser.parse_args())



def main(args=None):
    if args is None:
        args = get_args()
    try:
        os.makedirs(args.genomes_dir)
    except:
        print("genomes output directory already  exists!")
        sys.exit(1)

    if not os.path.exists(args.prokaryotes):

        subprocess.run("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt -O " + args.prokaryotes,
                       shell=sys.platform != "win32",
                       check=True)

    # column 9 has the nucc accession if it is a complete genome, or a "-" if empt
    # it starts with chromasome
    # here we select the lines for all the complete genomes with awk,
    # find the lines matching out query
    # and save a temp file with the results
    tmp_file = "/tmp/prok_subset_raw_outfile"
    if os.path.exists(tmp_file):
        os.remove(tmp_file)
    org_lines = []
    with open(args.prokaryotes, "r") as proks:
        for line in proks:
            splitline = line.strip().split("\t")
            if splitline[8].startswith("chrom"):
                if splitline[0].startswith(args.organism_name):
                    org_lines.append(splitline)
    # for line in org_lines:
    #    print(line[20])
    # # if file is empty, raise an error
    if len(org_lines) == 0:
        print("no " + args.organism_name +
              " matches in the prokaryotes.txt file")

    shuffle(org_lines)
    # # now we shuffle the file, get the top n lines, and use some sed to split apart the
    # # chromosome:NZ_CP013218.1/CP013218.1; plasmid paadD:NZ_CP014695.1/CP014695.1; plasmid plinE154:NZ_CP014694.1/CP014694.1
    # # to
    # # NZ_CP013218.1
    # # Note that we only get the first chromasome for a given entry. Sorry vibrioists
    for line in org_lines:
        print(line[8].split(":")[1].split(";")[0].split("/")[0])

    for line in org_lines[0:args.nstrains]:
        shortname = line[8].split(":")[1].split(";")[0].split("/")[0]
        name = os.path.basename(line[20])
        full_path = os.path.join(
            line[20],
            name + "_genomic.fna.gz")
        cmd = "wget " + full_path + " -O " + os.path.join(
            args.genomes_dir, shortname + ".fna.gz")
        print(cmd)
        subprocess.run(
            cmd,
            shell=sys.platform != "win32",
            check=True)

    unzip_cmd = "gunzip " + os.path.join(args.genomes_dir, "*")
    print(unzip_cmd)
    subprocess.run(
        unzip_cmd,
        shell=sys.platform != "win32",
        check=True)

if __name__ ==  "__main__":
    args = get_args()
    main(args)
