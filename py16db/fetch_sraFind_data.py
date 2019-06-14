#!/usr/bin/env python

import argparse
import sys
import os
import subprocess
import random

def get_args():
    parser = argparse.ArgumentParser(
    description = "given a genus species, downloads the raw reads for each SRA using " +
    "fastq-dump. Also downloads downsampled [1000000 reads] version of each SRA")
    parser.add_argument("-o", "--output_dir", help="path to output", required=True)
    parser.add_argument(
        "-p",
        "--prokaryotes",
        action="store",
        help="path_to_prokaryotes.txt",
        default="./prokaryotes.txt",
        required=False)
    return(parser.parse_args())


def main(args, output_dir):
    # download the file from github
    if os.path.exists(output_dir):
        print("Output folder already exists; exiting...")
        sys.exit(1)
    else:
        os.makedirs(output_dir)
    sraFind_results = "https://raw.githubusercontent.com/nickp60/sraFind/master/results/sraFind-CompleteGenome-biosample-with-SRA-hits.txt"
    # gets just the file name
    sraFind_results_path = os.path.join(
        output_dir, os.path.basename(sraFind_results))
    if not os.path.exists(sraFind_results_path):
        print("Downloading sraFind Dump")
        download_sraFind_cmd = str("wget " + sraFind_results + " -O " + sraFind_results_path)
        for cmd in [download_sraFind_cmd]:
            subprocess.run(
                download_sraFind_cmd,
                shell=sys.platform !="win32",
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True)
    return sraFind_results_path


if __name__ == "__main__":
    args = get_args()
    main(args, output_dir=args.output_dir)
