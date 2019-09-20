#!/usr/bin/env python3

import argparse
import subprocess
import sys
import shutil
import os
import subprocess
import math
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
from collections import Counter


def get_args():
    parser = argparse.ArgumentParser(
        description="Given a multiple sequence database (from " +
        "combine-focusdb-and-silva, generate an alignment with " +
        "mafft and trim to median sequence.  Requires mafft and TrimAl" )
    # REQUIRED
    parser.add_argument("-i", "--input", help="multifasta input file", required=True)
    parser.add_argument("-o", "--out_prefix", help="prefix for msa and trimmed msa", required=True)
    return(parser.parse_args())


def mafft(args, multifasta):
    ''' performs default mafft alignment
    '''
    msa = args.out_prefix + ".mafft"
    cmd = "mafft --retree 2 --reorder {multifasta} > {msa}".format(**locals())
    print(cmd)
    subprocess.run(cmd,
                   shell=sys.platform !="win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
    print("MSA complete")
    return(msa)


def run_TrimAl(msa):
    ''' performs default mafft alignment
    '''
    outmsa = msa + ".trimmed"
    cmd = "trimal -in {msa} -out {outmsa} -gappyout".format(**locals())
    print(cmd)
    subprocess.run(cmd,
                   shell=sys.platform !="win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
    print("Trimming complete")
    return(msa)



def main():
    args = get_args()
    if shutil.which("trimal") is None:
        print("This script requires TrimAl, which can be " +
              "installed from https://github.com/scapella/trimal/releases" +
              " . Exiting...")
        sys.exit(1)

    msa = mafft(args=args, multifasta = args.input)
    run_TrimAl(msa)


if __name__ == '__main__':
    main()
