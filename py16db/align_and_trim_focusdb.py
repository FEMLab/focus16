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
    parser.add_argument(
        "-i", "--input",
        help="Multifasta input file", required=True)
    parser.add_argument(
        "-o", "--out_prefix",
        help="Prefix for msa and trimmed msa", required=True)
    parser.add_argument(
        "--addcoli", action="store_true",
        help="Add cannonical MG1655 E coli sequence")
    return(parser.parse_args())


def add_cannonical(msa, tmp):
    coli_header = str(">NC_000913.3:4166659-4168200 Escherichia coli " +
                      "str. K-12 substr. MG1655, complete genome")
    coli_sequence = \
"""
AAATTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGT
AACAGGAAGAAGCTTGCTTCTTTGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGAT
GGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCTTC
GGGCCTCTTGCCATCGGATGTGCCCAGATGGGATTAGCTAGTAGGTGGGGTAACGGCTCACCTAGGCGAC
GATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAG
GCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCT
TCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCC
GCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAA
TTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAA
CTGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCG
TAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGC
GTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCGACTTGGAGGTTGTGCC
CTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTCGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAAC
TCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACC
TTACCTGGTCTTGACATCCACGGAAGTTTTCAGAGATGAGAATGTGCCTTCGGGAACCGTGAGACAGGTG
CTGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCC
TTTGTTGCCAGCGGTCCGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATG
ACGTCAAGTCATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCG
ACCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCAT
GAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACC
GCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCACT
TTGTGATTCATGACTGGGGTGAAGTCGTAACAAGGTAACCGTAGGGGAACCTGCGGTTGGATCACCTCCT
TA
"""
    # table from  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4802574/
    regions = {
        "V1": [8,    96],
        "V2": [97,   306],
        "V3": [307,  487],
        "V4": [488,  746],
        "V5": [747,  885],
        "V6": [886,  1029],
        "V7": [1030, 1180],
        "V8": [1181, 1372],
        "V9": [1373, 1468],
        "VT": [8,    1468]
    }
    with open(msa, "r") as inf, open(tmp,"w") as  outf:
        outf.write("%s\n%s\n" % (coli_header, coli_sequence.replace('\n', '').strip()))
        for line in inf:
            outf.write(line)
    return tmp


def mafft(pre, multifasta):
    ''' performs default mafft alignment
    '''
    msa = pre + ".mafft"
    cmd = "mafft --retree 2 --reorder {multifasta} > {msa}".format(**locals())
    return(msa, cmd)


def run_TrimAl(msa):
    ''' performs default mafft alignment
    '''
    outmsa = msa + ".trimmed"
    outhtml = msa + ".trimmed.html"
    cmd = str("trimal -in {msa} -htmlout {outhtml} -out {outmsa} " +
              # "-gappyout" +
              "-strictplus -keepheader").format(**locals())
    print(cmd)
    subprocess.run(cmd,
                   shell=sys.platform !="win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
    print("Trimming complete")
    return(outmsa)


def main():
    args = get_args()
    if shutil.which("trimal") is None:
        print("This script requires TrimAl, which can be " +
              "installed from https://github.com/scapella/trimal/releases" +
              " . Exiting...")
        sys.exit(1)
    if args.out_prefix.endswith(os.path.sep):
        print("Output prefix should not end with path separator")
        sys.exit(1)
    if args.addcoli:
        multifasta = add_cannonical(
            msa=args.input,
            tmp=args.out_prefix  + "with_coli")
    else:
        multifasta = args.input

    msa, msa_cmd = mafft(pre=args.output_pre, multifasta = multifasta)
    print(msa_cmd)
    subprocess.run(msa_cmd,
                   shell=sys.platform !="win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
    print("MSA complete")

    run_TrimAl(msa)


if __name__ == '__main__':
    main()
