#!/usr/bin/env python3

import argparse
import sys
import os
import math
from Bio import SeqIO


def get_args():
    parser = argparse.ArgumentParser(
        description="Given a trimmed multiple sequence alignment (from " +
        "align-and-trim-focusdb, calculate shannon entropy " )
    # REQUIRED
    parser.add_argument("-i", "--input", help="trimmed MSA", required=True)
    return(parser.parse_args())


def shannon_calc(values):
    """ calculate shannon entropy for a position
    """

    possible = set(values.split())
    entropy = 0
    for nuc in ["a", "t", "c", "g", "-"]:
        n = values.count(nuc)
        if n == 0:
            continue
        prob = n/len(values)
        ent = prob * math.log(prob)
        entropy = entropy + ent

    return(entropy)


def read_in_msa(path):
    """
    count occurances of each nucleotide in a sequence
    """
    seqs = []
    leaddash = []
    traildash = []
    with open(path, "r") as inf:
        for rec in SeqIO.parse(inf, "fasta"):
            this_lead_dash = 0
            this_trail_dash = 0
            for n in rec.seq:
                if n == "-":
                    this_lead_dash += 1

                else:
                    break
            leaddash.append(this_lead_dash)
            seqs.append(str(rec.seq).lower())
    return seqs


def main():
    args = get_args()
    seqs = read_in_msa(args.input)
    for i in range(len(seqs[0])):
        values_string = "".join([x[i] for x in seqs])
        entropy = shannon_calc(values_string)
        sys.stdout.write("{i}\t{entropy}\n".format(**locals()))
