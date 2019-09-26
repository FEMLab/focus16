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
    parser.add_argument(
        "-i", "--input",
        help="trimmed MSA",
        required=True)
    parser.add_argument(
        "--colifirst",
        help="use first for coli coords",
        action="store_true")
    parser.add_argument(
        "--pattern",
        help="pattern for grepping groups from seqeunce name",
        type=str)
    # parser.add_argument(
    #     "--region",
    #     help="Just analyse this hypervariable region",
    #     choices=["all", "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9"],
    #     default="all")
    args = parser.parse_args()
    # if args.region != "all" and not args.colifirst:
    #     print("Warning: cannot extract from region without having " +
    #           "ecoli first in alignment, using --colifirst")
    return(args)


def detect_ecoli_positions(path):
    # table from  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4802574/
    regions = {
        # "all": [8,    1468],  # VT
        1:  [8,    96],
        2:  [97,   306],
        3:  [307,  487],
        4:  [488,  746],
        5:  [747,  885],
        6:  [886,  1029],
        7:  [1030, 1180],
        8:  [1181, 1372],
        9:  [1373, 1468]
    }

    #this27F = "AGAGTTTGATCATGGCTCAG".lower()
    this357F =  "CTCCTACGGGAGGCAGCAG".lower()
    with open(path, "r") as inf:
        for i, rec in enumerate(SeqIO.parse(inf, "fasta")):
            if i == 0:
                aligned_indexes = []
                chars = ""
                for i, char in enumerate(list(str(rec.seq).lower())):
                    if char != "-":
                        chars = chars + char
                        aligned_indexes.append(i)
                # get where 357 is in the alignment
                try:
                    index = chars.index(this357F)
                    offset = aligned_indexes[index] - 357
                except ValueError:
                    sys.stderr.write("Unable to determine 16S coordinates for " +
                          "region analysis\n")
                    sys.exit(1)
                # get where region of interest is in relation to alignment
                # constraining our coords to the sequence length
                # ie nothing longer than max or less than 0

                # rewrite the dict with the official coords with the adjusted coords
                try:
                    for k, v in regions.items():
                        # if k != "all":
                        #     continue
                        # print(k, v)
                        # print(offset)
                        # print(len(aligned_indexes))
                        # # use offset from 357 to find how much the cannonical sequence was trimmed
                        # # constrain to range (1, max number of index positions)
                        start_in_trim = min(max(v[0] + offset, 1), len(aligned_indexes))
                        end_in_trim = min(v[1] + offset, len(aligned_indexes))
                        # print(start_in_trim, end_in_trim)
                        # # adjust for 0 index
                        start_in_aln = aligned_indexes[start_in_trim - 1]
                        end_in_aln = aligned_indexes[end_in_trim - 1]
                        regions[k] = [start_in_aln, end_in_aln]
                except Exception as e:
                    print(k, v)
                    print("calculated indexes")
                    print("length of aligned_indexes",len(aligned_indexes))
                    print(start_in_trim, end_in_trim)
                    print(e)
                    sys.exit(1)
                # sys.stderr.write(
                #     "Calculating entropy for {}: E. coli {}: {}\n".format(
                #         region, regions[region][0], regions[region][1]))
                # sys.stderr.write(
                #     "This corresponds to alignment {}: {}\n".format(
                #         start_in_aln, end_in_aln))
                # print(regions)
                return regions


def shannon_calc(values):
    """ calculate shannon entropy for a position
    """
    possible = set(values.split())
    entropy = 0
    # for nuc in ["a", "t", "c", "g", "-"]:
    for nuc in ["a", "t", "c", "g"]:
        n = values.lower().count(nuc)
        if n == 0:
            continue
        prob = n/len(values)
        ent = prob * math.log(prob)
        entropy = entropy + ent
    return - entropy



def read_in_msa(path, ignore_first, pattern):
    """
    count occurances of each nucleotide in a sequence
    """
    seqs, seqs_subset  = [], []
    with open(path, "r") as inf:
        for i, rec in enumerate(SeqIO.parse(inf, "fasta")):
            if i == 0 and ignore_first:
                continue
            if pattern is not None:
                if pattern in rec.id or pattern in rec.description:
                    # exclude seqeunces matching pattern from subset
                    pass
                else:
                    seqs_subset.append(str(rec.seq).lower())
            seqs.append(str(rec.seq).lower())
    if pattern is not None:
        sys.stderr.write("seqs in full vs subset:  %i, %i\n" %(len(seqs), len(seqs_subset)))
    return (seqs, seqs_subset)


def main():
    args = get_args()
    seqs, seqs_subset = read_in_msa(
        args.input, ignore_first=args.colifirst,
        pattern=args.pattern)
    coord_regions = None
    #if args.colifirst:
    if True:
        coord_regions = detect_ecoli_positions(args.input)
        for i in range(1, 10):
            sys.stderr.write("V{}: {} -- {}\n".format(i,coord_regions[i][0], coord_regions[i][1]))

    thisregion = "-"
    # go position by position
    for i in range(len(seqs[0])):  #coords:

        values_string = "".join([x[i] for x in seqs])
        entropy = shannon_calc(values_string)
        if args.pattern:
            values_string_subset = "".join([x[i] for x in seqs_subset])
            entropy_subset = shannon_calc(values_string_subset)
        else:
            entropy_subset = "-"
        if coord_regions is not None:
            # to get around dicts being unsorted
            for j in range(1, 10):
                if coord_regions[j][0] <= i <= coord_regions[j][1]:
                    thisregion = "V"  + str(j)
            assert thisregion is not None, "error finding region for position %i" % i
        else:
            thisregion = "-"
        sys.stdout.write("{i}\t{thisregion}\t{entropy}\t{entropy_subset}\n".format(**locals()))
