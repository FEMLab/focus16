#!/usr/bin/env python

import argparse
import math
import sys
from collections import Counter
from Bio import SeqIO
#shannon entropy calculator for sequences
#in a multiple-sequence alignment

#H(i)= -[(pi*log2*pi)]
#H>2: variable
#H<2: conserved
#Whichever H value is highest is the most conserved 
#nucleotide for that position

#5 possible outcomes for each base: [A,T,C,G,-]

def argParse():
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--mafft_alignment", required=True)
    return(parser.parse_args())

#calculate shannon entropy for a position
def shannon_calc(values):
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

#count occurances of each nucleotide in a sequence
def read_in_msa(path):
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

    print(Counter(leaddash))
    sys.exit(1)
    return seqs


def main():
    args=argParse()
    #shannon_for_msa(alignment=args.mafft_alignment, output=args.shannon_output)
    entropies = []
    seqs = read_in_msa(args.mafft_alignment)

    for i in range(len(seqs[0])):
        values_string = "".join([x[i] for x in seqs])
        entropy = shannon_calc(values_string)
        entropies.append(entropy)
    for i, value in enumerate(entropies):
        sys.stdout.write("{i}\t{value}\n".format(**locals()))

if __name__ == '__main__':
    main()
