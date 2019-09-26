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
        description="Given a path to Silva 16s database and the results " +
        "from 16db, build a species-specific database" )
    # REQUIRED
    parser.add_argument("-d", "--silva", help="silva database", required=True)
    #parser.add_argument("-o", "--out", help="Combined database", required=True)
    parser.add_argument("-n", "--org_name", help="organism name", required=True)
    parser.add_argument("-S", "--focus_seqs", help="path to focus sequences", required=False)
    # Optional
    parser.add_argument("--lower", help="lowecase sequence", action="store_true")
    parser.add_argument("--rna", help="transcribe DNA to RNA", action="store_true")
    parser.add_argument("--dna", help="RNA to DNA", action="store_true")
    return(parser.parse_args())


# def count_orgs_silva(org, silva):
#     ''' count the amount of sequences in the silva database for each species
#     '''
#     overall_count=0
#     count = 0
#     if os.path.splitext(silva)[-1] in ['.gz', '.gzip']:
#         open_fun = gzip.open
#     else:
#         open_fun = open
#     with open_fun(silva, "rt") as infile:
#         for line in infile:
#             if org in line:
#                 count += 1
#     print("{}: {}".format(org, count))

def new_silvadb_for_org(count, org, silva, lower, dna, rna):
    ''' write new silva database with just sequences from org
    '''
    nlines = 0
    write_next_line = False
    sys.stderr.write("extracting {org} sequences\n".format(**locals()))
    if os.path.splitext(silva)[-1] in ['.gz', '.gzip']:
        open_fun = gzip.open
    else:
        open_fun = open
    with open_fun(silva, "rt") as inf:
        for rec, seq in SimpleFastaParser(inf):
            if org in rec:
                if rna:
                    seq = Seq(seq).transcribe()
                if dna:
                    seq = Seq(seq).back_transcribe()
                if lower:
                    seq = seq.lower()
                sys.stdout.write(">%s\n%s\n" % (rec, seq))
                nlines += 1
    totalcount = count + nlines
    sys.stderr.write("SILVA seqeunces: {}\n".format(nlines))
    sys.stderr.write("Combined sequences: {totalcount}\n".format(**locals()))


def add_16db_seqs(focus_file, org, dna, rna, lower):
    ''' adds seqeunces to silva file as single-line fasta
    '''
    count=0
    with open(focus_file, "r") as inf:
        for rec, seq in SimpleFastaParser(inf):
            # for rec in SeqIO.parse(inf,"fasta"):
            if rna:
                seq = Seq(seq).transcribe()
            if dna:
                seq = Seq(seq).back_transcribe()
            if lower:
                seq = seq.lower()
            sys.stdout.write(">%s\n%s\n" % (rec, seq))
            count = count + 1
    sys.stderr.write("focusDB seqeunces: {}\n".format(count))
    # os.remove(seqs)
    return(count)


def rename_header_line(line, org):
    ''' for header line in file, rename to org name
    '''
    if line.startswith(">"):
        headerline = line.replace(":","_").replace(" ","_")
        return headerline
    else:
        return(line)
    sys.stderr.write("Renaming sequences\n")


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

    #print(Counter(leaddash))

    return seqs


def main():
    args = get_args()
    silva = args.silva
    #new_silva = args.out
    org=args.org_name

    #count occurances of organism in silva database
    # count_orgs_silva(org=org, silva=silva)

    # Adds sequences from a 16db ribo16s file to new silva database
    if args.focus_seqs is not None:
        sys.stderr.write("Adding results from focusDB\n")
        focus_seqs = args.focus_seqs
        count = add_16db_seqs(
            focus_file=focus_seqs,
            #new_silva=new_silva,
            rna=args.rna,
            dna=args.dna,
            lower=args.lower,
            org=org)
    else:
        count=0

    #Adds only sequences of that organism from silva database
    new_silvadb_for_org(
        org=org,
        silva=silva,
        #new_silva=new_silva,
        count=count,
        rna=args.rna,
        dna=args.dna,
        lower=args.lower
    )
    # msa = mafft(multifasta = new_silva)

    # entropies = []
    # seqs = read_in_msa(msa)
    # shannon_out = os.path.join(new_silva + ".shannon")

    # for i in range(len(seqs[0])):
    #     values_string = "".join([x[i] for x in seqs])
    #     entropy = shannon_calc(values_string)
    #     entropies.append(entropy)
    # for i, value in enumerate(entropies):
    #     with open(shannon_out, "a") as f:
    #         f.write("{i}\t{value}\n".format(**locals()))
    # print("Shannon entropy complete")


if __name__ == '__main__':
    main()
