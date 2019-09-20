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
    parser.add_argument("-o", "--silva_out", help="edited silva database", required=True)
    parser.add_argument("-n", "--org_name", help="organism name", required=True)
    parser.add_argument("-S", "--focus_seqs", help="path to focus sequences", required=False)
    # Optional
    parser.add_argument("--lower", help="lowecase sequence", action="store_true")
    parser.add_argument("--rna", help="transcribe DNA to RNA", action="store_true")
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

def new_silvadb_for_org(count, org, silva, new_silva, lower, transcribe):
    ''' write new silva database with just sequences from org
    '''
    nlines = 0
    write_next_line = False
    print("extracting {org} sequences".format(**locals()))
    if os.path.splitext(silva)[-1] in ['.gz', '.gzip']:
        open_fun = gzip.open
    else:
        open_fun = open
    with open_fun(silva, "rt") as inf, open(new_silva, "a") as outf:
        for rec, seq in SimpleFastaParser(inf):
            if org in rec:
                if transcribe:
                    seq = Seq(seq).transcribe()
                if lower:
                    seq = seq.lower()
                outf.write(">%s\n%s\n" % (rec, seq))
                nlines += 1
    totalcount = float(count + nlines)
    print("wrote %i lines entries" % nlines)
    print("Total sequences: {totalcount} ".format(**locals()))


def add_16db_seqs(focus_file, new_silva, org, transcribe, lower):
    ''' adds seqeunces to silva file as single-line fasta
    '''
    count=0
    with open(focus_file, "r") as inf, open(new_silva, "w") as f:
        for rec, seq in SimpleFastaParser(inf):
            # for rec in SeqIO.parse(inf,"fasta"):
            if transcribe:
                seq = Seq(seq).transcribe()
            if lower:
                seq = seq.lower()
            f.write(">%s\n%s\n" % (rec, seq))
            count = count + 1
    print("Adding {} sequences".format(count))
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
    print("Renaming sequences")

def mafft(multifasta):
    ''' performs default mafft alignment
    '''
    print("Aligning sequences with Mafft:")
    msa = os.path.join(multifasta + ".mafft")
    cmd = "mafft --retree 2 --reorder {multifasta} > {msa}".format(**locals())
    print(cmd)
    subprocess.run(cmd,
                   shell=sys.platform !="win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
    print("MSA complete")
    return(msa)


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

    #print(Counter(leaddash))

    return seqs


def main():
    args = get_args()
    silva = args.silva
    new_silva = args.silva_out
    org=args.org_name

    #count occurances of organism in silva database
    # count_orgs_silva(org=org, silva=silva)

    # Adds sequences from a 16db ribo16s file to new silva database
    if args.focus_seqs is not None:
        print("Adding results from focusdb")
        focus_seqs = args.focus_seqs
        count = add_16db_seqs(
            focus_file=focus_seqs,
            new_silva=new_silva,
            transcribe=args.rna,
            lower=args.lower,
            org=org)
    else:
        count=0

    #Adds only sequences of that organism from silva database
    new_silvadb_for_org(
        org=org,
        silva=silva,
        new_silva=new_silva,
        count=count,
        transcribe=args.rna,
        lower=args.lower
    )
    msa = mafft(multifasta = new_silva)

    entropies = []
    seqs = read_in_msa(msa)
    shannon_out = os.path.join(new_silva + ".shannon")

    for i in range(len(seqs[0])):
        values_string = "".join([x[i] for x in seqs])
        entropy = shannon_calc(values_string)
        entropies.append(entropy)
    for i, value in enumerate(entropies):
        with open(shannon_out, "a") as f:
            f.write("{i}\t{value}\n".format(**locals()))
    print("Shannon entropy complete")


if __name__ == '__main__':
    main()
