#!/usr/bin/env python

import argparse
import subprocess
import sys
import shutil
import os
import subprocess
import math
from Bio import SeqIO
from collections import Counter


def get_args():
    parser = argparse.ArgumentParser()
    ##REQUIRED
    parser.add_argument("-d", "--silva", help="silva database", required=True)
    parser.add_argument("-o", "--silva_out", help="edited silva database", required=True)
    parser.add_argument("-n", "--org_name", help="organism name", required=True)
    parser.add_argument("-S", "--focus_seqs", help="path to focus sequences", required=False)
    return(parser.parse_args())


def count_orgs_silva(org, silva):
    ''' count the amount of sequences in the silva database for each species
    '''
    overall_count=0
    count = 0
    with open(silva, "r") as infile:
        for line in infile:
            if org in line:
                count += 1
    print("{}: {}".format(org, count))
            
def new_silvadb_for_org(count, org, silva, new_silva):
    ''' write new silva database with just sequences from org
    ''' 
    nlines = 0
    write_next_line = False
    print("extracting {org} sequences".format(**locals()))
    with open(silva, "r") as inf, open(new_silva, "a") as outf:
        for line in inf:
            
            if write_next_line:
                outf.write("{line}".format(**locals()))
                write_next_line = False
            elif org in line:
                outf.write("{line}".format(**locals()))
                write_next_line = True
                nlines = nlines + 1 
            else:
                pass
    totalcount = float(count + nlines)
    print("wrote %i lines entries" % nlines)
    print("Total sequences: {totalcount} ".format(**locals()))
      

def add_16db_seqs(focus_file, new_silva, org):
    ''' given a file containing sequences, cats them to silva file
    '''
    seqs = "./tmp"
    #puts each sequence on single line
    cmd = "seqtk seq -S {focus_file} > {seqs}".format(**locals())
    subprocess.run(cmd,
                   shell=sys.platform !="win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
    count=0
    with open(seqs, "r") as e:
        lines = e.readlines()
    with open(new_silva, "a") as f:
            for line in lines:
                line = rename_header_line(line, org)
                f.write(line)
                count += 1
    count = count / 2
    print("Adding {} sequences".format(count))
    os.remove(seqs)
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
    msa = os.path.join(multifasta + ".mafft")
    cmd = "mafft --retree 2 --reorder {multifasta} > {msa}".format(**locals())
    
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
    count_orgs_silva(org=org, silva=silva)

    #Adds sequences from a 16db ribo16s file to new silva database
    if args.focus_seqs is not None:
        focus_seqs = args.focus_seqs
        count = add_16db_seqs(focus_file=focus_seqs, new_silva=new_silva, org=org)
    else:
        count=0
        
    #Adds only sequences of that organism from silva database
    new_silvadb_for_org(org=org, silva=silva, new_silva=new_silva, count=count)
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
