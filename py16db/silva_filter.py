#!/usr/bin/env python

import argparse
import subprocess
import sys
import shutil
import os

def get_args():
    parser = argparse.ArgumentParser()
    ##REQUIRED
    parser.add_argument("-d", "--silva", help="silva database", required=True)
    parser.add_argument("-o", "--silva_out", help="edited silva database", required=True)
    parser.add_argument("-n", "--org_name", help="organism name", required=True)
    parser.add_argument("-S", "--db16_seqs", help="path to 16db sequences", required=False)
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
      

def add_16db_seqs(db16_file, new_silva, org):
    ''' given a file containing sequences, cats them to silva file
    '''
    seqs = "./tmp"
    #puts each sequence on single line
    cmd = "seqtk seq -S {db16_file} > {seqs}".format(**locals())
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
        
def rename_headers(db16_file, org):
    ''' for header line in file, rename to org name
    '''
    with open(db16_file, "r") as e:
        lines = e.readlines()
        for line in lines:
            if line.startswith(">"):
                headerline = line.replace(">" , ">{org}_".format(**locals()))
                f.write(headerline)
            else:
                f.write(line)
    print("Renaming sequences")

def rename_header_line(line, org):
    ''' for header line in file, rename to org name
    '''
    if line.startswith(">"):
        headerline = line.replace(">" , ">{org}_".format(**locals())).replace(":","_").replace(" ","").replace("_","")
        
        return headerline
    else:
        return(line)
    
            

def main():
    args = get_args()
    silva = args.silva                                   
    new_silva = args.silva_out 
    
    org=args.org_name

    #count occurances of organism in silva database
    count_orgs_silva(org=org, silva=silva)

    #Adds sequences from a 16db ribo16s file to new silva database
    if args.db16_seqs is not None:
        db16_seqs = args.db16_seqs
        count = add_16db_seqs(db16_file=db16_seqs, new_silva=new_silva, org=org)
    else:
        count=0
        
    #Adds only sequences of that organism from silva database
    new_silvadb_for_org(org=org, silva=silva, new_silva=new_silva, count=count)

    
    
    
if __name__ == '__main__':
    main()
