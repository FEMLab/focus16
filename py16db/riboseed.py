#!/usr/bin/env python 

import subprocess
import os
import argparse

def get_args():
    parser = argparse.ArgumentParser():
    parser.add_argument("-r", "--reference", 
                        help="path to best_ref", required=True)
    parser.add_argument("-R", "--reads",
                        help="path to reads file", required=True)
    parser.add_argument("--cores", 
                        help="number of cores to use", required=True)
    parser.add_argument("--threads",
                        help="threads to use", required=True)
    parser.add_argument("-v",
                        help="???", required=True)
    parser.add_argument("--serialize",
                        help="???", required=True)
    parser.add_argument("--output",
                       help="desired output", required=True)
    return(parser.parse_args())


def riboseed(sra, Freads, Rreads, cores, threads, v,  serialize, output):
    best_ref = os.path(args.reference)
    for line in best_ref
         sra = []
         line.split('\t')
         sra.append([0])
    cmd = "ribo run -r" + sra + "-F" + Freads + "-R" + Rreads + "--cores" + cores +
    "--threads" + threads + "-v" + v + "--serialize" + serialize + "-o" + suboutput_dir
    subprocess.run(cmd, 
                   shell=sys.platform !="win32",
                   stdout=subprocess.PIPE,
                   sterr=subprocess.PIPE,
                   check=True)
    return()

if __name__ == '__main__':
    args=get_args()
    best_ref = os.path(args.reference)
    for line in best_ref
         sra = []
         line.split('\t')
         sra.append([0])
    Freads = os.path.join(args.reads, "*1.fq")
    Rreads = os.path.join(args.reads, "*2.fq")
    riboseed(sra=sra, Freads, Rreads, cores=args.cores, threads=args.threads, v=args.v,  serialize=args.serialize, output=args.output)
    
    
