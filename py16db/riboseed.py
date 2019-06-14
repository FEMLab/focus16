#!/usr/bin/env python 

import subprocess
import os
import argparse

def get_args():
    parser = argparse.ArgumentParser():
    description = "Run riboseed with reads obtained from " + 
    "get_sra_for_organism.py and reference genome obtained from " +
    "best_ref.py"
    parser.add_argument("-p", "--path", 
                        help="path to numbered directory, given for each SRA",
                        required=True)
    parser.add_argument("--threads",
                        help="threads to use", required=True)
    parser.add_argument("-v",
                        help="???", required=True)
    parser.add_argument("--serialize",
                        help="???", required=True)
    return(parser.parse_args())

#mash [best_ref.py] produces a file with "acc/tpercentage"
def riboseed(sra, Freads, Rreads, cores, threads, v,  serialize, output):
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
    for i in range(1, 11):
        best_ref = os.path.join(args.path, i, "pathtoref")
        with open(best_ref, "r") as infile:
            for line in best_ref:
                sra = []
                sraacc = line.strip().split('\t')
                sra.append(sraacc[0])
#this inputs the ref accession, the fasta for the ref is needed
        Freads = os.path.join(args.path, i, "downsampled/reads1.fq")
        Rreads = os.path.join(args.path, i, "downsampled/reads2.fq") 
        output = os.path.join(args.path, i)
        riboseed(sra=sra, Freads=Freads, Rreads=Rreads,
                 cores=args.cores, threads=args.threads, v=args.v,
                 serialize=args.serialize, output=output)

    
#'for ref in {1..6}; do 
#bestref=$(cat ~/FYP/data/ovatus/${ref}/plentyofbugs/best_reference | cut -f 1);
#get_genomes.py -q $bestref -o ~/FYP/results/ovatus/${ref}/;
#ribo run -r   ~/FYP/results/ovatus/${ref}/${bestref}.fasta -F ~/FYP/data/ovatus/${ref}/downsampled/reads1.fq -R ~/FYP/data/ovatus/${ref}/downsampled/reads2.fq  --cores 16 --threads 1 -v 1 --serialize -o ~/FYP/results/ovatus/${ref}/riboseed ; done'
#- Incorporate this into this script
