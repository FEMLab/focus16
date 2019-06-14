#!/usr/bin/env python 

import argparse
import os
import sys
import subprocess

#argument function
def get_args():
    parser = argparse.ArgumentParser()
    description = "Finds the best reference genome for a SRA with Mash"
    parser.add_argument(
        "-f",
        "--forward_reads",
        help="forward reads file",
        required=True)
    parser.add_argument(
        "-o",
        "--output",
        help="output directory",
        required=False)
    parser.add_argument(
        "-g",
        "--genomes_dir",
        help="path to output"
        required=True)
    #parser.add_argument(
        "-r",
        "--reverse_reads",
        help="reverse reads file",
        required=True)
   #parser.add_argument(
        "-n",
        "--nstrains",
        help="number of genomes that will be compared as the reference",
        required=True)
    #parser.add_argument(
        "--genus_species",
        help="organism name",
        required=True)
    #parser.add_argument(
        "--assembler",
        help="skesa or spades",
        required=True)
    #parser.add_argument("-m",
                        "--pyani_mash",
                        help="pyani or mash for genome selection",
                        required=True)
    
                                     
    return(parser.parse_args())

#function to form command to use mash
def pob(genomes_dir, readsf)  
    suboutput_dir = os.path.join(output, "mash")
    os.makedirs(suboutput_dir)
    pobcmd = "plentyofbugs -g" + genomes_dir +  "-f" + readsf 
    subprocess.run(mashcmd,
                   shell=sys.platform !="win32",
                   stdout=subprocess.PIPE,
                   sterr=subprocess.PIPE,
                   check=True)
    return()

#calling the above function 
if __name__ == '__main__':
 args=get_args()
 for genome in inputreads: 
     pob(genomes_dir=args.genomes_dir, readsf=forward_reads)
     
 
     
 
  
 
 
    
