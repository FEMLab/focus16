#!/usr/bin/env python                                                                                                                                                                

import sys
import os
import collections
import pickle

def singleSRA():
    org_dict = {}
    sraFind = "/Users/alexandranolan/Desktop/sraFind/results/sraFind-All-biosample-with-SRA-hits.txt"
    with open(sraFind, "r") as infile:
        for i,line in enumerate(infile):
            split_line = [x.replace('"', '').replace("'", "") for x in line.strip().split("\t")]
            gs = " ".join(split_line[11].split(" ")[0:2])
            if gs in org_dict.keys():
                org_dict[gs] += 1
            else:
                org_dict[gs] = 1
            if split_line[8].startswith("ILLUMINA"):
                singleSRA = []
                singleSRA.append(split_line[11]) #organism name
                print(org_dict)
        return()                    
                      
         
def main():
    singleSRA()


if __name__ == '__main__':
    main()
                 
