#!/usr/bin/env python

import argparse
import sys
import os
import subprocess
import random
import shutil

from Bio import SeqIO

import get_n_genomes as gng
import fetch_sraFind_data as fsd

def get_args():
    parser = argparse.ArgumentParser(
    description = "given a genus species, downloads the raw reads for each SRA using " +
    "fastq-dump. Also downloads downsampled [1000000 reads] version of each SRA")
    parser.add_argument("-o", "--output_dir", help="path to output", required=True)
    parser.add_argument("-n", "--organism_name", help="genus species in quotes", required=True)
    parser.add_argument("-l", "--approx_length", help="approximate genome length", required=True)
    parser.add_argument("-s", "--sraFind_path", dest="sra_path",
                        default="srapure",
                        help="path to sraFind file", required=False)
    parser.add_argument("--single_SRA",
                        default=None,
                        help="run pipeline on this SRA accession only",
                        required=False)
    parser.add_argument(
        "-g",
        "--genomes_dir",
        help="path to genomes directory",
        required=True)
    parser.add_argument(
        "-p",
        "--prokaryotes",
        action="store",
        help="path_to_prokaryotes.txt",
        default="./prokaryotes.txt",
        required=False)
    parser.add_argument("-S", "--nstrains", help="number of strains",
                        type=int, required=True)
    parser.add_argument("--get_all", help="get both sras if organism has two",
                        action="store_true",
                        required=False)
    parser.add_argument("--cores", help="number of cores for riboSeed",
                        required=False)
    return(parser.parse_args())


#srapure = file containing only [0]accession and [1]genuspecies for 16S sequences
#remove " and ' from the lines

def filter_srapure(path, organism_name, strains, get_all):
    results = []
    with open(path, "r") as infile:
        for line in infile:
            split_line = [x.replace('"', '').replace("'", "") for x in line.strip().split("\t")]
            if split_line[11].startswith(organism_name):
                results.append(split_line[17])
    random.shuffle(results)
    print(results)
    if strains != 0:
        results = results[0:strains]

    sras = []
    for result in results:
        these_sras = result.split(",")
        if get_all:
            for sra in these_sras:
                sras.append(sra)
        else:
            sras.append(these_sras[0])
    return(sras)

#create two suboutputs, raw and downampled
#create a variable that is the fastq-dump command, with relevant SRA and output

def download_SRA(SRA, destination):
    suboutput_dir_raw = os.path.join(destination, "raw", "")
    suboutput_dir_downsampled = os.path.join(destination, "downsampled", "")
    os.makedirs(suboutput_dir_raw)
    os.makedirs(suboutput_dir_downsampled)
    cmd = "fastq-dump --split-files " + SRA + " -O " + suboutput_dir_raw
    subprocess.run(cmd,
                   shell=sys.platform !="win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
    downpath = os.path.join(suboutput_dir_downsampled, "downsampledreadsf.fastq")
    downcmd = "seqtk sample -s100 " + suboutput_dir_raw + SRA + "_1.fastq 0.5 > " + downpath
    
    subprocess.run(downcmd,
                  shell=sys.platform !="win32",
                  stdout=subprocess.PIPE,
                  stderr=subprocess.PIPE,
                   check=True)
    return(downpath)


def pob(genomes_dir, readsf, output_dir):
    pobcmd = "plentyofbugs -g " + genomes_dir +  " -f " + readsf + " -o " + output_dir
    subprocess.run(pobcmd,
                   shell=sys.platform !="win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
    best_ref = os.path.join(output_dir, "best_reference")
    sra = []
    with open(best_ref, "r") as infile:
        for line in infile:
            print(line)
            sraacc = line.strip().split('\t')
            return(sraacc)

#taken from github.com/nickp60/riboSeed/riboSeed/classes.py
def get_ave_read_len_from_fastq(fastq1, N=50):
    """from LP; return average read length in fastq1 file from first N reads
    """
    count, tot = 0, 0
    if os.path.splitext(fastq1)[-1] in ['.gz', '.gzip']:
        open_fun = gzip.open
    else:
        open_fun = open
    with open_fun(fastq1, "rt") as file_handle:
        data = SeqIO.parse(file_handle, "fastq")
        for read in data:
            count += 1
            tot += len(read)
    ave_read_len = float(tot / count)  
    return(ave_read_len, count)

def get_coverage(approx_length, fastq1):
    count, tot = 0, 0
   
    if os.path.splitext(fastq1)[-1] in ['.gz', '.gzip']:
        open_fun = gzip.open
    else:
        open_fun = open
    with open_fun(fastq1, "rt") as file_handle:
        data = SeqIO.parse(file_handle, "fastq")
        for read in data:
            count += 1
            tot += len(read)

    ave_read_len = float(tot / count)  
    coverage = float((count * ave_read_len) / approx_length)
    print(coverage)
    return(coverage)

def downsample(approx_length, fastq1, fastq2, maxcoverage, destination):
    suboutput_dir_raw = os.path.join(destination, "raw", "")
    suboutput_dir_downsampled = os.path.join(destination, "downsampled", "")
    downpath1 = os.path.join(suboutput_dir_downsampled, "downsampledreadsf.fastq") 
    downpath2 = os.path.join(suboutput_dir_downsampled, "downsampledreadsr.fastq")
    coverage = get_coverage(approx_length, fastq1)
    
    covfraction = round(float(maxcoverage / coverage), 3)
    print(covfraction)
    if (coverage > maxcoverage):
         downcmd = "seqtk sample -s100 {fastq1} {covfraction} > {downpath1}".format(**locals())
         downcmd2 = "seqtk sample -s100 {fastq2} {covfraction} > {downpath2}".format(**locals())
         for command in [downcmd, downcmd2]:
             print(command)
             subprocess.run(command,
                            shell=sys.platform !="win32",
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            check=True)
         return(downpath1, downpath2)
    else:
        return(fastq1, fastq2)
    

def riboseed(sra, readsf, readsr, cores, threads, v, output):
    cmd = "ribo run -r " + sra + " -F " + readsf + " -R " + readsr + " --cores " + cores + " --threads " + threads + " -v " + v + " --serialize -o " + output
    subprocess.run(cmd,
                   shell=sys.platform !="win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
    return()


if __name__ == '__main__':
    # main()
    args=get_args()
    os.makedirs(args.output_dir)
    if not os.path.exists(args.sra_path):
        sraFind_output_dir = os.path.join(args.output_dir, "sraFind")
        args.sra_path = fsd.main(args, output_dir=sraFind_output_dir)
    if args.single_SRA is None:
        filtered_sras = filter_srapure(
            path=args.sra_path,
            organism_name=args.organism_name,
            strains=args.nstrains,
            get_all=args.get_all)
    else:
        filtered_sras = [args.single_SRA]
    if os.path.exists(args.genomes_dir):
        if len(os.listdir(args.genomes_dir)) == 0:
            print("Warning: genome directory exists but is " +
                  "empty: downloading genomes")
            shutil.rmtree(args.genomes_dir)
            gng.main(args)
        else:
            pass
    else:
        gng.main(args)
    if filtered_sras == []:
        print("Organism not found on sraFind")
    for i, accession in enumerate(filtered_sras):
        this_output=os.path.join(args.output_dir, str(i))
        os.makedirs(this_output)
        print("Downloading " + accession)
        download_SRA(SRA=accession, destination=this_output)
        
        rawreadsf=os.path.join(this_output, "raw/*_1.fastq")
        rawreadsr=os.path.join(this_output, "raw/*_2.fastq")
        downreadsf=os.path.join(this_output, "downsampled/downsampledreadsf.fastq")
        downreadsr=os.path.join(this_output, "downsampled/downsampledreadsr.fastq")
        pob_dir=os.path.join(this_output, "plentyofbugs")
        ribo_dir=os.path.join(this_output, "riboSeed")
        
        get_ave_read_len_from_fastq(fastq1=downreadsf, N=50)
        pob(genomes_dir=args.genomes_dir, readsf=downreadsf, output_dir=pob_dir)
        best_reference=os.path.join(this_output, "plentyofbugs/best_reference")
        sra=[]
        for sra in best_reference: 
            best_ref = sra.strip().split('/t')
            sra = sra.append(best_ref[0])
            for acc in sra:
                if acc == os.path.basename(genomes_dir): 
                    print(os.path(genomes_dir)) 
                    
                    coverage(approx_length=args.approx_length,
                             fastq1=rawreadsf)
                    
                    downsample(approx_length=args.approx_length, fastq1=rawreadsf,
                               maxcoverage=50, destination=this_output)
                    
                    if os.path.exists(downreadsr):    
                        riboseed(sra=best_reference, readsf=downreadsf,
                                 readsr=downreadsr, cores=args.cores,
                                 threads=1, v=1, output=ribo_dir)
                    else:
                        riboseed(sra=best_reference, readsf=rawreadsf,
                                 readsr=rawreadsr, cores=args.cores,
                                 threads=1, v=1, output=ribo_dir)
                    
