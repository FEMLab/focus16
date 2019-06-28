#!/usr/bin/env python

import argparse
import sys
import os
import subprocess
import random
import shutil
import re
import gzip

from . import get_n_genomes as gng
from . import fetch_sraFind_data as fsd
from . import run_sickle
from Bio import SeqIO

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output_dir", help="path to output", required=True)
    parser.add_argument("-n", "--organism_name", help="genus species in quotes", required=True)
    parser.add_argument("-l", "--approx_length", help="approximate genome length", required=True, type=int)
    parser.add_argument("-s", "--sraFind_path", dest="sra_path", default="srapure",
                        help="path to sraFind file", required=False)
    parser.add_argument("--single_SRA", default=None,
                        help="run pipeline on this SRA accession only", required=False)
    parser.add_argument( "-g", "--genomes_dir",help="path to genomes directory",
                         required=True)
    parser.add_argument("-p", "--prokaryotes", action="store",
                        help="path_to_prokaryotes.txt", default="./prokaryotes.txt",
                        required=False)
    parser.add_argument("-S", "--nstrains", help="number of strains",
                        type=int, required=True)
    parser.add_argument("--get_all", help="get both sras if organism has two",
                        action="store_true", required=False)
    parser.add_argument("--cores", help="number of cores for riboSeed", default=1,
                        required=False, type=int)
    parser.add_argument("--maxcov", help="maximum coverage of reads", default=50,
                        required=False, type=int)
    parser.add_argument("--example_reads", help="input of example reads", nargs='+',
                        required=False, type=str)
    return(parser.parse_args())


def check_programs():
    required_programs = ["ribo", "barrnap", "fastq-dump", "mash", "skesa", "plentyofbugs"]
    for program in required_programs:
        if shutil.which(program) is None:
            print("{program} not installed".format(**locals()))
            sys.exit(1)
   
           
#srapure = file containing only [0]accession and [1]genuspecies for 16S sequences

def filter_srapure(path, organism_name, strains, get_all):
    """For srapure or sraFind, take specified number of organisms and SRAs"""
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

def download_SRA(SRA, destination):
    """Download sras from filter_srapure into raw reads file, downsample forward reads to 1000000"""
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
    downcmd = "seqtk sample -s100 " + suboutput_dir_raw + SRA + "_1.fastq 1000000 > " + downpath
    
    subprocess.run(downcmd,
                  shell=sys.platform !="win32",
                  stdout=subprocess.PIPE,
                  stderr=subprocess.PIPE,
                   check=True)
    return(downpath)


def pob(genomes_dir, readsf, output_dir):
    """Uses plentyofbugs, a package that uses mash to find the best reference genome for draft genome """

    ds_reads = os.path.join(os.path.dirname(output_dir), "reads_for_pob.fastq")
    downcmd2 = "seqtk sample -s100 {readsf} 1000000 > {ds_reads}".format(**locals())
        
    pobcmd = "plentyofbugs -g " + genomes_dir +  " -f " + ds_reads + " -o " + output_dir
    for command in [downcmd2, pobcmd]:
        print(command)
        subprocess.run(command,
                   shell=sys.platform !="win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
    best_ref = os.path.join(output_dir, "best_reference")
    with open(best_ref, "r") as infile:
        for line in infile:
            print(line)
            sraacc = line.strip().split('\t')
            return(sraacc)


def get_ave_read_len_from_fastq(fastq1, N=50):
    """from LP: taken from github.com/nickp60/riboSeed/riboSeed/classes.py; return average read length in fastq1 file from first N reads"""
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
    return(ave_read_len)

def get_coverage(approx_length, fastq1):
    """Obtains the coverage for a read set, when given the estimated genome size"""
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
    """Given the coverage from coverage(), downsamples the reads if over the max coverage set"""
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
    

def run_riboseed(sra, readsf, readsr, cores, threads, output):
    """Runs riboSeed to reassemble reads """
    cmd = "ribo run -r {sra} -F {readsf} -R {readsr} --cores {cores} --threads {threads} -v 1 --serialize -o {output} --subassembler skesa --additional_libs '--only-assembler'".format(**locals())
    if readsr is None:
        cmd = "ribo run -r {sra} -F {readsf} --cores {cores} --threads {threads} -v 1 --serialize -o {output} --subassembler skesa".format(**locals())
    print(cmd)
    subprocess.run(cmd,
                   shell=sys.platform !="win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
    return()

def  extract_16s_from_contigs(input_contigs, barr_out, output):
    """Uses barrnap to identify rRNA operons within the riboSeed assembled contigs, then uses extractRegion to extract the 16S sequences """
    barrnap = "barrnap {input_contigs} > {barr_out}".format(**locals())
    subprocess.run(barrnap,
                   shell=sys.platform !="win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
    
    results16s = []   # [chromosome, start, end, reverse complimented]
    with open(barr_out, "r") as rrn:
        for rawline in rrn:
            line = rawline.strip().split('\t')
            if line[0].startswith("##"):
               continue
            if line[8].startswith("Name=16S"):
                if line[6] == "-":
                    suffix = 'strep-RC@'
                else:
                    suffix = ''
                chrom=line[0]
                start=line[3]
                end=line[4]
                results16s = [chrom, start, end, suffix]
                print(results16s)    
                cmd = "extractRegion \'{results16s[3]}{results16s[0]} :{results16s[1]}:{results16s[2]}\' -f {input_contigs} -v 1 >> {output}".format(**locals())
                subprocess.run(cmd,
                               shell=sys.platform !="win32",
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               check=True)
                
    return(output)

def alignment(fasta, output):
    output = os.path.join(output, "alignment", "")
    os.makedirs(output)
    seqout = os.path.join(output, "16soneline.fasta")
    mafftout = os.path.join(output, "MSA.fasta")
    #iqtreeout = os.path.join(output, "iqtree")
    seqcmd = "seqtk seq -S {fasta} > {seqout}".format(**locals())
    mafftcmd = "mafft {seqout} > {mafftout}".format(**locals())
    iqtreecmd = "iqtree -s {mafftout} -nt AUTO".format(**locals())
    for cmd in [seqcmd, mafftcmd, iqtreecmd]:
        subprocess.run(cmd,
                       shell=sys.platform !="win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)
    return(mafftout)
            

def process_strain(rawreadsf, rawreadsr, this_output, args):
    pob_dir=os.path.join(this_output, "plentyofbugs")
    ribo_dir=os.path.join(this_output, "riboSeed")
    sickle_out=os.path.join(this_output, "sickle")
    get_ave_read_len_from_fastq(fastq1=rawreadsf, N=50)
    pob(genomes_dir=args.genomes_dir, readsf=rawreadsf, output_dir=pob_dir)
    best_reference=os.path.join(this_output, "plentyofbugs/best_reference")
    with open(best_reference, "r") as infile:
        for line in infile:
            best_ref_fasta = line.split('\t')[0]
    trimmed_fastq1, trimmed_fastq2 = run_sickle(fastq1=rawreadsf, fastq2=rawreadsr, output_dir=sickle_out)        
    print("Trimmed f reads: ", trimmed_fastq1)
    print("Trimmed r reads: ", trimmed_fastq2)
                       
    downsampledf, downsampledr = downsample(
        approx_length=args.approx_length, fastq1=trimmed_fastq1, fastq2=trimmed_fastq2, maxcoverage=args.maxcov, destination=this_output)
    print("Downsampled f reads: ", downsampledf)
    print("Downsample r reads: ", downsampledr)
        
    run_riboseed(sra=best_ref_fasta, readsf=downsampledf,
                 readsr=downsampledr, cores=args.cores,
                 threads=1, output=ribo_dir)
    
    barr_out=os.path.join(ribo_dir, "barrnap")
    sixteens_extracted=os.path.join(ribo_dir, "ribo16s")
    extract_16s_from_contigs(input_contigs=ribo_contigs, barr_out=barr_out, output=sixteens_extracted)
    alignment(fasta=sixteens_extracted, output=this_output)
    

if __name__ == '__main__':
    # main()
    args=get_args()
    check_programs()
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
    if args.example_reads is not None:
        this_output=os.path.join(args.output_dir, "example", "")
        os.makedirs(this_output)
        rawreadsf = args.example_reads[0]
        try:
            rawreadsr = args.example_reads[1]
        except IndexError:
            rawreadsr = None
        
        process_strain(rawreadsf, rawreadsr, this_output, args)
    else:
        for i, accession in enumerate(filtered_sras):
            this_output=os.path.join(args.output_dir, str(i))
            os.makedirs(this_output)
            print("Downloading " + accession)
            download_SRA(SRA=accession, destination=this_output)
        
            rawreadsf=os.path.join(this_output, "raw/", accession + "_1.fastq")
            rawreadsr=os.path.join(this_output, "raw/", accession + "_2.fastq")
            process_strain(rawreadsf, rawreadsr, this_output, args)


