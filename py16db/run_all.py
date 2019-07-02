#!/usr/bin/env python

import argparse
import sys
import os
import subprocess
import random
import shutil
import re
import gzip
import logging

from . import get_n_genomes as gng
from . import fetch_sraFind_data as fsd
from py16db.run_sickle import run_sickle
from Bio import SeqIO

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output_dir", help="path to output", required=True)
    parser.add_argument("-n", "--organism_name", help="genus species in quotes",
                        required=True)
    parser.add_argument("-l", "--approx_length", help="approximate genome length",
                        required=True, type=int)
    parser.add_argument("-s", "--sraFind_path", dest="sra_path", default="srapure",
                        help="path to sraFind file", required=False)
    parser.add_argument("--single_SRA", default=None,
                        help="run pipeline on this SRA accession only",
                        required=False)
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
    """exits if the following programs are not installed"""
    
    required_programs = ["ribo", "barrnap",
                         "fastq-dump", "mash", "skesa", "plentyofbugs"]
    for program in required_programs:
        if shutil.which(program) is None:
            logging.error ( '%s is not installed: exiting.', program) 
            sys.exit(1)
   
           
def filter_srapure(path, organism_name, strains, get_all):
    """sraFind [github.com/nickp60/srafind], contains"""
    results = []
    with open(path, "r") as infile:
        for line in infile:
            split_line = [x.replace('"', '').replace("'", "") for x in line.strip().split("\t")]
            if split_line[11].startswith(organism_name):
                results.append(split_line[17])
    random.shuffle(results)
    logging.debug('Found SRAs: %s', results)
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
        logging.debug('Finding best reference genome: %s', command)
        subprocess.run(command,
                   shell=sys.platform !="win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
    best_ref = os.path.join(output_dir, "best_reference")
    with open(best_ref, "r") as infile:
        for line in infile:
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
    logging.debug('Read coverage is : %s', coverage)
    return(coverage)

def downsample(approx_length, fastq1, fastq2, maxcoverage, destination):
    """Given the coverage from coverage(), downsamples the reads if over the max coverage set by args.maxcov. Default 50."""
    suboutput_dir_raw = os.path.join(destination, "raw", "")
    suboutput_dir_downsampled = os.path.join(destination, "downsampled", "")
    downpath1 = os.path.join(suboutput_dir_downsampled, "downsampledreadsf.fastq") 
    downpath2 = os.path.join(suboutput_dir_downsampled, "downsampledreadsr.fastq")
    coverage = get_coverage(approx_length, fastq1)
    
    covfraction = round(float(maxcoverage / coverage), 3)
    if (coverage > maxcoverage):
        logging.debug('Downsampling to %s X coverage', maxcoverage)
        downcmd = "seqtk sample -s100 {fastq1} {covfraction} > {downpath1}".format(**locals())
        downcmd2 = "seqtk sample -s100 {fastq2} {covfraction} > {downpath2}".format(**locals())
        for command in [downcmd, downcmd2]:
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
    cmd = "ribo run -r {sra} -F {readsf} -R {readsr} --cores {cores} --threads {threads} -v 1 --serialize -o {output} --subassembler skesa --stages score".format(**locals())

    if readsr is None:
        cmd = "ribo run -r {sra} -F {readsf} --cores {cores} --threads {threads} -v 1 --serialize -o {output} --subassembler skesa --stages score".format(**locals())

    logging.debug('Running riboSeed: %s', cmd)
    return(cmd)

def extract_16s_from_contigs(input_contigs, barr_out, output):
    """Uses barrnap to identify rRNA operons within the riboSeed assembled contigs, then uses extractRegion to extract the 16S sequences """

    barrnap = "barrnap {input_contigs} > {barr_out}".format(**locals())
    logging.debug('Extracting 16S sequences: %s', barrnap)
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
                
                cmd = "extractRegion \'{results16s[3]}{results16s[0]} :{results16s[1]}:{results16s[2]}\' -f {input_contigs} -v 1 >> {output}".format(**locals())
                subprocess.run(cmd,
                               shell=sys.platform !="win32",
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               check=True)
                
    return(output)

def alignment(fasta, output):
    output = os.path.join(output, "alignment")
    os.makedirs(output)
    seqout = os.path.join(output, "16soneline.fasta")
    mafftout = os.path.join(output, "MSA.fasta")
    iqtreeout = os.path.join(output, "iqtree")
    seqcmd = "seqtk seq -S {fasta} > {seqout}".format(**locals())
    mafftcmd = "mafft {seqout} > {mafftout}".format(**locals())
    iqtreecmd = "iqtree -s {mafftout} -nt AUTO > iqtreeout".format(**locals())
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
    best_reference=os.path.join(pob_dir, "best_reference")
    with open(best_reference, "r") as infile:
        for line in infile:
            best_ref_fasta = line.split('\t')[0]
    logging.debug('Quality trimming reads')
    trimmed_fastq1, trimmed_fastq2 = run_sickle(fastq1=rawreadsf,
                                                fastq2=rawreadsr,
                                                output_dir=sickle_out) 
       
    logging.debug('Trimmed forward reads: %s', trimmed_fastq1)
    logging.debug('Trimmed reverse reads: %s', trimmed_fastq2)
    
    downsampledf, downsampledr = downsample(approx_length=args.approx_length,
                                            fastq1=trimmed_fastq1,
                                            fastq2=trimmed_fastq2, 
                                            maxcoverage=args.maxcov,
                                            destination=this_output)

    logging.debug('Downsampled forward reads: %s', downsampledf)
    logging.debug('Downsample reverse reads: %s', downsampledr)
       
    riboseed_cmd = run_riboseed(sra=best_ref_fasta, readsf=downsampledf,
                                readsr=downsampledr, cores=args.cores,
                                threads=1, output=ribo_dir)
    subprocess.run(riboseed_cmd,
                   shell=sys.platform !="win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
    
    dir_for_16s=os.path.join(this_output, "16s", "")
    barr_out=os.path.join(dir_for_16s, "barrnap")
    sixteens_extracted=os.path.join(dir_for_16s, "ribo16s")
    ribo_contigs = os.path.join(this_output, "riboSeed", "seed",
                                "final_long_reads", "riboSeedContigs.fasta")
    extract_16s_from_contigs(input_contigs=ribo_contigs, 
                             barr_out=barr_out, output=sixteens_extracted)
    alignment(fasta=sixteens_extracted, output=this_output)

    
def main():
    args=get_args()
    logging.basicConfig(
        level=logging.DEBUG,
        filename= "16db.log",
        format="%(asctime)s - %(levelname)s - %(message)s",
        )
    console = logging.StreamHandler()
    console.setLevel(logging.DEBUG)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)
    logger = logging.getLogger(__name__)

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
            logging.debug('Warning: genome directory exists but is ' +
                  'empty: downloading genomes')
            shutil.rmtree(args.genomes_dir)
            gng.main(args)
        else:
            pass
    else:
        gng.main(args)
    if filtered_sras == []:
        logging.debug('No complete genomes found on NCBI by sraFind')
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
            logging.debug('Downloading SRA: %s', accession)
            download_SRA(SRA=accession, destination=this_output)
        
            rawreadsf=os.path.join(this_output, "raw/", accession + "_1.fastq")
            rawreadsr=os.path.join(this_output, "raw/", accession + "_2.fastq")
            process_strain(rawreadsf, rawreadsr, this_output, args)


if __name__ == '__main__':
    main()
