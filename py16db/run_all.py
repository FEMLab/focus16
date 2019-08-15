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
from . import __version__

def setup_logging(args):
    logging.basicConfig(
        level=logging.DEBUG,
        filename= os.path.join(args.output_dir, "16db.log"),
        format="%(asctime)s - %(levelname)s - %(message)s",
        )
    logger = logging.getLogger(__name__)
    console_err = logging.StreamHandler(sys.stderr)
    console_err.setLevel(logging.DEBUG)
    console_err_format = logging.Formatter(
        str("%(asctime)s \u001b[3%(levelname)s\033[1;0m  %(message)s"), "%H:%M:%S")
    console_err.setFormatter(console_err_format)
   
    logging.addLevelName(logging.DEBUG,    "4m --")
    logging.addLevelName(logging.INFO,     "2m ==")
    logging.addLevelName(logging.WARNING,  "3m !!")
    logging.addLevelName(logging.ERROR,    "1m xx")
    logging.addLevelName(logging.CRITICAL, "1m XX")
    logger.addHandler(console_err)
    return logger

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output_dir", help="path to output", required=True)
    parser.add_argument("-n", "--organism_name", help="genus or genus species in quotes",
                        required=True)
    parser.add_argument("--sra_list", help="Input a list of sras for assembly[one column]",
                        required=False)
    parser.add_argument("--version", action='version',
                        version='%(prog)s {version}'.format(version=__version__))
    parser.add_argument("-l", "--approx_length", help="Integer forapproximate genome length",
                        required=True, type=int)
    parser.add_argument("-s", "--sraFind_path", dest="sra_path", 
                        help="path to sraFind file", required=False)
    parser.add_argument("--single_SRA", default=None,
                        help="run pipeline on this SRA accession only",
                        required=False)
    parser.add_argument( "-g", "--genomes_dir",help="path to directory containing, or empty, candidate genomes for reference",
                         required=True)
    parser.add_argument("-p", "--prokaryotes", action="store",
                        help="path_to_prokaryotes.txt", default="./prokaryotes.txt",
                        required=False)
    parser.add_argument("-S", "--nstrains", help="number of SRAs to be run",             
                        type=int, required=True)
    parser.add_argument("--get_all", help="get both SRAs if organism has two",
                        action="store_true", required=False)
    parser.add_argument("--cores", help="integer for how many cores you wish to use", default=1,
                        required=False, type=int)
    parser.add_argument("--maxcov", help="integer for maximum coverage of reads", default=50,
                        required=False, type=int)
    parser.add_argument("--example_reads", help="input of example reads", nargs='+',
                        required=False, type=str)
    parser.add_argument("--subassembler", 
                        help="which program should riboseed use for sub assemblies",
                        choices=["spades", "skesa"], 
                        required=False, default="spades")
    parser.add_argument("--memory", 
                        help="amount of RAM to be used",
                        default=4, 
                        required=False, type=int)
    
    return(parser.parse_args())


def check_programs(logger):
    """exits if the following programs are not installed"""
    
    required_programs = ["ribo", "barrnap",
                         "fasterq-dump", "mash", "skesa", "plentyofbugs", "iqtree"]
    for program in required_programs:
        if shutil.which(program) is None:
            logger.error ( '%s is not installed: exiting.', program) 
            sys.exit(1)
   
           
def filter_SRA(sraFind, organism_name, strains, get_all, logger):
    """sraFind [github.com/nickp60/srafind], contains"""
    results = []
    with open(sraFind, "r") as infile:
        for line in infile:
            split_line = [x.replace('"', '').replace("'", "") for x in line.strip().split("\t")]
            if split_line[11].startswith(organism_name):
                if split_line[8].startswith("ILLUMINA"):
                    results.append(split_line[17])
    random.shuffle(results)
    
    if strains != 0:
        results = results[0:strains]
        logger.debug('Found SRAs: %s', results)

    sras = []
    for result in results:
        these_sras = result.split(",")
        if get_all:
            for sra in these_sras:
                sras.append(sra)
        else:
            sras.append(these_sras[0])

    return(sras)

def sralist(list):
    """ takes a list of SRAs as input, for if you wish to use SRAs that are very recent and ahead of sraFind
    """
    sralistn = []
    with open(list, "r") as infile:
        for sra in infile:
            sra.split()
            sralistn.append(sra)
            sralist = [sra.strip() for sra in sralistn]
    return(sralist)
    

def download_SRA(cores, SRA, destination, logger):
    """Download SRAs, downsample forward reads to 1000000"""
    suboutput_dir_raw = os.path.join(destination, "")
    os.makedirs(suboutput_dir_raw)

    cmd = "fasterq-dump {SRA} --threads {cores} -O {suboutput_dir_raw} --split-files ".format(**locals())
  
    try:
        subprocess.run(cmd,
                       shell=sys.platform !="win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)
    except subprocess.CalledProcessError():
        logger.critical("Error with fasterq-dump")
      
def pob(genomes_dir, readsf, output_dir, logger):
    """Uses plentyofbugs, a package that useqs mash to find the best reference genome for draft genome """
       
    pobcmd = "plentyofbugs -g {genomes_dir} -f {readsf} -o {output_dir} --downsampling_ammount 1000000".format(**locals())
    logger.debug('Finding best reference genome: %s', pobcmd)
    
    for command in [pobcmd]:
        try:
            subprocess.run(command,
                           shell=sys.platform !="win32",
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           check=True)
            best_ref = os.path.join(output_dir, "best_reference")
        except:
            raise bestreferenceError("Error running the following command: %s", command)
    
    
    with open(best_ref, "r") as infile:
        for line in infile:
            sim = float(line.split('\t')[1])
            percentSim = float(100.0 - sim)
            logger.debug("Reference genome similarity: {percentSim}%".format(**locals()))
            sraacc = line.strip().split('\t')            
            
            return(sraacc)

def check_rDNA_copy_number(ref, output, logger):
    """ Using barrnap to check that there are multiple rDNA copies in the reference genome   """
    barroutput = os.path.join(output, "barrnap_reference")
    cmd = "barrnap {ref} > {barroutput}".format(**locals())
    subprocess.run(cmd,
                   shell=sys.platform !="win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
    rrn_num = 0
    with open(barroutput, "r") as rrn:
        for rawline in rrn:
            line = rawline.strip().split('\t')
            if line[0].startswith("##"):
                continue
            if line[8].startswith("Name=16S"):
                rrn_num += 1
    if rrn_num > 1:
        logger.debug('%s rrn operons detected in reference genome', rrn_num)
    else:
        logger.critical('SRA does not have multiple rrn operons')
        sys.exit(1)
    return(rrn_num)            


def get_and_check_ave_read_len_from_fastq(fastq1, logger=None):
    """from LP: taken from github.com/nickp60/riboSeed/riboSeed/classes.py; return average read length in fastq1 file from first N reads"""
    count, tot = 0, 0
    if os.path.splitext(fastq1)[-1] in ['.gz', '.gzip']:
        open_fun = gzip.open
    else:
        open_fun = open


    with open_fun(fastq1, "rt") as file_handle:
        data = SeqIO.parse(file_handle, "fastq")
        logger.debug("Obtaining average read length from first 30 reads")
        for read in data:
            count += 1
            tot += len(read)
            if count == 30:
                break

    
    ave_read_len = float(tot / 30)
    if ave_read_len < 65:
        logger.critical("Average read length is too short: %s", ave_read_len)
        sys.exit(1)
    if ave_read_len > 300:
        logger.critical("Average read length is too long: %s", ave_read_len)
        sys.exit(1)
    logger.debug("Average read length: %s", ave_read_len)
    return(ave_read_len)

def get_coverage(read_length, approx_length, fastq1, fastq2, logger):
    """Obtains the coverage for a read set, when given the estimated genome size"""
        
    if os.path.splitext(fastq1)[-1] in ['.gz', '.gzip']:
        open_fun = gzip.open
    else:
        open_fun = open
    logger.debug("Counting reads")
    
    with open_fun(fastq1, "rt") as data:
        for count, line in enumerate(data):
            pass
    
    if fastq2 is not None:
        read_length = read_length * 2
       
    coverage = float((count * read_length) / (approx_length * 4))
    logger.debug('Read coverage is : %s', coverage)
    return(coverage)

def downsample(read_length, approx_length, fastq1, fastq2, maxcoverage, destination, logger):
    """Given the coverage from coverage(), downsamples the reads if over the max coverage set by args.maxcov. Default 50."""

    suboutput_dir_downsampled = destination
    os.makedirs(suboutput_dir_downsampled)
    downpath1 = os.path.join(suboutput_dir_downsampled, "downsampledreadsf.fastq") 
    downpath2 = os.path.join(suboutput_dir_downsampled, "downsampledreadsr.fastq")
    coverage = get_coverage(read_length, approx_length, fastq1, fastq2, logger=logger)
    covfraction = round(float(maxcoverage / coverage), 3)
    downcmd = "seqtk sample -s100 {fastq1} {covfraction} > {downpath1}".format(**locals())
    downcmd2 = "seqtk sample -s100 {fastq2} {covfraction} > {downpath2}".format(**locals())     
    # at least downsample the forward/single reads, but add the other command if
    # using paired reads
    commands = [downcmd]
    if (coverage > maxcoverage):
        logger.debug('Downsampling to %s X coverage', maxcoverage)
        if fastq2 is not None:
            commands.append(downcmd2)
        else:
            downpath2 = None
        for command in commands:
           try:
               subprocess.run(command,
                              shell=sys.platform !="win32",
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
           except:
               raise downsamplingError("Error running following command ", command) 
        return(downpath1, downpath2)
    else:
        logger.debug('Skipping downsampling as max coverage is < %s', maxcoverage)
        return(fastq1, fastq2)
            
    
def run_riboseed(sra, readsf, readsr, cores, subassembler, threads, output, memory, logger):
    """Runs riboSeed to reassemble reads """
    cmd = "ribo run -r {sra} -F {readsf} -R {readsr} --cores {cores} --threads {threads} -v 1 --serialize -o {output} --subassembler {subassembler} --stages score --memory {memory}".format(**locals())

    if readsr is None:
        cmd = "ribo run -r {sra} -S1 {readsf} --cores {cores} --threads {threads} -v 1 --serialize -o {output} --subassembler {subassembler} --stages score --memory {memory}".format(**locals())

    logger.debug('Running riboSeed: %s', cmd)
    return(cmd)

def extract_16s_from_contigs(input_contigs, barr_out, output, logger):
    """Uses barrnap to identify rRNA operons within the riboSeed assembled contigs, then uses extractRegion to extract the 16S sequences """

    barrnap = "barrnap {input_contigs} > {barr_out}".format(**locals())
    logger.debug('Extracting 16S sequences: %s', barrnap)
    try:
        subprocess.run(barrnap,
                       shell=sys.platform !="win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)
    except:
        raise extracting16sError("Error running the following command %s", barrnap)

    results16s = []   # [chromosome, start, end, reverse complimented]
    with open(barr_out, "r") as rrn:
        rrn_num = 0
        for rawline in rrn:
            line = rawline.strip().split('\t')
            if line[0].startswith("##"):
               continue
            if line[8].startswith("Name=16S"):
                rrn_num += 1
                if line[6] == "-":
                    suffix = 'chromosome-RC@'
                else:
                    suffix = ''
                chrom=line[0]
                start=line[3]
                end=line[4]
                results16s = [chrom, start, end, suffix]
        

                cmd = "extractRegion \'{results16s[3]}{results16s[0]} :{results16s[1]}:{results16s[2]}\' -f {input_contigs} -v 1 >> {output}".format(**locals())
                try:
                    subprocess.run(cmd,
                                   shell=sys.platform !="win32",
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE,
                                   check=True)
                except:
                    raise extracting16sError("Error running following command ", cmd)
    if rrn_num == 1:
        logger.debug('%s 16S sequence extracted from genome', rrn_num)
    elif rrn_num > 1:
        logger.debug('%s 16S sequences extracted from genome', rrn_num)
    else:
        logger.critical('NO 16S sequences extracted', rrn_num)
 
    return(output)
            

def process_strain(rawreadsf, rawreadsr, this_output, args, logger):
    pob_dir=os.path.join(this_output, "plentyofbugs")
    ribo_dir=os.path.join(this_output, "riboSeed")
    sickle_out=os.path.join(this_output, "sickle")
    read_length = get_and_check_ave_read_len_from_fastq(fastq1=rawreadsf, logger=logger)
    
    pob(genomes_dir=args.genomes_dir, readsf=rawreadsf, output_dir=pob_dir, logger=logger)
      
    best_reference=os.path.join(pob_dir, "best_reference")
    with open(best_reference, "r") as infile:
        for line in infile:
            best_ref_fasta = line.split('\t')[0]
            check_rDNA_copy_number(ref=best_ref_fasta, output=this_output, logger=logger) 
            
    logger.debug('Quality trimming reads')
    trimmed_fastq1, trimmed_fastq2 = run_sickle(fastq1=rawreadsf,
                                                fastq2=rawreadsr,
                                                output_dir=sickle_out)
    logger.debug('Quality trimmed f reads: %s', trimmed_fastq1)
    logger.debug('Quality trimmed r reads: %s', trimmed_fastq2)
    
    
    logger.debug('Downsampling reads')
    downsampledf, downsampledr = downsample(approx_length=args.approx_length,
                                            fastq1=trimmed_fastq1,
                                            fastq2=trimmed_fastq2, 
                                            maxcoverage=args.maxcov,
                                            destination=os.path.join(this_output, "downsampled"),
                                            read_length=read_length,
                                            logger=logger)
    logger.debug('Downsampled f reads: %s', downsampledf)    
    logger.debug('Downsampled r reads: %s', downsampledr)    
    
    
    riboseed_cmd = run_riboseed(sra=best_ref_fasta, readsf=downsampledf,
                                readsr=downsampledr, cores=args.cores,
                                memory=args.memory,
                                subassembler=args.subassembler,
                                threads=1, output=ribo_dir, logger=logger)
    
    status = os.path.join(this_output, "status")

    #file that will contain riboseed contigs
    ribo_contigs = os.path.join(this_output, "riboSeed", "seed",
                                "final_long_reads", "riboSeedContigs.fasta")
    if not "RIBOSEED COMPLETE" in parse_status_file(status):            
        try:
            subprocess.run(riboseed_cmd,
                           shell=sys.platform !="win32",
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           check=True)
                
            with open(status, "a") as statusfile:
                statusfile.write("RIBOSEED COMPLETE")
                
        except:
            if os.path.exists(ribo_contigs):
                pass
            else:
                raise riboSeedError("Error running the following command: %s", riboseed_cmd)
    else:
        logger.debug("Skipping riboSeed")
        
    
    barr_out=os.path.join(this_output, "barrnap")
    extract16soutput=os.path.join(args.output_dir, "ribo16s")
    ribo_contigs = os.path.join(this_output, "riboSeed", "seed",
                                "final_long_reads", "riboSeedContigs.fasta")
    
    if not os.path.exists(ribo_contigs):
        logger.critical("Assembled contigs are not long enough")
    else:
        extract_16s_from_contigs(input_contigs=ribo_contigs, 
                                 barr_out=barr_out, output=extract16soutput, logger=logger)
    

def alignment(fasta, output, logger):
    """ Performs multiple sequence alignment with mafft and contructs a tree with iqtree
    """
    output = os.path.join(output, "alignment")
    os.makedirs(output)
    seqout = os.path.join(output, "16soneline.fasta")
    mafftout = os.path.join(output, "MSA.fasta")
    iqtreeout = os.path.join(output, "iqtree")
    seqcmd = "seqtk seq -S {fasta} > {seqout}".format(**locals())
    mafftcmd = "mafft {seqout} > {mafftout}".format(**locals())
    iqtreecmd = "iqtree -s {mafftout} -nt AUTO > {iqtreeout}".format(**locals())
    for cmd in [seqcmd, mafftcmd, iqtreecmd]:
        logger.debug('Performing alignment: %s', cmd)
        subprocess.run(cmd,
                       shell=sys.platform !="win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)
    return(output)


class bestreferenceError(Exception):
    pass
class downsamplingError(Exception):
    pass
class riboSeedError(Exception):
    pass
class extracting16sError(Exception):
    pass



def parse_status_file(path):
    if not os.path.exists(path):
        return []
    status = [] 
    with open(path, "r") as statusfile:
        for line in statusfile:
            status.append(line.strip())
    return(status)


def main():
    args=get_args()
    
    genomesdir = args.genomes_dir

    if not args.genomes_dir.endswith("/"):
        args.genomes_dir = args.genomes_dir + "/"
        
    if os.path.exists(args.output_dir):
        pass
    else:
        os.makedirs(args.output_dir)

    logger = setup_logging(args)
    
        
    check_programs(logger)


    if args.sra_path is None:
        args.sra_path  = os.path.join(args.output_dir, "sraFind", "sraFind-All-biosample-with-SRA-hits.txt")

    if not os.path.exists(os.path.dirname(args.sra_path)):
        args.sra_path = fsd.main(args, output_dir=os.path.dirname(args.sra_path))

    if args.single_SRA is not None:
        filtered_SRA = [args.single_SRA]
        
    elif args.sra_list is not None:
        filtered_sras = sralist(list=args.sra_list)
    else:
        filtered_sras = filter_SRA(
            sraFind=args.sra_path,
            organism_name=args.organism_name,
            strains=args.nstrains,
            logger=logger,
            get_all=args.get_all)

    if os.path.exists(args.genomes_dir):
        if len(os.listdir(args.genomes_dir)) == 0:
            logger.debug('Warning: genome directory exists but is ' +
                         'empty: downloading genomes')
            shutil.rmtree(args.genomes_dir)
            gng.main(args, logger)
        else:
            pass
    else:
        gng.main(args, logger)
    if filtered_sras == []:
        logger.debug('No SRAs found on NCBI by sraFind')
        
    
    if args.example_reads is not None:
        this_output=os.path.join(args.output_dir, "example", "")
        os.makedirs(this_output)        
        rawreadsf = args.example_reads[0]
        try:
            rawreadsr = args.example_reads[1]
        except IndexError:
            rawreadsr = None    
        process_strain(rawreadsf, rawreadsr, this_output, args, logger)
       
    else:
        for i, accession in enumerate(filtered_sras):
            this_output = os.path.join(args.output_dir, accession)
            this_data = os.path.join(this_output, "data")
            this_results = os.path.join(this_output, "results")
            os.makedirs(this_output, exist_ok=True)
            status = os.path.join(this_output, "status")
            logger.debug("Organism: %s", args.organism_name)
            # check status file for SRA COMPLETE
            
            if "SRA COMPLETE" not in parse_status_file(path=status) :
                logger.debug('Downloading SRA: %s', accession)
                # a fresh start
                if os.path.exists(this_data):
                    shutil.rmtree(this_data)

                download_SRA(cores=args.cores, SRA=accession, destination=this_data,logger=logger)
                with open(status, "a") as statusfile:
                    statusfile.write("SRA COMPLETE\n")
            else:
                logger.debug("Skipping SRA download: %s", accession)

            rawreadsf=os.path.join(this_data, accession + "_1.fastq")
            rawreadsr=os.path.join(this_data, accession + "_2.fastq")
            if not os.path.exists(rawreadsr):
                rawreadsr = None
            if not os.path.exists(rawreadsf):
                logger.critical('Forward reads not detected')
                continue
            try:
                if "PROCESSED" not in parse_status_file(path=status):
                    if os.path.exists(this_results):
                        shutil.rmtree(this_results)
                    process_strain(rawreadsf, rawreadsr, this_results, args, logger)
                    with open(status, "a") as statusfile:
                        statusfile.write("PROCESSED\n")
                            
                else:
                    logger.debug("Already processed: %s", accession)

            except subprocess.CalledProcessError:
                logger.error('Unknown subprocess error')
                continue
            except bestreferenceError as e:
                logger.error(e)
                continue
            except downsamplingError as e:
                logger.error(e)
                continue
            except riboSeedError as e:
                logger.error(e)
                continue
            except extracting16sError as e:
                logger.error(e)
                continue

    extract16soutput = os.path.join(args.output_dir, "ribo16s")
    alignoutput = os.path.join(args.output_dir, "allsequences")
    pathtotree = os.path.join(alignoutput, "MSA.fasta.tree")

    if os.path.exists(alignoutput):
        os.remove(alignoutput)
        alignment(fasta=extract16soutput, output=alignoutput, logger=logger)        
        logger.debug('Maximum-likelihood tree available at: %s', pathtotree) 
        

if __name__ == '__main__':
    main()
