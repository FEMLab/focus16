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
import glob

from plentyofbugs import get_n_genomes as gng
from py16db.run_sickle import run_sickle
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

#from Bio.Seq import Seq
from . import __version__


class bestreferenceError(Exception):
    pass


class downsamplingError(Exception):
    pass


class riboSeedError(Exception):
    pass


class extracting16sError(Exception):
    pass


def setup_logging(args):
    logging.basicConfig(
        level=logging.DEBUG,
        filename=os.path.join(args.output_dir, "16db.log"),
        format="%(asctime)s - %(levelname)s - %(message)s",
        )
    logger = logging.getLogger(__name__)
    console_err = logging.StreamHandler(sys.stderr)
    console_err.setLevel(logging.DEBUG)
    console_err_format = logging.Formatter(
        str("%(asctime)s \u001b[3%(levelname)s\033[1;0m  %(message)s"),
        "%H:%M:%S")
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
    parser.add_argument("-o", "--output_dir",
                        help="path to output", required=True)
    parser.add_argument("-n", "--organism_name",
                        help="genus or genus species in quotes",
                        required=True)
    parser.add_argument("--SRA_list",
                        help="path to file containing list of sras " +
                        "for assembly[one column]",
                        required=False)
    parser.add_argument("--version", action='version',
                        version='%(prog)s {version}'.format(
                            version=__version__))
    parser.add_argument("-l", "--approx_length",
                        help="Integer for approximate genome length",
                        required=False, type=int)
    parser.add_argument("-s", "--sraFind_path", dest="sra_path",
                        help="path to sraFind file",
                        default="sraFind-All-biosample-with-SRA-hits.txt",
                        required=False)
    parser.add_argument("--SRAs", default=None, nargs="+",
                        help="run pipeline on this (these) SRA(s) only",
                        required=False)
    parser.add_argument("-g", "--genomes_dir",
                        help="path to directory containing, or empty, " +
                        "candidate genomes for reference",
                        required=True)
    #  Note that this agument doesn't get called, but is inheirited by get_n_genomes
    parser.add_argument("-p", "--prokaryotes", action="store",
                        help="path_to_prokaryotes.txt",
                        default="./prokaryotes.txt",
                        required=False)
    parser.add_argument("-S", "--n_SRAs", help="max number of SRAs to be run",
                        type=int, required=False)
    parser.add_argument("-R", "--n_references",
                        help="max number of reference strains to consider",
                        type=int, required=False)
    parser.add_argument("--get_all",
                        help="if a biosample is associated with " +
                        "multiple libraries, default behaviour is to " +
                        "download the first only.  Use --get_all to " +
                        "analyse each library",
                        action="store_true", required=False)
    parser.add_argument("--cores",
                        help="integer for how many cores you wish to use",
                        default=1,
                        required=False, type=int)
    parser.add_argument("--maxcov",
                        help="integer for maximum coverage of reads",
                        default=50,
                        required=False, type=int)
    parser.add_argument("--example_reads",
                        help="input of example reads", nargs='+',
                        required=False, type=str)
    parser.add_argument("--subassembler",
                        help="which program should riboseed " +
                        "use for sub assemblies",
                        choices=["spades", "skesa"],
                        required=False, default="spades")
    # this is needed for plentyofbugs, should not be user set
    parser.add_argument("--nstrains", help=argparse.SUPPRESS,
                        type=int, required=False)
    parser.add_argument("--memory",
                        help="amount of RAM to be used",
                        default=4,
                        required=False, type=int)
    args = parser.parse_args()
    # plentyofbugs uses args.nstrains, but wecall it args.n_references for clarity
    args.nstrains = args.n_references
    if args.SRAs is None:
        if args.n_SRAs is None:
            print("if not running with --SRAs, " +
                  "then --n_SRAs must be provided!")
            sys.exit(1)
    return(args)


def check_programs(logger):
    """exits if the following programs are not installed"""

    required_programs = [
        "ribo", "barrnap", "fasterq-dump", "mash",
        "skesa", "plentyofbugs", "iqtree"]
    for program in required_programs:
        if shutil.which(program) is None:
            logger.error('%s is not installed: exiting.', program)
            sys.exit(1)


def filter_SRA(sraFind, organism_name, strains, get_all, logger):
    """sraFind [github.com/nickp60/srafind], contains"""
    results = []
    with open(sraFind, "r") as infile:
        for line in infile:
            split_line = [x.replace('"', '').replace("'", "") for x
                          in line.strip().split("\t")]
            if split_line[11].startswith(organism_name):
                if split_line[8].startswith("ILLUMINA"):
                    results.append(split_line[17])
    #  arbitrary seed number
    SEED = 8
    random.seed(SEED)
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
    """ takes a file list of  of SRAs, return list
    for if you wish to use SRAs that are very recent and ahead of sraFind
    """
    sras = []
    with open(list, "r") as infile:
        for sra in infile:
            sras.append(sra.strip())
    return sras


def download_SRA(cores, SRA, destination, logger):
    """Download SRAs, downsample forward reads to 1000000"""
    suboutput_dir_raw = os.path.join(destination, "")
    os.makedirs(suboutput_dir_raw)

    cmd = str("fasterq-dump {SRA} --threads {cores} -O " +
              "{suboutput_dir_raw} --split-files").format(**locals())

    try:
        subprocess.run(cmd,
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)
    except subprocess.CalledProcessError:
        logger.critical("Error with fasterq-dump")


def pob(genomes_dir, readsf, output_dir, logger):
    """use plentyofbugs to identify best reference
    Uses plentyofbugs, a package that useqs mash to
    find the best reference genome for draft genome
    """

    pobcmd = str("plentyofbugs -g {genomes_dir} -f {readsf} -o {output_dir} " +
                 "--downsampling_ammount 1000000").format(**locals())
    logger.debug('Finding best reference genome: %s', pobcmd)

    for command in [pobcmd]:
        try:
            subprocess.run(command,
                           shell=sys.platform != "win32",
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           check=True)
            best_ref = os.path.join(output_dir, "best_reference")
        except:
            raise bestreferenceError("Error running the following command: %s",
                                     command)

    with open(best_ref, "r") as infile:
        for line in infile:
            sraacc = line.strip().split('\t')
            sim = float(sraacc[1])
            ref = sraacc[0]
            percentSim = float(100.0 - sim)
            #  this isn't really true, but is sort of helpful
            logger.debug("Reference genome similarity: %s", percentSim)

            length_path = os.path.join(output_dir, "genome_length")
            cmd = "wc -c {ref} > {length_path}".format(**locals())
            subprocess.run(cmd,
                           shell=sys.platform != "win32",
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           check=True)
            return(ref, sim)


def download_SRA(cores, SRA, destination, logger):
    """Download SRAs, downsample forward reads to 1000000"""
    suboutput_dir_raw = os.path.join(destination, "")
    os.makedirs(suboutput_dir_raw)

    cmd = str("fasterq-dump {SRA} --threads {cores} -O " +
              "{suboutput_dir_raw} --split-files ").format(**locals())

    try:
        subprocess.run(cmd,
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)
    except subprocess.CalledProcessError:
        logger.critical("Error with fasterq-dump")


def check_rDNA_copy_number(ref, output, logger):
    """ensure reference has multiple rDNAs
    Using barrnap to check that there are multiple rDNA copies in the reference genome
    """
    os.makedirs(os.path.join(output, "barrnap_reference"), exist_ok=True)
    barroutput = os.path.join(output, "barrnap_reference",
                              os.path.basename(ref) + ".gff")
    cmd = "barrnap {ref} > {barroutput}".format(**locals())
    subprocess.run(cmd,
                   shell=sys.platform != "win32",
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
    return rrn_num
    # if rrn_num > 1:
    #     logger.debug('%s rrn operons detected in reference genome', rrn_num)
    # else:
    #     logger.critical('SRA does not have multiple rrn operons')
    #     sys.exit(1)
    # return(rrn_num)


def get_and_check_ave_read_len_from_fastq(fastq1, logger=None):
    """return average read length in fastq1 file from first N reads
    from LP: taken from github.com/nickp60/riboSeed/riboSeed/classes.py;
    """
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
    """Obtains the coverage for a read set given the estimated genome size"""

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


def downsample(read_length, approx_length, fastq1, fastq2,
               maxcoverage, destination, logger):
    """downsample for optimal assembly
    Given the coverage from coverage(), downsamples the reads if over
    the max coverage set by args.maxcov. Default 50.
    """

    suboutput_dir_downsampled = destination
    os.makedirs(suboutput_dir_downsampled)
    downpath1 = os.path.join(suboutput_dir_downsampled,
                             "downsampledreadsf.fastq")
    downpath2 = os.path.join(suboutput_dir_downsampled,
                             "downsampledreadsr.fastq")
    coverage = get_coverage(read_length, approx_length,
                            fastq1, fastq2, logger=logger)
    covfraction = round(float(maxcoverage / coverage), 3)
    downcmd = "seqtk sample -s100 {fastq1} {covfraction} > {downpath1}".format(**locals())
    downcmd2 = "seqtk sample -s100 {fastq2} {covfraction} > {downpath2}".format(**locals())
    # at least downsample the forward/single reads, but add the
    # other command if using paired reads
    commands = [downcmd]
    if (coverage > maxcoverage):
        logger.debug('Downsampling to %s X coverage', maxcoverage)
        if fastq2 is not None:
            commands.append(downcmd2)
        else:
            downpath2 = None
        for command in commands:
            try:
                # logger.debug(command)
                subprocess.run(command,
                               shell=sys.platform != "win32",
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               check=True)
            except:
                raise downsamplingError("Error running following command ", command)
        return(downpath1, downpath2)
    else:
        logger.debug('Skipping downsampling as max coverage is < %s', maxcoverage)
        return(fastq1, fastq2)


def make_riboseed_cmd(sra, readsf, readsr, cores, subassembler, threads,
                      output, memory, logger):
    """Runs riboSeed to reassemble reads """
    cmd = str("ribo run -r {sra} -F {readsf} -R {readsr} --cores {cores} " +
              "--threads {threads} -v 1 --serialize -o {output} " +
              "--subassembler {subassembler} --stages score " +
              "--memory {memory}").format(**locals())

    if readsr is None:
        cmd = str("ribo run -r {sra} -S1 {readsf} --cores {cores} " +
                  "--threads {threads} -v 1 --serialize -o {output} " +
                  "--subassembler {subassembler} --stages score " +
                  "--memory {memory}").format(**locals())
    logger.debug('Running riboSeed: %s', cmd)
    return(cmd)


def process_strain(rawreadsf, rawreadsr, this_output, args, logger):
    pob_dir = os.path.join(this_output, "plentyofbugs")

    ribo_dir = os.path.join(this_output, "riboSeed")
    sickle_out = os.path.join(this_output, "sickle")
    read_length = get_and_check_ave_read_len_from_fastq(fastq1=rawreadsf,
                                                        logger=logger)
    best_reference = os.path.join(pob_dir, "best_reference")
    if not os.path.exists(best_reference):
        pob(genomes_dir=args.genomes_dir, readsf=rawreadsf,
            output_dir=pob_dir, logger=logger)

    with open(best_reference, "r") as infile:
        for line in infile:
            best_ref_fasta=line.split('\t')[0]

    if args.approx_length is None:
        genome_length=os.path.join(pob_dir, "genome_length")
        with open(genome_length, "r") as infile:
            for line in infile:
                approx_length = float(line.split()[0])
                logger.debug("Using genome length: %s", approx_length)
    else:
        approx_length = args.approx_length

    logger.debug('Quality trimming reads')
    trimmed_fastq1, trimmed_fastq2 = run_sickle(fastq1=rawreadsf,
                                                fastq2=rawreadsr,
                                                output_dir=sickle_out)
    logger.debug('Quality trimmed f reads: %s', trimmed_fastq1)
    logger.debug('Quality trimmed r reads: %s', trimmed_fastq2)


    logger.debug('Downsampling reads')

    downsampledf, downsampledr = downsample(
        approx_length=approx_length,
        fastq1=trimmed_fastq1,
        fastq2=trimmed_fastq2,
        maxcoverage=args.maxcov,
        destination=os.path.join(this_output, "downsampled"),
        read_length=read_length,
        logger=logger)
    logger.debug('Downsampled f reads: %s', downsampledf)
    logger.debug('Downsampled r reads: %s', downsampledr)


    riboseed_cmd = make_riboseed_cmd(sra=best_ref_fasta, readsf=downsampledf,
                                     readsr=downsampledr, cores=args.cores,
                                     memory=args.memory,
                                     subassembler=args.subassembler,
                                     threads=1, output=ribo_dir, logger=logger)

    status = os.path.join(this_output, "status")

    # file that will contain riboseed contigs
    ribo_contigs = os.path.join(this_output, "riboSeed", "seed",
                                "final_long_reads", "riboSeedContigs.fasta")
    if "RIBOSEED COMPLETE" not in parse_status_file(status):
        try:
            subprocess.run(riboseed_cmd,
                           shell=sys.platform != "win32",
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           check=True)

            with open(status, "a") as statusfile:
                statusfile.write("RIBOSEED COMPLETE")

        except:
            if os.path.exists(ribo_contigs):
                pass
            else:
                raise riboSeedError("Error running the following command: %s",
                                    riboseed_cmd)
    else:
        logger.debug("Skipping riboSeed")


    if not os.path.exists(ribo_contigs):
        logger.critical(
            "riboSeed was not successful; for details, see log file at %s",
            os.path.join(this_output, "riboSeed", "run_riboSeed.log")
        )


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
    logger.debug("Running the following commands: ")

    for cmd in [seqcmd, mafftcmd, iqtreecmd]:
        logger.debug("   " + cmd)
        subprocess.run(cmd,
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)
    logger.debug('Performing alignment with mafft')
    logger.debug('Building tree with iqtree')
    return(output)


def parse_status_file(path):
    # because downloading and assembling can fail for many reasons,
    # we write out status to a file.  this allows for easier restarting of
    # incomplete runs
    if not os.path.exists(path):
        return []
    status = []
    with open(path, "r") as statusfile:
        for line in statusfile:
            status.append(line.strip())
    return(status)


def fetch_sraFind_data(dest_path):
    sraFind_results = "https://raw.githubusercontent.com/nickp60/sraFind/master/results/sraFind-All-biosample-with-SRA-hits.txt"
    # gets just the file name
    if not os.path.exists(dest_path):
        print("Downloading sraFind Dump")
        download_sraFind_cmd = str("wget " + sraFind_results + " -O " + dest_path)
        subprocess.run(
            download_sraFind_cmd,
            shell=sys.platform != "win32",
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True)

def extract_16s_from_assembly_list(all_assemblies, args, logger):

    extract16soutput = os.path.join(args.output_dir, "ribo16s.fasta")
    if os.path.exists(extract16soutput):
        os.remove(extract16soutput)
    results16s = {}  # [sra_#, chromosome, start, end, reverse complimented]
    nseqs = 0
    for assembly in all_assemblies:
        sra = assembly.split(os.path.sep)[1]
        barr_out = os.path.join(os.path.join(args.output_dir, sra, "barrnap"))
        with open(assembly, "r") as inf:
            barrnap = "barrnap {assembly} > {barr_out}".format(**locals())
            logger.debug('Identifying 16S sequences with barnap: %s', barrnap)
            try:
                subprocess.run(barrnap,
                               shell=sys.platform != "win32",
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               check=True)
            except:
                raise extracting16sError("Error running the following command %s", barrnap)
            with open(barr_out, "r") as rrn, open(assembly, "r")  as asmb, open(extract16soutput, "a") as outf:
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
                        chrom = line[0]
                        ori = line[6]
                        start = int(line[3])
                        end = int(line[4])
                        thisid = "{}_{}".format(sra, rrn_num)
                        results16s[thisid] = [chrom, start, end, line[6]]
                        for rec in SeqIO.parse(asmb, "fasta"):
                            if rec.id  == chrom:
                                seq = rec.seq[start + 1: end + 1]
                                if  ori == "-":
                                    seq = seq.reverse_complement()
                                thisdesc = "{chrom}:{start}:{end}({ori})".format(**locals())
                                SeqIO.write(SeqRecord(seq, id=thisid, description=thisdesc), outf,  "fasta")
                                nseqs += 1
    return(nseqs, extract16soutput)


def write_pass_fail(args, stage, status, note):
    """
    format fail messages in tabular fomat:
    organism\tstage\tmessage
    """
    path = os.path.join(args.output_dir, "SUMMARY")
    org = args.organism_name
    with open(path, "a") as failfile:
        failfile.write("{}\t{}\t{}\t{}\n".format(org, status, stage, note))


def main():
    args = get_args()

    genomesdir = os.path.join(args.genomes_dir, "")

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    if os.path.exists(os.path.join(args.output_dir, "SUMMARY")):
        os.remove(os.path.join(args.output_dir, "SUMMARY"))

    logger = setup_logging(args)

    check_programs(logger)

    fetch_sraFind_data(args.sra_path)

    if args.SRAs is not None:

        filtered_sras = args.SRAs
    elif args.SRA_list is not None:
        filtered_sras = sralist(list=args.SRA_list)
    else:
        filtered_sras = filter_SRA(
            sraFind=args.sra_path,
            organism_name=args.organism_name,
            strains=args.n_SRAs,
            logger=logger,
            get_all=args.get_all)

    sra_num = len(filtered_sras)
    logger.debug("%s SRAs found:", sra_num)

    # check basic integrity of genomes dir
    # it might exist, but making it here simplifies the control flow later
    # it seems counter intuitive, but checking if the dir we might have just
    # created is easier than checking if it exists/is intact twice
    os.makedirs(args.genomes_dir, exist_ok=True)
    try:
        if len(glob.glob(os.path.join(args.genomes_dir, "*.fna"))) == 0:
            logger.debug('Warning: genome directory exists but is ' +
                         'empty: downloading genomes')
            shutil.rmtree(args.genomes_dir)
            gng.main(args)  # remake and fill with genomes
        if len(glob.glob(os.path.join(args.genomes_dir, "*.fna.gz"))) != 0:
            logger.debug('Warning: genome downloading may have been interupted; ' +
                         'downloading fresh')
            shutil.rmtree(args.genomes_dir)
            gng.main(args)  # remake and fill with genomes
    except Exception as e:
        logger.error("Error downloading genomes")
        logger.error(e)
        write_pass_fail(args, status="FAIL",
                        stage="global", note="Downloading_references")
        sys.exit(1)

    logger.debug("checking reference genomes for rDNA counts")
    for potential_reference in glob.glob(os.path.join(args.genomes_dir, "*.fna")):
        rDNAs = check_rDNA_copy_number(ref=potential_reference,
                                       output=args.genomes_dir,
                                       logger=logger)
        if rDNAs < 2:
            logger.warning(
                "reference %s does not have multiple rDNAs; excluding", potential_reference)
            os.remove(potential_reference)
    if len(glob.glob(os.path.join(args.genomes_dir, "*.fna"))) == 0:
        logger.error("No usable reference genome found!")
        write_pass_fail(args, status="FAIL",
                        stage="global",
                        note="filtering_references: no references had more than 1 rDNA")

    if filtered_sras == []:
        if args.example_reads is None:
            logger.error('No SRAs found on NCBI by sraFind')
            write_pass_fail(args, status="FAIL",
                            stage="global",
                            note="identifying SRAs: double check SRAfind sheet to ensure data exists for your organism")

    if args.example_reads is not None:
        this_output = os.path.join(args.output_dir, "example", "")
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

            if "SRA COMPLETE" not in parse_status_file(path=status):
                logger.debug('Downloading SRA: %s', accession)
                # a fresh start
                if os.path.exists(this_data):
                    shutil.rmtree(this_data)

                download_SRA(cores=args.cores, SRA=accession,
                             destination=this_data, logger=logger)
                with open(status, "a") as statusfile:
                    statusfile.write("SRA COMPLETE\n")
            else:
                logger.debug("Skipping SRA download: %s", accession)

            rawreadsf = os.path.join(this_data, accession + "_1.fastq")
            rawreadsr = os.path.join(this_data, accession + "_2.fastq")
            if not os.path.exists(rawreadsr):
                rawreadsr = None
            if not os.path.exists(rawreadsf):
                logger.info('Forward reads not detected')
                continue
            try:
                if "PROCESSED" not in parse_status_file(path=status):
                    if os.path.exists(this_results):
                        shutil.rmtree(this_results)
                    process_strain(rawreadsf, rawreadsr, this_results,
                                   args, logger)
                    with open(status, "a") as statusfile:
                        statusfile.write("PROCESSED\n")

                else:
                    logger.debug("Already processed: %s", accession)
                write_pass_fail(args, status="PASS",
                                stage=accession,
                                note="Processed strain")

            except subprocess.CalledProcessError:
                write_pass_fail(args, status="FAIL",
                                stage=accession,
                                note="unknown failure")
                logger.error('Unknown subprocess error')
                continue
            except bestreferenceError as e:
                write_pass_fail(args, status="FAIL",
                                stage=accession,
                                note="unknown error selecting reference")
                logger.error(e)
                continue
            except downsamplingError as e:
                write_pass_fail(args, status="FAIL",
                                stage=accession,
                                note="unknown error downsampling")
                logger.error(e)
                continue
            except riboSeedError as e:
                write_pass_fail(args, status="FAIL",
                                stage=accession,
                                note="unknown failure running riboSeed")
                logger.error(e)
                continue
            except extracting16sError as e:
                write_pass_fail(args, status="FAIL",
                                stage=accession,
                                note="unknown extracting 16s sequences")
                logger.error(e)
                continue

    alignoutput = os.path.join(args.output_dir, "allsequences",  "")
    pathtotree = os.path.join(alignoutput, "MSA.fasta.tree")
    all_assemblies  =  glob.glob(
        os.path.join(args.output_dir, "*", "results", "riboSeed",
                     "seed", "final_long_reads", "riboSeedContigs.fasta"))
    n_extracted_seqs, extract16soutput  =  extract_16s_from_assembly_list(
        all_assemblies, args, logger)
    logger.debug("Wrote out  %i  sequences", n_extracted_seqs)
    if n_extracted_seqs == 0:
        write_pass_fail(args, status="FAIL",
                        stage="global",
                        note="no 16s sequences detected in re-assemblies")
        logger.warning("No  16s Sequences  recovered. exiting")

        sys.exit()


    if os.path.exists(alignoutput):
        shutil.rmtree(alignoutput)
    #alignment(fasta=extract16soutput, output=alignoutput, logger=logger)
    logger.debug('Maximum-likelihood tree available at: %s', pathtotree)
    write_pass_fail(args, status="PASS",
                    stage="global",
                    note="")


if __name__ == '__main__':
    main()
